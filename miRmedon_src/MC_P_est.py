import os
import pickle
import numpy as np
from scipy.stats import beta
import pysam


def phred_to_prob(asciiChar):
    """This function takes an ascii character and as an input and returns 
    the error probability"""
    return 10**(-(ord(asciiChar)-33)/10)


def get_error_probs(final_filt_bam, e_miRbase_dict):
    bam_input = pysam.AlignmentFile(final_filt_bam, 'rb')
    error_probs_dict = {}
    for alignment in bam_input.fetch(until_eof=True):
        ref_name = alignment.reference_name
        if ref_name.find('_e') != -1:
            basename = ref_name[:ref_name.find('_')]
            error_probs = list(map(phred_to_prob, list(alignment.qual)))
            src_editing_sites = e_miRbase_dict[basename][ref_name]['editingSites']
            ref_seq = alignment.get_reference_sequence()
            src_ref_seq = e_miRbase_dict[basename][ref_name]['sequence']
            match_index = src_ref_seq.find(ref_seq.upper())
            edited_sites = list(filter(lambda x: x >= 0 and x < len(ref_seq), [i - match_index for i in src_editing_sites]))
            if basename not in error_probs_dict.keys():
                error_probs_dict[basename] = {}
                for query_editing_site in edited_sites:
                    ref_editing_site = query_editing_site + match_index
                    error_probs_dict[basename][ref_editing_site] = [error_probs[query_editing_site]]
            elif basename in error_probs_dict.keys():
                for query_editing_site in edited_sites:
                    ref_editing_site = query_editing_site + match_index
                    if ref_editing_site not in error_probs_dict[basename].keys():
                        error_probs_dict[basename][ref_editing_site] = [error_probs[query_editing_site]]
                    elif ref_editing_site in error_probs_dict[basename].keys():
                        error_probs_dict[basename][ref_editing_site] += [error_probs[query_editing_site]]
    return error_probs_dict


def get_haplotypes_counts():
    haplotypes_counts = {}
    with open('counts.txt', 'r') as file:
        for line in file:
            line = line.rstrip('\n').split('\t')
            haplotypes_counts[line[0]] = float(line[1])
    return haplotypes_counts


def get_editing_data(final_filt_bam, path_to_emiRbase):
    e_miRbase_dict = pickle.load(open(path_to_emiRbase, 'rb'))
    haplotypes_counts = get_haplotypes_counts()
    error_probs_dict = get_error_probs(final_filt_bam=final_filt_bam, e_miRbase_dict=e_miRbase_dict)
    unedited_sites_counts = {}
    edited_sites_counts = {}
    for haplotype, count in haplotypes_counts.items():
        basename = haplotype[:haplotype.find('_')]
        sequence = e_miRbase_dict[basename][haplotype]['sequence']
        unedited_sites = [index for index, base in enumerate(sequence) if base == 'A']
        if basename not in unedited_sites_counts.keys():
            unedited_sites_counts[basename] = {}
            for i in unedited_sites:
                unedited_sites_counts[basename][i] = count
        elif basename in unedited_sites_counts.keys():
            for i in unedited_sites:
                if i not in unedited_sites_counts[basename].keys():
                    unedited_sites_counts[basename][i] = count
                elif i in unedited_sites_counts[basename].keys():
                    unedited_sites_counts[basename][i] += count
        if haplotype.find('_e') != -1:
            edited_sites = e_miRbase_dict[basename][haplotype]['editingSites']
            if basename not in edited_sites_counts.keys():
                edited_sites_counts[basename] = {}
                for i in edited_sites:
                    edited_sites_counts[basename][i] = count
            elif basename in edited_sites_counts.keys():
                for i in edited_sites:
                    if i not in edited_sites_counts[basename].keys():
                        edited_sites_counts[basename][i] = count
                    elif i in unedited_sites_counts[basename].keys():
                        edited_sites_counts[basename][i] += count
    return edited_sites_counts, unedited_sites_counts, haplotypes_counts, error_probs_dict


def monte_catlo_p_estimation(final_filt_bam, path_to_emiRbase, resamples):
    edited_sites_counts, unedited_sites_counts, haplotypes_counts, error_probs_dict = get_editing_data(final_filt_bam=final_filt_bam,
                                                                                                       path_to_emiRbase=path_to_emiRbase)
    editing_data_file = open('editing_info.txt', 'w')
    editing_data_file.write(
            'miRNA' + '\t' + 'position' + '\t' + 'edited' + '\t' + 'unedited' + '\t' + 
            'editing_level' + '\t' + 'LCI' + '\t' 'UCI' + '\t' + 'p_value' + '\n')
    edited_miRs = list(set(edited_sites_counts.keys()))
    editing_rates_dict = {}
    p_values_dict = {}
    for miRNA in edited_miRs:
        p_values_dict[miRNA] = {}
        editing_rates_dict[miRNA] = {}
        for position in list(edited_sites_counts[miRNA].keys()):
            edited_count = edited_sites_counts[miRNA][position]
            try:
                unedited_count = unedited_sites_counts[miRNA][position]
            except KeyError:
                unedited_count = 0
            beta_dist = beta(edited_count+1, unedited_count+1)
            LCI, UCI = beta_dist.ppf(0.025), beta_dist.ppf(0.975) 
            bootstrap_beta = beta.rvs(edited_count+1, unedited_count+1, size=resamples)
            error_probs = np.array(error_probs_dict[miRNA][position])
            ext_error_probs = np.random.choice(error_probs, resamples)
            arraySubtraction = bootstrap_beta - ext_error_probs
            positive = len(arraySubtraction[arraySubtraction <= 0])
            p_value = (positive+1)/(resamples+1)
            p_values_dict[miRNA][position] = p_value
            editing_rate = beta_dist.mean()
            editing_rates_dict[miRNA][position] = editing_rate
            editing_data_file.write(
                    miRNA + '\t' + str(position+1) + '\t' + str(edited_count) + '\t' + str(unedited_count) + '\t' +
                    str(editing_rate) + '\t' + str(LCI) + '\t' + str(UCI) + '\t' + str(p_values_dict[miRNA][position]) + '\n')
