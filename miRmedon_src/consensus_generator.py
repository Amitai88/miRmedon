from functools import reduce
from decimal import Decimal
import numpy as np
import random
import subprocess
import os
import pysam

################################################### Auxilary functions #############################################


def phred_to_prob(asciiChar):
    """This function takes an ascii character and as an input and returns 
    the error probability"""
    return 10**(-(ord(asciiChar)-33)/10)


def length_correction(sequences, phred_seqs):
    corrected_phred_seqs = []
    for i in range(len(sequences)):
        startGap = ''
        endGap = ''
        for j in range(len(sequences[i])):
            if sequences[i][j] not in ['A', 'T', 'C', 'G', 'N']:
                startGap += '~'
            else:
                break
        for j in range(len(sequences[i]))[::-1]:
            if sequences[i][j] not in ['A', 'T', 'C', 'G', 'N']:
                endGap += '~'
            else:
                break
        corrected_phred_seqs.append(startGap + phred_seqs[i] + endGap) 
    return corrected_phred_seqs


def consensus_calculator(sequences, phred_seqs):
    sequences = [seq.upper() for seq in sequences]
    phred_seqs = length_correction(sequences=sequences, phred_seqs=phred_seqs)
    index_to_base = {0: 'A', 1: 'T', 2: 'C', 3: 'G', 4: 'N'}
    base_to_index = {'A': 0, 'T': 1, 'C': 2, 'G': 3, 'N':4}
    max_seq_length = len(sequences[0])
    N_sequences = len(sequences)
    consensus = ''
    for j in range(max_seq_length):
        conflation_vector = np.zeros((5,), dtype=object)
        position_matrix = np.zeros((N_sequences, 5), dtype=object)
        for i in range(N_sequences):
            base = sequences[i][j]
            error_prob = phred_to_prob(phred_seqs[i][j])
            try:
                position_matrix[i, 4] = Decimal(1e-17)
                position_matrix[i, base_to_index[base]] = Decimal(1) - Decimal(error_prob)
                position_matrix[i][position_matrix[i] == 0.0] = Decimal(error_prob)-Decimal(1e-17)/3
            except KeyError:
                position_matrix[i] = Decimal(1e-17)/4
                position_matrix[i, 4] = Decimal(1)-Decimal(1e-17)
        for base_index in range(0, 5):
            n = np.prod(position_matrix[:, base_index])
            m = np.prod(1-position_matrix[:, base_index])
            conflation_vector[base_index] = n/(n+m)
        max_indices = np.argwhere(conflation_vector == np.amax(conflation_vector))
        max_bases = [index_to_base[array_index] for array_index in max_indices[0]]
        if len(max_bases) == 1:
            consensus += max_bases[0]
        else:
            consensus += 'N'
    return consensus


def MSAligner(list_of_sequences, path_to_mafft):
    """ This function takes a list of sequences an returns gapped sequences list as
    calculated by MAFFT mutliple sequence aligner """
    mafft_param = '--op'
    mafft_param_value = '10'
    tmp_file = 'tmp' + ''.join(map(str, random.sample(range(0, 10), 9))) + '.fasta'
    with open(tmp_file, 'w') as inputFastaFile:
        for i in range(len(list_of_sequences)):
            inputFastaFile.write('>' + str(i) + '\n' + list_of_sequences[i] + '\n')
    process = subprocess.Popen([path_to_mafft, mafft_param, mafft_param_value, tmp_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    alignedSequences = [x.upper() for x in stdout.decode('utf-8').split('\n')[1:-1:2]]
    os.system('rm -rf ' + tmp_file)
    return alignedSequences


####################################################################################################################


def consensus_generator(filt_bam, path_to_mafft, path_to_samtools):
    only_edited_bam = filt_bam[:filt_bam.find('.bam')] + '.only.edited.bam'
    os.system(
        "{} view -h {} | grep _e | samtools view -Sb > {}".format(path_to_samtools, filt_bam, only_edited_bam))
    only_edited_sorted_bam = only_edited_bam[:only_edited_bam.find('.bam')] + '.sort.by.miR.bam'
    os.system('{} sort  {} -o {} -O bam'.format(path_to_samtools, only_edited_bam, only_edited_sorted_bam))
    sortedBamFile = pysam.AlignmentFile(only_edited_sorted_bam, 'rb')
    fasta_file = open('consensus.fasta', 'w')
    unvalid_haplotypes = open('unvalids.txt', 'a')
    alignments_list = []
    for alignment in sortedBamFile.fetch(until_eof=True):
        if alignment.reference_name.find('_e') != -1:
            if len(alignments_list) == 0:
                alignments_list.append(alignment)
            else:
                if alignment.reference_name == alignments_list[-1].reference_name:
                    alignments_list.append(alignment)
                else:
                    ref_name = alignments_list[0].reference_name
                    class_sequences = []
                    class_phreds = []
                    for sub_alignment in alignments_list:
                        class_sequences.append(sub_alignment.query_sequence)
                        class_phreds.append(sub_alignment.qual)
                    aligned_class_sequences = MSAligner(class_sequences, path_to_mafft=path_to_mafft)
                    consensus = consensus_calculator(aligned_class_sequences, class_phreds).strip('N')
                    # Write consensus file into a fasta file with all consensus sequences
                    if 'N' not in consensus:
                        fasta_file.write('>' + ref_name + '\n' + consensus + '\n')
                    else:
                        unvalid_haplotypes.write(ref_name + '\t' + 'Consensus sequence contains ambiguous characters' + '\n')
                    alignments_list = [alignment]
    ref_name = alignments_list[0].reference_name
    class_sequences = []
    class_phreds = []
    for sub_alignment in alignments_list:
        class_sequences.append(sub_alignment.query_sequence)
        class_phreds.append(sub_alignment.qual)
    aligned_class_sequences = MSAligner(class_sequences, path_to_mafft=path_to_mafft)
    consensus = consensus_calculator(aligned_class_sequences, class_phreds).strip('N')
    # Write consensus file into a fasta file with all consensus sequences
    if 'N' not in consensus:
        fasta_file.write('>' + ref_name + '\n' + consensus + '\n')
    else:
        unvalid_haplotypes.write(ref_name + '\t' + 'Consensus sequence contains ambiguous characters' + '\n')
    fasta_file.close()
    os.system('rm -rf ' + only_edited_bam)
    os.system('rm -rf ' + only_edited_sorted_bam)

