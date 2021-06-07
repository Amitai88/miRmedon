from functools import reduce
from decimal import Decimal
import numpy as np
import random
import subprocess
import os
import pysam

################################################### Auxilary functions #############################################
    
    
def qual_to_prob(qual_score):
    """This function takes an ascii character and as an input and returns 
    the error probability"""
    return 10**(-qual_score/10)


def length_correction(sequences, qual_lists):
    if all(isinstance(x, str) for x in qual_lists):
        qual_lists = [[ord(c)-33 for c in seq] for seq in qual_lists] 
    corrected_qual_lists = []
    for i in range(len(sequences)):
        startGapLen = 0
        endGapLen = 0
        for j in range(len(sequences[i])):
            if sequences[i][j] not in ['A', 'T', 'C', 'G', 'N']:
                startGapLen += 1
            else:
                break
        for j in range(len(sequences[i]))[::-1]:
            if sequences[i][j] not in ['A', 'T', 'C', 'G', 'N']:
                endGapLen += 1
            else:
                break 
        corrected_qual_lists.append([99]*startGapLen + qual_lists[i] + [99]*endGapLen)
    return corrected_qual_lists


def consensus_calculator(sequences, qual_lists):
    sequences = [seq.upper() for seq in sequences]
    qual_lists = length_correction(sequences=sequences, qual_lists=qual_lists)
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
            error_prob = qual_to_prob(qual_lists[i][j])
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
        "{} view -h {} | grep _e | {} view -Sb > {}".format(path_to_samtools, filt_bam, path_to_samtools, only_edited_bam))
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
                    class_quals = []
                    for sub_alignment in alignments_list:
                        class_sequences.append(sub_alignment.query_alignment_sequence)
                        class_quals.append(list(sub_alignment.query_alignment_qualities))
                        #class_sequences.append(sub_alignment.query_sequence)
                        #class_quals.append(list(sub_alignment.query_qualities))
                    if len(class_sequences) > 1e6:
                        random_indices = np.random.choice(np.arange(len(class_sequences)), 1000)
                        class_sequences = list(np.array(class_sequences)[random_indices])
                        class_quals = list(np.array(class_quals)[random_indices])
                    aligned_class_sequences = MSAligner(class_sequences, path_to_mafft=path_to_mafft)
                    consensus = consensus_calculator(aligned_class_sequences, class_quals).strip('N')
                    # Write consensus file into a fasta file with all consensus sequences
                    if 'N' not in consensus:
                        fasta_file.write('>' + ref_name + '\n' + consensus + '\n')
                    else:
                        unvalid_haplotypes.write(ref_name + '\t' + 'Consensus sequence contains ambiguous characters' + '\n')
                    alignments_list = [alignment]
    ref_name = alignments_list[0].reference_name
    class_sequences = []
    class_quals = []
    for sub_alignment in alignments_list:
        class_sequences.append(sub_alignment.query_alignment_sequence)
        class_quals.append(list(sub_alignment.query_alignment_qualities))
        #class_sequences.append(sub_alignment.query_sequence)
        #class_quals.append(sub_alignment.qual)
    if len(class_sequences) > 1e6:
        random_indices = np.random.choice(np.arange(len(class_sequences)), 100)
        class_sequences = list(np.array(class_sequences)[random_indices])
        class_quals = list(np.array(class_quals)[random_indices])
    aligned_class_sequences = MSAligner(class_sequences, path_to_mafft=path_to_mafft)
    consensus = consensus_calculator(aligned_class_sequences, class_quals).strip('N')
    # Write consensus file into a fasta file with all consensus sequences
    if 'N' not in consensus:
        fasta_file.write('>' + ref_name + '\n' + consensus + '\n')
    else:
        unvalid_haplotypes.write(ref_name + '\t' + 'Consensus sequence contains ambiguous characters' + '\n')
    fasta_file.close()
    os.system('rm -rf ' + only_edited_bam)
    os.system('rm -rf ' + only_edited_sorted_bam)

