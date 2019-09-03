import os
from collections import Counter
import pysam


def haplotype_counter_func(st_filt_bam):
    bam_input = pysam.AlignmentFile(st_filt_bam, 'rb')
    haplotypes_counter = Counter()
    alignments_list = []
    for alignment in bam_input.fetch(until_eof=True):
        if len(alignments_list) == 0:
            alignments_list.append(alignment)
        else:
            if alignment.query_name == alignments_list[-1].query_name:
                alignments_list.append(alignment)
            else:
                tmp_counter = Counter([x.reference_name for x in alignments_list])
                tmp_counter = {k: v/len(tmp_counter) for k, v in tmp_counter.items()}
                haplotypes_counter += Counter(tmp_counter)
                alignments_list = [alignment]
    tmp_counter = Counter([x.reference_name for x in alignments_list])
    tmp_counter = {k: v / len(tmp_counter) for k, v in tmp_counter.items()}
    haplotypes_counter += Counter(tmp_counter)
    haplotypes_counter = dict(haplotypes_counter)
    return haplotypes_counter


def write_counts(counts_dict, count_thresh):
    with open('unvalids.txt', 'w') as unvalids_file:
        for haplotype, count in counts_dict.items():
            if count < count_thresh:
                unvalids_file.write(haplotype + '\t' + 'low count' + '\n')


def BAM_filter_counts(st_filt_bam, path_to_samtools):
    filt_bycount_bam = st_filt_bam[:st_filt_bam.find('.filt.st.bam')] + '.filt.bycount.bam'
    os.system('cut -f 1 unvalids.txt > unvalids_k1.txt')
    os.system(
        "{} view -h {} | grep -v -w -F -f unvalids_k1.txt | samtools view -Sb > {}".format(path_to_samtools, st_filt_bam, filt_bycount_bam))
    os.system('rm -rf unvalids_k1.txt')


def main_counts_filter(st_filt_bam, path_to_samtools, count_thresh):
    counts_dict = haplotype_counter_func(st_filt_bam)
    write_counts(counts_dict=counts_dict, count_thresh=count_thresh)
    BAM_filter_counts(st_filt_bam, path_to_samtools=path_to_samtools)
