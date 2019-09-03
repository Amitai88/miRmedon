import pysam
from collections import Counter

def counter(bam):
    bam_input = pysam.AlignmentFile(bam, 'rb')
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
    with open('counts.txt', 'w') as counts_file:
        for haplotype, count in list(haplotypes_counter.items()):
            counts_file.write(haplotype + '\t' + str(count) + '\n')


