import os
import pickle


def length_filter(e_miRbase_path):
    e_miRbase = pickle.load(open(e_miRbase_path, 'rb'))
    unvalids_file = open('unvalids.txt', 'a')
    consensus_file = open('consensus.fasta', 'r')
    consensus = {line.lstrip('>').rstrip('\n'):next(consensus_file).rstrip('\n') for line in consensus_file 
                 if line.startswith('>')}
    consensus_file = open('consensus.fasta', 'w')
    for title, seq in consensus.items():
        miR = title[:title.find('_')]
        ref_seq = e_miRbase[miR][title]['sequence']
        rel_len = len(seq)/len(ref_seq)
        if 0.95 <= rel_len <= 1.05:
            consensus_file.write('>' + title + '\n' + seq + '\n')
        else:
            unvalids_file.write(title + '\t' + 'Consensus sequence too short/long' + '\n')


def alignment_to_genome(path_to_bowtie, path_to_genome):
    os.system('{} 1 consensus.fasta {} genome.eland /output_all_matches'.format(path_to_bowtie, path_to_genome))
    return None


def get_miRs_coors(path_to_miRs_bed):
    miRs_bed_file = open(path_to_miRs_bed, 'r')
    miRs_coors = {}
    for line in miRs_bed_file:
        l = line.split('\t')
        miR = l[3].rstrip('\n') 
        chromosome = l[0]
        startIndex = int(l[1])
        endIndex = int(l[2])
        if miR not in miRs_coors.keys():
            miRs_coors[miR] = {chromosome: [(startIndex, endIndex)]}
        elif miR in miRs_coors.keys():
            if chromosome not in miRs_coors[miR].keys():
                miRs_coors[miR].update({chromosome: [(startIndex, endIndex)]})
            elif chromosome in miRs_coors[miR].keys():
                miRs_coors[miR][chromosome] += [(startIndex, endIndex)]
    return miRs_coors


def get_chr_name(string):
    chromosome = string.split(' ')
    if chromosome[0].startswith('chr'):
        chromosome = chromosome[0]
    elif chromosome[0] == chromosome[1]:
        chromosome = chromosome[0]
    else:
        chromosome = ' '.join(chromosome)
    return chromosome


def filter_genome_aligned_haplotypes(path_to_bowtie,
                                     path_to_genome,
                                     path_to_miRs_bed):
    alignment_to_genome(path_to_bowtie=path_to_bowtie, path_to_genome=path_to_genome)
    miRs_coordinates = get_miRs_coors(path_to_miRs_bed=path_to_miRs_bed)
    genome_aligned_haplotypes = set()
    unvalid_genome_aligned = []
    os.system("tail -n +2 genome.eland | tr '\t' ';' | sort -t ';' -k 4 > genome.sorted.eland")
    with open('genome.sorted.eland', 'r') as sorted_eland_file:
        for line in sorted_eland_file:
            line = line.rstrip('\n').split(';')
            haplotype = line[3]
            if haplotype not in unvalid_genome_aligned:
                mismatches = int(line[5])
                if mismatches > 0:
                    basename = haplotype[:haplotype.find('_')]
                    #chromosome = line[0].split(' ')[0]
                    chromosome = get_chr_name(line[0])
                    start_index = int(line[1])
                    seq = line[2]
                    end_index = start_index + len(seq)
                    coordinates = range(start_index, end_index)
                    bool_list = []
                    try:
                        for source_coordinates in miRs_coordinates[basename][chromosome]:
                            if set(coordinates).intersection(set(range(source_coordinates[0], source_coordinates[1]))):
                                bool_list.append(True)
                            else:
                                bool_list.append(False)
                        if False in bool_list:
                            unvalid_genome_aligned.append(haplotype)
                    except KeyError:
                        unvalid_genome_aligned.append(haplotype)
                elif mismatches == 0:
                    unvalid_genome_aligned.append(haplotype)
            genome_aligned_haplotypes.add(haplotype)

    return genome_aligned_haplotypes, unvalid_genome_aligned


def alignment_to_transcriptome(path_to_bowtie, path_to_transcriptome):
    os.system('{} 1 reduced_consensus.fasta {} transcriptome.eland /output_all_matches'.format(path_to_bowtie, path_to_transcriptome))
    return None


def filter_transcriptome_aligned_haplotypes(path_to_bowtie,
                                            path_to_transcriptome,
                                            haplotypes_to_skip):
    with open('consensus.fasta', 'r') as in_file, open('reduced_consensus.fasta', 'w') as out_file:
        for line in in_file:
            if line.startswith('>'):
                haplotype = line.lstrip('>').rstrip('\n')
                if haplotype not in haplotypes_to_skip:
                    out_file.write(line + next(in_file))
    alignment_to_transcriptome(path_to_bowtie=path_to_bowtie, path_to_transcriptome=path_to_transcriptome)
    try:
        transcriptome_aligned_haplotypes = list(set([line.split('\t')[3] for line in 
                                                     open('transcriptome.eland', 'r')
                                                     if not line.startswith('trans_id')]))
    except FileNotFoundError:
        transcriptome_aligned_haplotypes = []
    return transcriptome_aligned_haplotypes


def write_aligned_haplotypes(path_to_bowtie,
                             path_to_genome,
                             path_to_transcriptome,
                             path_to_miRs_bed):
    genome_aligned_haplotypes, unvalid_genome_aligned = filter_genome_aligned_haplotypes(path_to_bowtie=path_to_bowtie,
                                                                                         path_to_genome=path_to_genome,
                                                                                         path_to_miRs_bed=path_to_miRs_bed)
    unvalid_transcriptome_aligned = filter_transcriptome_aligned_haplotypes(path_to_bowtie=path_to_bowtie,
                                                                            path_to_transcriptome=path_to_transcriptome,
                                                                            haplotypes_to_skip=genome_aligned_haplotypes)
    unvalid_haplotypes = unvalid_genome_aligned + unvalid_transcriptome_aligned
    with open('unvalids.txt', 'a') as file:
        file.writelines([haplotype + '\t' + 'Consensus sequence was aligned to genome/transcriptome' + '\n' for haplotype in unvalid_haplotypes])
    #return notvalid_haplotypes
    return None


def BAM_filter_aligned(filt_bycount_bam, path_to_samtools):
    final_filt_bam = filt_bycount_bam[:filt_bycount_bam.find('.filt.bycount.bam')] + '.filt.final.bam'
    os.system('cut -f 1 unvalids.txt > unvalids_k1.txt')
    os.system(
        "{} view -h {} | grep -v -w -F -f unvalids_k1.txt | {} view -Sb > {}".format(path_to_samtools, filt_bycount_bam, path_to_samtools, final_filt_bam))
    os.system('rm -rf unvalids_k1.txt')



def main_clusters_filter_seqmap(filt_bycount_bam, path_to_emiRbase, path_to_bowtie, path_to_genome, path_to_transcriptome, path_to_miRs_bed, path_to_samtools):
    length_filter(path_to_emiRbase)
    write_aligned_haplotypes(path_to_bowtie=path_to_bowtie,
                             path_to_genome=path_to_genome,
                             path_to_transcriptome=path_to_transcriptome,
                             path_to_miRs_bed=path_to_miRs_bed)
    BAM_filter_aligned(filt_bycount_bam=filt_bycount_bam, path_to_samtools=path_to_samtools)

