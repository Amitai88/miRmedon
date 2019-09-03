import miRmedon_src
from miRmedon_src.alignment_to_emiRbase import alignment_to_emiRbase
from miRmedon_src.alignments_filter import main_alignments_filter
from miRmedon_src.count_filter import main_counts_filter
from miRmedon_src.consensus_generator import consensus_generator
from miRmedon_src.clusters_filter_bowtie import main_clusters_filter_bowtie
from miRmedon_src.clusters_filter_seqmap import main_clusters_filter_seqmap
from miRmedon_src.counter import counter
from miRmedon_src.MC_P_est import monte_catlo_p_estimation
from miRmedon_src.recount import recount
import argparse
import os


miRmedon_dir = miRmedon_src.__file__.replace('__init__.py', '')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '-fastq', type=str,
                        help='Path to fastq file')
    parser.add_argument('-star', '-star', type=str,
                        help='Path to Star')
    parser.add_argument('-t', '-threads', type=str,
                        help='Number of threads for Star alignment')
    parser.add_argument('-s', '-soft_clipping_threshold', type=int, default=2,
                        help='Soft-cliiped bases threshold')
    parser.add_argument('-x', '-total_modifications_threshold', type=int, default=2,
                        help='Total number of modifications threshold (sofft-clipped+mismatches')
    parser.add_argument('-c', '-counts_threshold', type=int, default=10,
                        help='Counts threshold')
    parser.add_argument('-samtools', '-samtools', type=str,
                        help='Path to samtools.')
    parser.add_argument('-mafft', '-mafft', type=str,
                        help='Path to Mafft.')
    parser.add_argument('-seqmap', '-seqmap', type=str, required=False, default=None,
                        help='Path to Seqmap.')
    parser.add_argument('-bowtie', '-bowtie', type=str, required=False, default=None,
                        help='Path to Bowtie.')
    parser.add_argument('-G', '-genome_fasta', type=str,
                        help='Genome FASTA file or bowtie index.')
    parser.add_argument('-T', '-transcriptome_fasta', type=str,
                        help='Transcriptome FASTA file or bowtie index.')
    parser.add_argument('-r', '-resamples', type=int, default=10000,
                        help='Number of resamples for Monte-Carlo p-value estimation')
    args = parser.parse_args()

    vars_args = vars(args)

    if (vars_args['seqmap'] != None and vars_args['bowtie'] != None):
        raise Exception('You could only choose one of the following aligners:\n '
                        '          - Bowtie\n '
                        '          - Seqmap\n '
                        '          Both were given to miRmedon')

    fastq_file = args.f
    path_to_star = args.star
    threads = args.t
    star_ref_dir = miRmedon_dir + 'star_ref/'
    soft_clip_thresh = args.s
    total_mod_thresh = args.x
    counts_thresh = args.c
    path_to_samtools = args.samtools
    path_to_mafft = args.mafft
    e_miRbase = miRmedon_dir + 'e_miRbase.pkl'
    path_to_seqmap = args.seqmap
    path_to_bowtie = args.bowtie
    genome_index = args.G
    transcriptome_index = args.T
    resamples = args.r
    bam = '_Aligned.out.bam'
    st_filt_bam = '_Aligned.out.filt.st.bam'
    bycount_filt_bam = '_Aligned.out.filt.bycount.bam'
    final_filt_bam = '_Aligned.out.filt.final.bam'

    alignment_to_emiRbase(fastq_file, path_to_star, threads, star_ref_dir)
    main_alignments_filter(bam, e_miRbase, soft_clip_thresh, total_mod_thresh)
    main_counts_filter(st_filt_bam, path_to_samtools, counts_thresh)
    consensus_generator(bycount_filt_bam, path_to_mafft, path_to_samtools)
    if vars_args['seqmap'] != None:
        miRs_bed = miRmedon_dir + 'miRbase.seqmap.bed'
        main_clusters_filter_seqmap(bycount_filt_bam, e_miRbase, path_to_seqmap, genome_index, transcriptome_index, miRs_bed, path_to_samtools)
    elif vars_args['bowtie'] != None:
        miRs_bed = miRmedon_dir + 'miRbase.bowtie.bed'
        main_clusters_filter_bowtie(bycount_filt_bam, e_miRbase, path_to_bowtie, genome_index, transcriptome_index, miRs_bed, path_to_samtools)
    counter(final_filt_bam)
    monte_catlo_p_estimation(final_filt_bam, e_miRbase, resamples)
    recount(e_miRbase)

    os.system('rm -rf _Log.final.out _Log.out _Log.progress.out _Unmapped.out.mate1 _SJ.out.tab')
    os.system('rm -rf _Aligned.out.*')
    os.system('rm -rf *.sam')
    os.system('rm -rf *.eland')
    os.system('rm -rf consensus.fasta reduced_consensus.fasta')
    os.system('rm -rf counts.txt')
    os.system('mv crt_counts.txt counts.txt')
