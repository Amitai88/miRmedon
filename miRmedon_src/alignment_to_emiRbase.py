import os

def alignment_to_emiRbase(fastq_file, path_to_star, threads, star_ref_dir):
    params = '''--runThreadN {}
                --alignIntronMin 1
                --outFilterMultimapNmax 200
                --outFilterMatchNmin 12
                --outFilterMatchNminOverLread 0.66
                --outFilterMismatchNoverLmax 0.08
                --seedSearchStartLmax 6
                --winAnchorMultimapNmax 2000
                --outFilterMultimapScoreRange 0
                --outSAMtype BAM Unsorted
                --outReadsUnmapped Fastx
                --outFilterMismatchNmax 1
                --outSAMprimaryFlag AllBestScore
                --outWigType None
                --outSAMattributes NH AS NM MD'''.format(threads)
    os.system('{} --genomeDir {} --readFilesIn {} --outFileNamePrefix {} {}'.format(path_to_star, star_ref_dir, fastq_file, '_', params.replace('\n', ' ')))

