import os

def alignment_to_emiRbase(fastq_file, path_to_star, threads, star_ref_dir, path_to_samtools):
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
    os.system('mv _Aligned.out.bam _Aligned.out.FCRC.bam')
    os.system('{} view -b -h -F 16 _Aligned.out.FCRC.bam > _Aligned.out.bam'.format(path_to_samtools))
