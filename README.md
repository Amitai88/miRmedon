# miRmedon
Confident detection of A-to-I miRNA editing events in small-RNA seq data

# Description
miRmedon it's a python3 code which takes fastq file of trimmed small-RNA reads, and reports on A-to-I editing events 
detected in the data. In addition, miRmedon reports on read counts and sequence of both edited and non-edited miRNAs. 
miRmedon require several extrenal softwares and python libraries, as descirbed below.

For detailed description see _link to bioRxiv_

# Installation
No installation is required for miRmedon. Please directly use miRmedon.py source file and make sure it is located within the  
same directory as miRmedon_src directory.  

# Dependencies
**External softwares:**
- STAR - https://github.com/alexdobin/STAR 
- MAFFT - https://mafft.cbrc.jp/alignment/software/source.html
- Samtools - http://www.htslib.org/download/ 
- Seqmap - http://www-personal.umich.edu/~jianghui/seqmap/ \
or
- Bowtie - http://bowtie-bio.sourceforge.net/index.shtml 

**Python libraries:** 
- pysam - https://pysam.readthedocs.io/en/latest/api.html 
- numpy 
- scipy

**Reference genome and transcriptome:**
- GRCh38.p12 - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh38.p12.genome.fa.gz
- gencode.v31.transcripts - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.transcripts.fa.gz

# Usage
**Parameters** \
**miRmedon takes the following parameters:**
```
-f - path to fastq file
-star – path to STAR 
-samtools – path to samtools
-mafft – path to mafft 
-seqmap – path to seqmap
-bowtie – path to bowtie
-G – path to GRCh38.p12 fasta file or GRCh38.p12 Bowtie index
-T – path to gencode.v31.transcripts fasta file or gencode.v31.transcripts Bowtie index
-t – number of threads to run STAR
-s – number of soft-clipped bases threshold (default = 2)
-x – total number of modifications (soft-clipping and mismatches) threshold (default = 3)
-c – read counts threshold (default = 10)
-r – number of resamples for Monte-Carlo p-value estimation (default = 10000)
```
**Note**
- Only one tool is required for consensus sequence alignment (see paper). Currently, Seqmap and Bowtie are implemented in miRmedon. 
Please indicate the path for one alignment tool, and leave the additional path parameter empty. If Seqmap paramter was indicated,
-G and -T paramters should point on the required fasta files (GRCh38.p12.genome.fa and gencode.v31.transcripts.fa).
However, if bowtie paramter was indicated, -G and -T should point on bowtie indices for GRCh38.p12.genome.fa and gencode.v31.transcripts.fa
respectively.

**Command line examples** 

In order to apply miRmedon with default paramters and Bowtie, use the following command line example: 
```
python3 miRmedon.py -f path_to_fastq_file -star path_to_star -t number_of_threads -samtools path_to_samtools
-mafft path_to_mafft -bowtie path_to_bowtie -G path_to_GRCh38.p12.genome_bowtie_index -T path_to_bowtie_gencode.v31.transcripts_bowtie_index
```

In order to apply miRmedon with default paramters and Seqmap, use the following command line example: 
```
python3 miRmedon.py -f path_to_fastq_file -star path_to_star -t number_of_threads -samtools path_to_samtools
-mafft path_to_mafft -seqmap path_to_seqmap -G path_to_GRCh38.p12.genome.fa -T path_to_gencode.v31.transcripts.fa
```

# Output
miRmedon will generate two tab-delimited report files: 
- editing_info.txt - contain inforamtion on editing sites detected, including: mature miRNA name, position, editing levels, 
lower and upper confidence interval for point estimated editing levels and a p-value (see paper).

- counts.txt - contain information on read counts of both edited and non-edited miRNAs, toghether with reference sequence and editing sites
within the reference sequence. Non-edited miRNA will have the extension "_ne", while edited miRNAs will have the extension "_eX" with X as
an arbitrary number.

# Contact 
Please contact me with regard to any problem or suggestion at mordeami@post.bgu.ac.il
