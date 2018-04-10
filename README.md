# rna-pipeline
RNA-Seq analysis pipeline for the [Lee Lab](http://leelab.stanford.edu/) at Stanford University. The pipeline is designed to mimic dnanexus function on Stanford's Sherlock HPC cluster. The pipeline receives as input paired-end FASTQ files and performs alignment as well as several analysis functions.

## Pipeline Design
The pipeline consists of 3 different steps:
1. Alignment
2. RNA Expression Quantification
3. Gene Expression Analysis

Alignment is first completed using the [STAR Aligner](https://github.com/alexdobin/STAR). Its settings mimic those found in ENCODE's [long-rna-seq-pipeline](https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/f9ff54ddf1d955382a1f0aa50b55c8627702f6e1/dnanexus/align-star-pe/resources/usr/bin/lrna_align_star_pe.sh). RNA Expression Quantification and Gene Expression Analysis are performed in parallel. RNA Expression Quantification is computed using [RSEM](https://github.com/deweylab/RSEM) software package. Gene Expression Analysis is performed first by counting RNA reads with [HTSeq-count](http://htseq.readthedocs.io/en/master/count.html) and then analyzing counts with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). 

## Individual Pipeline Components
### STAR Aligner
`star-aligner.py` performs sequence alignment with the STAR aligner

`star-aligner.py` acts as a wrapper for `star.sh` the bash script that executes the star aligner. `star-aligner.py` submits a slurm sbatch job to run the STAR aligner.

#### Inputs

`-wd`|`--workingDirectory` Working Directory: Where all output files will be located

`-f1`|`--fastq1` FASTQ file 1: 1 of the two paired end fastq files to align

`-f2`|`--fast21` FASTQ file 2: 1 of the two paired end fastq files to align

`-gDir`|`--genomeDirectory` Genome Directory: Location of reference files for star aligner

`-prefix`|`--outFilePrefix` Output File Prefix: file prefix that will be appended to start of output files. Optional. By default, the prefix is set to the filename of FASTQ1

#### Outputs

(Files are prefixed with specified or default prefix above)

* `_Aligned.sortedByCoord.out.bam` Aligned BAM file sorted by coordinate. Used for programs such as HTSeq-count
* `_Aligned.toTranscriptome.out.bam` Aligned BAM of translated coordinates. Used for programs such as RSEM
* `_Log.out` Main file with information about run
* `_Log.progress.out` Job progress statistics

Example command:
`python star-aligner.py -wd /scratch/users/tbencomo/RNA_seq/pipeline-tests -gDir /scratch/users/tbencomo/RNA_seq/refs/out -f1 /scratch/users/tbencomo/RNA_seq/input_files/SG13_004_004_CGCTCATT-ATAGAGGC_R1.fastq -f2 /scratch/users/tbencomo/RNA_seq/input_files/SG13_004_004_CGCTCATT-ATAGAGGC_R2.fastq`
