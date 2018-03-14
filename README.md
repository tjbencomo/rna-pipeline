# rna-pipeline
RNA-Seq analysis pipeline for the Lee Lab at Stanford University's School of Medicine. The pipeline consists of several individual bash scripts for each step in the pipeline that are wrapped in a python script. The python script glues the scripts together to allow for analysis of several files in a single command.
## Scripts
Note: Inputs must be in the precise order they are listed in the documentation. Named parameters will be supported in a future version.
### STAR Aligner
`star-aligner.sh` Performs sequence alignment with the STAR aligner

It takes the following inputs:
1. Working directory
2. genomeDir directory
3. Input FASTQ file 1
4. Input FASTQ file 2
5. Output file prefix

It outputs the following files (All files will be prefixed with the Output file prefix parameter):
* `_Aligned.sortedByCoord.out.bam` Aligned BAM file sorted by coordinate. Used for programs such as HTSeq-count
* `_Aligned.toTranscriptome.out.bam` Aligned BAM of translated coordinates. Used for programs such as RSEM
* `_Log.out` Main file with information about run
* `_Log.progress.out` Job progress statistics

Example command:
`sbatch star-aligner.sh /scratch/users/tbencomo/RNA_seq/pipeline-tests /scratch/users/tbencomo/RNA_seq/refs/out /scratch/users/tbencomo/RNA_seq/input_files/SG13_004_004_CGCTCATT-ATAGAGGC_R1.fastq /scratch/users/tbencomo/RNA_seq/input_files/SG13_004_004_CGCTCATT-ATAGAGGC_R2.fastq SG14`

This command enqueues the star-aligner script onto Sherlock with the following options:
1. `/scratch/users/tbencomo/RNA_seq/pipeline-tests` Working directory
2. `/scratch/users/tbencomo/RNA_seq/refs/out` genomeDir directory required by STAR
3. `/scratch/users/tbencomo/RNA_seq/input_files/SG13_004_004_CGCTCATT-ATAGAGGC_R1.fastq` FASTQ file
4. `/scratch/users/tbencomo/RNA_seq/input_files/SG13_004_004_CGCTCATT-ATAGAGGC_R2.fastq` FASTQ file
5. `SG14` The output prefix. All output files will be prefixed with SG14_
