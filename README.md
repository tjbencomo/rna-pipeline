# rna-pipeline
RNA-Seq analysis pipeline for the [Lee Lab](http://leelab.stanford.edu/) at Stanford University. The pipeline is designed to mimic dnanexus function on Stanford's Sherlock HPC cluster. The pipeline receives as input paired-end FASTQ files and performs alignment as well as several analysis functions.

## Pipeline Design
The pipeline consists of 3 different steps:
1. [Alignment](https://github.com/tjbencomo/rna-pipeline/blob/master/README.md#star-aligner)
2. [RNA Expression Quantification](https://github.com/tjbencomo/rna-pipeline/blob/master/README.md#rsem-expression-quantification)
3. [Differential Expression Analysis](https://github.com/tjbencomo/rna-pipeline/blob/master/README.md#deseq2-differential-expression-analysis)

Alignment is first completed using the [STAR Aligner](https://github.com/alexdobin/STAR). Its settings mimic those found in ENCODE's [long-rna-seq-pipeline](https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/f9ff54ddf1d955382a1f0aa50b55c8627702f6e1/dnanexus/align-star-pe/resources/usr/bin/lrna_align_star_pe.sh). RNA Expression Quantification and Gene Expression Analysis are performed in parallel. RNA Expression Quantification is computed using [RSEM](https://github.com/deweylab/RSEM) software package. Gene Expression Analysis is performed first by counting RNA reads with [HTSeq-count](http://htseq.readthedocs.io/en/master/count.html) and then analyzing counts with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). 

## Individual Pipeline Components
### STAR Aligner
`star-aligner.py` performs sequence alignment with the STAR aligner

`star-aligner.py` acts as a wrapper for `star.sh` the bash script that executes the star aligner. `star-aligner.py` submits a slurm sbatch job to run the STAR aligner.

#### Inputs

`-wd`|`--workingDirectory` Working Directory: Where all output files will be located

`-f1`|`--fastq1` FASTQ file 1: 1 of the two paired end fastq files to align

`-f2`|`--fastq2` FASTQ file 2: 1 of the two paired end fastq files to align

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

### RSEM Expression Quantification
`rsem-calculate.py` performs rna expression quantification with the RSEM software package. 

`rsem-calculate.py` is a wrapper for rsem.sh, which actually calls the rsem-calculate-expression program. `rsem-calculate.py` submits a SLURM sbatch job of rsem.sh

#### Inputs

`-I`|`--input-bam-file` Input Bam: RNA-Seq aligned transcriptome bam file

`-rDir`|`--rsem-ref-directory` RSEM Reference Directory: Prebuilt files RSEM needs to run

`wd`|`--working-directory` Working Directory: Where all output files are located 

#### Outputs

(Files are prefixed with prefix computated from input bam file - contains filename until 'Aligned')

* `.genes.results` Gene level expression estimates
* `.isoforms.results` Isoform level expression estimates
* `.stat` Folder containing model statistics
* `.transcript.bam` Read alignments in transcript coordinates

Example command: 
`python rsem-calculate.py -wd /scratch/users/tbencomo/RNA_seq/pipeline-tests/ -I /scratch/users/tbencomo/RNA_seq/pipeline-tests/SG13_004_004_CGCTCATT-ATAGAGGC_Aligned.toTranscriptome.out.bam  -rDir /scratch/users/tbencomo/RNA_seq/refs/out/rsem`

### HTSeq Read Counts
`htseq-launch.py` performs a count of all the read sequences from the aligned sorted by coordinate bam input file via the HTSeq software package. 

`htseq-launch.py` is a wrapper for `htseq.sh,` which actually calls the htseq-count program. `htseq-launch.py` submits a SLURM sbatch job of htseq.sh

#### Inputs

`-I`|`--input-bam-file` Input Bam: RNA-Seq aligned sorted by coordinates bam file

`wd`|`--working-directory` Working Directory: Where all output files are located 

#### Outputs

(Files are prefixed with prefix computated from input bam file - contains filename until 'Aligned')


* `_counts.txt` Read counts for each gene

Example command: 
`python htseq-launch.py -I SG13_004_004_CGCTCATT-ATAGAGGC_Aligned.sortedByCoord.out.bam -wd /scratch/users/tbencomo/RNA_seq/pipeline-tests/`

### DESeq2 Differential Expression Analysis
`DESeq-launch.py` performs differential expression analysis on several samples, reporting which genes are upregulated and downregualted. 

`DESeq-launch.py` is a wrapper for desq.sh. deseq.sh is a shell script to submit a SLURM sbatch job. The actual analysis is performed by an R script `DESeq-analysis.R`

NOTE: As of now, the R script determines if a sample is tumor or normal based on its filename. If the file contains 'NS' or 'ctrl', it is classified as normal. Otherwise, it is classified as a tumor. In the future, the user will be able to enter strings for the program to classify each condition.
#### Inputs

`-counts`|`--counts-Directory` Directory containing HTSeq-count counts files
`results`|`--results-Directory` Directory to store results files

#### Outputs

`unfiltered_normal_tumor_results.csv` CSV file containing all genes. Not filtered for pvalues or log2FoldChange levels
`upregulated_genes.csv` CSV file containing all genes with log2FoldChange scores > 0 and with adjusted p-values <= .05
`downregulated_genes.csv` CSV file containing all genes with log2FoldChange scores < 0 and with adjusted p-values <= .05

Example command:
`python DESeq-launch.py -counts /scratch/PI/carilee/NatComm-Analysis/counts/ -results /scratch/PI/carilee/NatComm-Analysis/results/`
