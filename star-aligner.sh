#!/bin/bash
#SBATCH --job-name=star-aligner
#SBATCH --output=/home/users/$USER/out/star_aligner.%j.out
#SBATCH --error=/home/users/$USER/errout/star_aligner.%j.err
#SBATCH --nodes=1
#SBATCH --mem=60000
#SBATCH --time=0-05:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=tbencomo@stanford.edu
#SBATCH --workdir=/scratch/users/tbencomo/RNA_seq/input_files

date
STAR --genomeDir /scratch/users/tbencomo/RNA_seq/refs/out --readFilesIn SG14_006_006_GAGATTCC-ATAGAGGC_R1.fastq SG14_006_006_GAGATTCC-ATAGAGGC_R2.fastq --runThreadN 1 --genomeLoad NoSharedMemory --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --sjdbScore 1 --limitBAMsortRAM 10000000000 --outFileNamePrefix /scratch/users/tbencomo/RNA_seq/align-files/SG14_



