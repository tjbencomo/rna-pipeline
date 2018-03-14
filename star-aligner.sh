#!/bin/bash
#SBATCH --job-name=star-aligner
#SBATCH --output=/home/users/tbencomo/out/star_aligner.%j.out
#SBATCH --error=/home/users/tbencomo/errout/star_aligner.%j.err
#SBATCH --nodes=1
#SBATCH --mem=60000
#SBATCH --time=0-03:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=tbencomo@stanford.edu
#SBATCH --workdir=/home/users/tbencomo/rna-pipeline

echo $USER
echo "round 2"
echo "$USER"
echo $1
echo "$1"



STAR --genomeDir $2 --readFilesIn $3 $4 --runThreadN $SLURM_CPUS_ON_NODE --genomeLoad NoSharedMemory --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --sjdbScore 1 --limitBAMsortRAM 10000000000 --outFileNamePrefix $1/$5_



