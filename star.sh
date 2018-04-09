#!/bin/bash
#SBATCH --job-name=star-internal
#SBATCH --output=/home/users/tbencomo/out/star-internal.%j.out
#SBATCH --error=/home/users/tbencomo/errout/star-internal.%j.err
#SBATCH --nodes=1
#SBATCH --mem=60000
#SBATCH --time=0-03:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=tbencomo@stanford.edu


POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -f1|--fastq1)
    FASTQ1="$2"
    shift # past argument
    shift # past value
    ;;
    -f2|--fastq2)
    FASTQ2="$2"
    shift # past argument
    shift # past value
    ;;
    -gDir|--genomeDir)
    GENOME_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    -prefix|--outFilePrefix)
    PREFIX="$2"
    shift # past argument
    shift # past argument
    ;;
    -wd|--workingDirectory)
    WORK_DIR="$2"
    shift # past argument
    shift # past argument
    ;;
    --default)
    DEFAULT=YES
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

echo FASTQ1  = "${FASTQ1}"
echo FASTQ2     = "${FASTQ2}"
echo GENOME DIR    = "${GENOME_DIR}"
echo PREFIX	= "${PREFIX}"
echo WORK DIR	= "${WORK_DIR}"
echo DEFAULT         = "${DEFAULT}"


cd $WORK_DIR

STAR --genomeDir $GENOME_DIR --readFilesIn $FASTQ1 $FASTQ2 --runThreadN $SLURM_CPUS_ON_NODE --genomeLoad NoSharedMemory --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --sjdbScore 1 --limitBAMsortRAM 10000000000 --outFileNamePrefix $WORK_DIR/$PREFIX_

