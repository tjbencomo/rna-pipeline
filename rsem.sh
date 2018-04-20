#!/bin/bash
#SBATCH --job-name=rsem-internal
#SBATCH --nodes=1
#SBATCH --mem=30000
#SBATCH --time=0-12:00:00
#SBATCH --mail-type=END

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -b|--aligned-transcriptome-bam)
    BAM="$2"
    shift # past argument
    shift # past value
    ;;
    -prefix|--filePrefixName)
    PREFIX="$2"
    shift # past argument
    shift # past value
    ;;
    -rDir|--referenceDirectory)
    REF_DIR="$2"
    shift # past argument
    shift # past value
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

echo BAM        = "${BAM}"
echo PREFIX	= "${PREFIX}"
echo REF_DIR	= "${REF_DIR}"
echo DEFAULT         = "${DEFAULT}"

rsem-calculate-expression --alignments --paired-end --ci-memory 8000 --num-threads $SLURM_CPUS_ON_NODE $BAM $REF_DIR $PREFIX

