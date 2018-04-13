#!/bin/bash
#SBATCH --job-name=htseq-internal
#SBATCH --nodes=1
#SBATCH --mem=1000
#SBATCH --time=03:00:00
#SBATCH --mail-type=END




POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -b|--aligned-sortedCoordinate-bam)
    BAM="$2"
    shift # past argument
    shift # past value
    ;;
    -prefix|--filePrefixName)
    PREFIX="$2"
    shift # past argument
    shift # past value
    ;;
    -pipeDir|--pipeline-directory)
    PIPE_DIR="$2"
    shift
    shift
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
echo PREFIX     = "${PREFIX}"
echo PIPE_DIR	= "${PIPE_DIR}"
echo DEFAULT         = "${DEFAULT}"

source $PIPE_DIR/python/bin/activate

READS=$BAM
ANNOTS=$PIPE_DIR/refs/gencode.v27.annotation.gtf.gz

$PIPE_DIR/python/bin/htseq-count -f bam $READS $ANNOTS > "${PREFIX}"_counts.txt

deactivate


