#!/bin/bash
#SBATCH --job-name=shell-executer
#SBATCH --output=/home/users/tbencomo/out/star_executer.%j.out
#SBATCH --error=/home/users/tbencomo/errout/star_executer.%j.err
#SBATCH --nodes=1
#SBATCH --mem=1000
#SBATCH --time=0-00:05:00
#SBATCH --mail-type=END
#SBATCH --mail-user=tbencomo@stanford.edu

cd $HOME/rna-pipeline

echo "HELLO DID THIS PRINT OUT TO A STANDARD CONSOLE WINDOW? THIS WAS LAUNCHED FROM PYTHON"
