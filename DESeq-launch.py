'''
Author: Tomas Bencomo

DESeq-launch.py wrapper script for DESeq-analysis.R. Final step of the rna-pipeline

Script takes as input: 
1) Filepath for directory with counts data
2) Filepath for directory to store results

The following files are written to the results directory designated above:
1) Unfiltered output file, containing all results
2) Upregulated csv file containing all genes that show upregulated expression and have adjusted pvalues <= .05
3) Downregulated csv file containing all genes that show downregulated expression and have adjusted pvalues <= .05

'''


import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='Launch DESeq-analysis.R')

parser.add_argument('-counts', '--counts-Directory', metavar = '--counts-Directory', type=str, nargs=1, help='Directory containing the read counts files', dest='counts_directory')

parser.add_argument('-results', '--results-Directory', metavar = '--results-Directory', type=str, nargs=1, help='Directory containing the results files from DESeq2 analysis', dest= 'results_directory')

args = parser.parse_args()

#Get environment variable to be user independent
USER = os.environ['USER']

#Fix formating to convert to strings
args.counts_directory = ''.join(args.counts_directory)
args.results_directory = ''.join(args.results_directory)

workdir = "--workdir=" + args.results_directory
output = "--output=/home/users/" + USER + "/out/deseq.%j.out"
error = "--error=/home/users/" + USER + "/errout/deseq.%j.err"
mail = "--mail-user=" + USER + "@stanford.edu" #assumes user's email is structured sherlock_username@stanford.edu

subprocess.call(['sbatch', workdir, output, error, mail, 'deseq.sh', '-counts', args.counts_directory, '-results', args.results_directory])






