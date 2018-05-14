'''
Author: Tomas Bencomo

Script to execute the RSEM software package to quantify rna expression. Acts as a wrapper for the rsem.sh bash script. Submits sbatch call for rsem.sh with SBATCH directives for --workdir, --output, and --error.

This script:
Parses command line options
submits sbatch job for rsem.sh script

Command line parameters:
input bam file
RSEM reference directory
working directory

Determines based off command line parameters:
Sample prefix name based on input file

'''

import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description='Launch rsem-calculate-expression')

parser.add_argument('-I', '--input-bam-file', metavar='--input-bam-file', type=str, nargs=1, help='Aligned transcriptome bam file to quanitify reads', dest='input_bam')

parser.add_argument('-rDir', '--rsem-ref-directory', metavar='-rsem-ref-directory', type=str, nargs=1, help='location of reference files rsem needs to operate', dest='ref_directory')

parser.add_argument('-wd', '--working-directory', metavar='--working-directory', type=str, nargs=1, help='location to store output files', dest='working_directory')

args = parser.parse_args()


#Fix formating for parameters after parsing to strings
args.working_directory = ''.join(args.working_directory)
args.input_bam = ''.join(args.input_bam)
args.ref_directory = ''.join(args.ref_directory)

file_name = args.input_bam

while '/' in file_name:
	file_name = file_name[file_name.find('/')+1:]

file_prefix = file_name[0:file_name.find('_Aligned')]

#print args.working_directory
#print args.input_bam
#print args.ref_directory
#print file_prefix


#Get environment username so the module can be user independent
USER = os.environ['USER']
PIPE_DIR = os.getcwd() # Get path to rna-pipeline directory so we can run the binary


#set SBATCH command line directives and submit sbatch rsem.sh script to SLURM queue
workdir = "--workdir=" + args.working_directory
output = "--output=/home/users/" + USER + "/out/rsem.%j.out"
error = "--error=/home/users/" + USER + "/errout/rsem.%j.err"
mail = "--mail-user=" + USER + "@stanford.edu" #assumes user's email is structured sherlock_username@stanford.edu

subprocess.call(['sbatch', workdir, output, error, mail, 'rsem.sh', '-b', args.input_bam, '-prefix', file_prefix, '-rDir', args.ref_directory, '-pipeDir', PIPE_DIR])









