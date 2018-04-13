'''
Author: Tomas Bencomo

Script to execute HTSeq-count software package to count the number of RNA reads per sequence. This script acts as a wrapper function for the htseq.sh script, which actually calls the HTSeq-count software. 

This script:
Parses the command line options
Submits sbatch job for htseq.sh script

Command line parameters:
input bam file
working directory

Internally this script
Determines the prefix name based on the input file
'''


import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description='Launch htseq-count')

parser.add_argument('-I', '--input-bam-file', metavar='--input-bam-file', type=str, nargs=1, help='Sorted by coordinate bam file to count reads', dest='input_bam')

parser.add_argument('-wd', '--working-directory', metavar='--working-directory', type=str, nargs=1, help='location to store output files', dest='working_directory')

args = parser.parse_args()


#Fix formating to convert to strings
args.working_directory = ''.join(args.working_directory)
args.input_bam = ''.join(args.input_bam)

file_name = args.input_bam

while '/' in file_name:
	file_name = file_name[file_name.find('/')+1:]

file_prefix = file_name[0:file_name.find('_Aligned')]

print args.working_directory
print args.input_bam
print file_prefix

#Get environment username so the module can be user independent
USER = os.environ['USER']
PIPE_DIR = os.getcwd()


print "pipeline dir : " + PIPE_DIR


#set SBATCH command line directives and submit sbatch htseq.sh script to SLURM queue
workdir = "--workdir=" + args.working_directory
output = "--output=/home/users/" + USER + "/out/htseq.%j.out"
error = "--error=/home/users/" + USER + "/errout/htseq.%j.err"
mail = "--mail-user=" + USER + "@stanford.edu" #assumes user's email is structured sherlock_username@stanford.edu

subprocess.call(['sbatch', workdir, output, error, mail, 'htseq.sh', '-b', args.input_bam, '-prefix', file_prefix, '-pipeDir', PIPE_DIR])


