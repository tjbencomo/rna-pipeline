'''
Author: Tomas Bencomo
Script to fully automate pipeline operations. Takes as input a directory containing paired-end .gz FASTQ files, 
determines the pairs, and processes these samples through the pipeline.

Samples are processed in the following order:
1) STAR Alignment
2) RSEM and HTSeq-count quantification occurs in parallel

The script performs the following steps:
1) Determine the pairings for each sample within the directory
2) Launch a SLURM job using the STAR alignment program on each sample
3) Launch waiting RSEM and HTSeq-count programs for each sample once the sample's STAR job has finished
'''

import os
import subprocess
import argparse
import sys

parser = argparse.ArgumentParser(description = 'Launch complete pipeline')
parser.add_argument('-I', '--inputDirectory', metavar='--inputDirectory', type = str, nargs = 1, help = 'Directory containing FASTQ files to process', dest='directory')

args = parser.parse_args()

directory = ''.join(args.directory)

if directory == None:
	print "No directory specified. Exiting"
	sys.exit()


files = os.listdir(directory)
files = [file for file in files if '.fastq.gz' in file]

#Next - compute sample pairs by finding string between the _ and .fastq pattern
#Add the additional command line arguments for the other pipeline steps stuch as the star reference directory, rsem reference directory and other stuff
# write launch steps
# Turn all of these steps into functions
