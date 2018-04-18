'''
Author: Tomas Bencomo
Script to execute star alignment. Acts as a wrapper script for star.sh, which actually calls the star executable.
Submits sbatch call, with command line option  SBATCH directives for --workdir, --output, and --error.
This script performs the following tasks:
Parse command line options
submit sbatch job for the star.sh star aligner
'''


import argparse
import subprocess
import os

#Build parser
parser = argparse.ArgumentParser(description='Launch star aligner')

parser.add_argument('-wd','--workingDirectory', metavar='--working-directory', type=str, nargs=1, help='file path for where star output files are stored', dest='working_directory')

parser.add_argument('-f1', '--fastq1', metavar='fastq1', type=str, nargs=1, help='first fastq input file', dest='fastq1')

parser.add_argument('-f2', '--fastq2', metavar='fastq2', type=str, nargs=1, help='second fastq input file', dest='fastq2')

parser.add_argument('-gDir', '--genomeDirectory', metavar='--genomeDirectory', type=str, nargs=1, help='reference files for star aligner', dest='genome_directory')

parser.add_argument('-prefix', '--outFilePrefix', metavar='--outFilePrefix', type=str, nargs=1, help='Prefix for output files for star aligner', dest='prefix')

args = parser.parse_args()



#If no output prefix specified, use filename as prefix and convert other parameters to strings
args.working_directory = ''.join(args.working_directory)
args.fastq1 = ''.join(args.fastq1)
args.fastq2 = ''.join(args.fastq2)
args.genome_directory = ''.join(args.genome_directory)

if args.prefix == None:
	prefix_string1 = ''.join(args.fastq1)
	prefix_string2 = ''.join(args.fastq2)
	while '/' in prefix_string1 or '/' in prefix_string2:
		prefix_string1 = prefix_string1[prefix_string1.find('/')+1:]
		prefix_string2 = prefix_string2[prefix_string2.find('/')+1:]
		#print(prefix_string);
	#print prefix_string

#	args.prefix = prefix_string[0:prefix_string.find('R1')]
	args.prefix = ""
	if len(prefix_string1) != len(prefix_string2):
		print 'ERROR ERROR PROBLEM WITH PREFIX NAMING. FASTQ FILES DIFF LENGTHS'
	else:
		for i in xrange(len(prefix_string1)):
			if prefix_string1[i] != prefix_string2[i]:
				break
			else:
				args.prefix += prefix_string1[i]
	print args.prefix
else:
	args.prefix = ''.join(args.prefix)


#We now need to sbatch star.sh with some additional batch directives
#--workdir
#--output
#--error
#All three options are dependent on the user - different users will have different filepaths - must make user indepenedent by modifying filepath according to user

USER = os.environ['USER']

workdir = "--workdir=" + args.working_directory
output = "--output=/home/users/" + USER + "/out/star-aligner.%j.out"
error = "--error=/home/users/" + USER + "/errout/star-aligner.%j.err"
mail = "--mail-user=" + USER + "@stanford.edu" #assumes user's email is structured sherlock_username@stanford.edu - add option later in command line parameters

#subprocess.call(['sbatch', 'star.sh', '-gDir', args.genome_directory, '-prefix', args.prefix, '-f1', args.fastq1, '-f2', args.fastq2, '-wd', '/scratch/users/tbencomo/RNA_seq/pipeline-tests'])

subprocess.call(['sbatch', workdir, output, error, mail, 'star.sh', '-gDir', args.genome_directory, '-prefix', args.prefix, '-f1', args.fastq1, '-f2', args.fastq2])


