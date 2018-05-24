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

#parseArgs
#getSamplePairs
#launchPipeline
#main


def findParentDirectory(directory):
	if directory[len(directory)-1] == '/':
		splitter = directory.rfind('/', 0, len(directory) - 1)
	else:
		splitter = directory.rfind('/', 0, len(directory))
	parent_directory = directory[0:splitter] + '/'
	
	return parent_directory

def extractDirectoryName(directory):
	if directory[len(directory)-1] == '/':
		splitter = directory.rfind('/', 0, len(directory) - 1)
		directoryName = directory[splitter+1 : len(directory)-1]
	else:
		splitter = directory.rfind('/', 0, len(directory))
		directoryName = directory[splitter+1: ]
	return directoryName

def parseArgs():
	parser = argparse.ArgumentParser(description = 'Launch complete pipeline')
	parser.add_argument('-I', '--inputDirectory', metavar='--inputDirectory', type = str, nargs = 1, 
				help = 'Directory containing FASTQ files to process', dest='directory')
	parser.add_argument('-gDir', '--genomeDirectory', metavar='--genomeDirectory', type=str, nargs=1, 
				help = 'Directory containing reference genome files for STAR', dest= 'genomeDirectory')
	parser.add_argument('-rDir', '--rsemReferenceDirectory', metavar='--rsemReferenceDirectory', type=str, nargs=1, 
				help = 'Directory containing RSEM reference files', dest='rsemDirectory')
	parser.add_argument('-output', '--outputDirectory', metavar='--outputDirectory', type=str, nargs=1, dest='outputDirectory')

	args = parser.parse_args()

	directory = ''.join(args.directory)
	genomeDirectory = ''.join(args.genomeDirectory)
	rsemDirectory = ''.join(args.rsemDirectory)

	if args.outputDirectory == None:
		outputDirectory = findParentDirectory(''.join(args.directory)) + extractDirectoryName(directory) + '-output'
	else:
		outputDirectory = ''.join(args.outputDirectory)


	return (directory, genomeDirectory, rsemDirectory, outputDirectory) 


def getSamplePairs(directory):
	if directory == None:
		print "No directory specified. Exiting"
		sys.exit(1)

	files = os.listdir(directory)
	files = [file for file in files if '.fastq.gz' in file]
	
	samples = {}

	for file in files:
		file_extension_pos = file.find('.fastq.gz')
		paired_end_indicator_pos = file.rfind('_', 0, file_extension_pos)
		pair_version = file[paired_end_indicator_pos:file_extension_pos]
		prefix = file[0:paired_end_indicator_pos]
		
		if prefix in samples:
			samples[prefix].append(file)
		else:
			samples[prefix] = [file]
	
	for sample in samples:
		if len(samples[sample]) != 2:
			print("Not a pair of files. Either too less or too many files to be considered a pair. Something went wrong with pair processing. Sample " + prefix + " has been removed from analysis")
			samples.pop(sample, None)

	return samples


def launchPipeline(genome_directory, rsem_directory, output_directory, samples):
	PIPE_DIR = os.getcwd()
	for sample in samples:
		paired_files = samples[sample]
		output = subprocess.check_output(['python', PIPE_DIR + '/star-aligner.py', '-wd', output_directory, 
						'-f1', paired_files[0], '-f2', paired_files[1], '-gDir', genome_directory])
		phrase = 'Submitted batch job'
		index = output.find(phrase)
		end = index + len(phrase)
		jobID = output[end + 1:output.find('\n')]

def main():
	directory, genomeDirectory, rsemDirectory, outputDirectory = parseArgs()
	samples = getSamplePairs(directory)
	launchPipeline(genomeDirectory, rsemDirectory, outputDirectory, samples)
	print("All jobs launched")


if __name__ == "__main__": main()


# finish creating a directory-pipeline-output option if no output directory inputted
# write launch steps
# Turn all of these steps into functions
