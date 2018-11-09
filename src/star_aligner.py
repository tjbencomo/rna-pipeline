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
import rna_pipeline
import difflib

def createParser():
	parser = argparse.ArgumentParser(description='Pipeline module to align fastq files')
	parser.add_argument('-I', '--input', type=str, nargs='*', help='Specify input file(s). Should be gzipped', dest='input')
	parser.add_argument('-gDir', '--genomeDirectory', type=str, nargs=1, 
						help='Optional. Location of reference directory for STAR. If not used, uses directory from rna-pipeline-config.ini',
						dest='genomeDirectory')
	parser.add_argument('-outLog', '--outputLogDirectory', type=str, nargs=1,
                        help='Optional. Location to store output log file. If not used, uses directory from rna-pipeline-config.ini',
                        dest='outputLogDirectory')
	parser.add_argument('-errorLog', '--errorLogDirectory', type=str, nargs=1, 
                        help='Optional. Location to store error log file. If not used, uses directory from rna-pipeline-config.ini',
                        dest='errorLogDirectory')
	parser.add_argument('-outDir', '--outputDirectory', type=str, nargs=1,
                        help='Directory to store pipeline results files. Creates directory if doesnt exist',
                        dest='outputDirectory')
	
	args = parser.parse_args()
	return args

def parseArgs():
	args = createParser()

	configurationSettings = rna_pipeline.loadConfigData()

	# clean up arguments in dict
	arguments = {}

	if len(args.input) > 2:
		raise ValueError("More than 2 input files specified! Only 1 or 2 files can be set as input")
	else:
		arguments['pipelineInput'] = args.input
	
	if args.outputDirectory is None:
		arguments['outputDirectory'] = os.getcwd()
	else:
		arguments['outputDirectory'] = ''.join(args.outputDirectory)
	
	if args.genomeDirectory is None:
		arguments['genomeDirectory'] = configurationSettings['rna-pipeline']['genomeDirectory']
	else:
		arguments['genomeDirectory'] = args.genomeDirectory[0]
	
	if args.outputLogDirectory is None:
		arguments['outputLogDirectory'] = configurationSettings['sherlock']['outputLogDirectory']
	else:
		arguments['outputLogDirectory'] = args.outputLogDirectory[0]
	
	if args.errorLogDirectory is None:
		arguments['errorLogDirectory'] = configurationSettings['sherlock']['errorLogDirectory']
	else:
		arguments['errorLogDirectory'] = args.errorLogDirectory[0]
	
	if os.path.isdir(arguments['outputDirectory']) is False:
		os.makedirs(arguments['outputDirectory'])

	return arguments

def getPrefix(fileNames):
	'''
	Expects 2 element list containing filenames
	Returns the common prefix between the filenames - if full paths are given, truncates path to return only file name
	'''
	sequence = difflib.SequenceMatcher(None, fileNames[0], fileNames[1])
	results = sequence.find_longest_match(0, len(fileNames[0]), 0, len(fileNames[1]))
	return os.path.basename(fileNames[0][results[0] : results[0] + results[2]])

def runSTAR(args):
	'''
	Required args keys:
		-outputLogDirectory
		-errorLogDirectory
		-pipelineInput
		-outputDirectory
		-genomeDirectory
	Returns SLURM ID for submitted job
	'''
	
	outputLogFile = os.path.join(args['outputLogDirectory'], 'star_aligner.%j.out')
	errorLogFile = os.path.join(args['errorLogDirectory'], 'star_aligner.%j.err')
	# memory = '80000' # This is in MB
	memory = '1000'
	# time = '0-07:00:00'
	time = '0-00:05:00'
	nodes = '1'
	cpuPerNode = '$SLURM_CPUS_ON_NODE'
	limitBAMsortRAM = '20000000000' # in KB I think - check STAR documentation
	
	if len(args['pipelineInput']) != 1:
		# 2 files for paired end sequencing
		outFilePrefix = os.path.join(args['outputDirectory'], getPrefix(args['pipelineInput']))
		starInput = ' '.join(args['pipelineInput'])
	else:
		# single file
		print(args['pipelineInput'])
		starInput = ''.join(args['pipelineInput'])
		outFilePrefix = os.path.join(args['outputDirectory'], os.path.basename(starInput[0:starInput.index('.fast')]))
		
	
	starCommand = 'STAR --genomeDir {} --readFilesIn {} --readFilesCommand zcat --runThreadN {}' \
					'--genomeLoad NoSharedMemory --outFilterMultimapNmax 20 --alignSJoverhangMin 8' \
					'--alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04' \
					'--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000' \
					'--outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate' \
					'--outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD' \
					'--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --sjdbScore 1' \
					'--limitBAMsortRAM {} --outFileNamePrefix {}'.format(args['genomeDirectory'], starInput, cpuPerNode, limitBAMsortRAM, os.path.join(args['outputDirectory'], outFilePrefix))
	sbatchCommand = 'sbatch --job-name=star_aligner --output={} --error={} --nodes={} --mem={} --time={} --wrap="{}"'.format(outputLogFile, errorLogFile, nodes, memory, time, starCommand)

	# return sbatchCommand

	output = subprocess.run(sbatchCommand, shell=True, stdout=subprocess.PIPE)
	slurmJobID = output.stdout.split()[-1].decode("utf-8")
	
	return slurmJobID
	
def test():
	args = {'outputDirectory' : 'outputdir', 'pipelineInput' : ['file1_2.fq', 'file1_1.fq'], 'genomeDirectory' : 'somegenomepath',
			'outputLogDirectory' : 'log_path'}
	print(runSTAR(args))

def main():
	# test()
	args = parseArgs()
	results = runSTAR(args)
	print(results)

if __name__ == '__main__':
	main()


# #Build parser
# parser = argparse.ArgumentParser(description='Launch star aligner')

# parser.add_argument('-wd','--workingDirectory', metavar='--working-directory', type=str, nargs=1, help='file path for where star output files are stored', dest='working_directory')

# parser.add_argument('-f1', '--fastq1', metavar='fastq1', type=str, nargs=1, help='first fastq input file', dest='fastq1')

# parser.add_argument('-f2', '--fastq2', metavar='fastq2', type=str, nargs=1, help='second fastq input file', dest='fastq2')

# parser.add_argument('-gDir', '--genomeDirectory', metavar='--genomeDirectory', type=str, nargs=1, help='reference files for star aligner', dest='genome_directory')

# parser.add_argument('-prefix', '--outFilePrefix', metavar='--outFilePrefix', type=str, nargs=1, help='Prefix for output files for star aligner', dest='prefix')

# args = parser.parse_args()



# #If no output prefix specified, use filename as prefix and convert other parameters to strings
# args.working_directory = ''.join(args.working_directory)
# args.fastq1 = ''.join(args.fastq1)
# args.fastq2 = ''.join(args.fastq2)
# args.genome_directory = ''.join(args.genome_directory)
# '''
# if args.prefix == None:
# 	prefix_string1 = ''.join(args.fastq1)
# 	prefix_string2 = ''.join(args.fastq2)
# 	while '/' in prefix_string1 or '/' in prefix_string2:
# 		prefix_string1 = prefix_string1[prefix_string1.find('/')+1:]
# 		prefix_string2 = prefix_string2[prefix_string2.find('/')+1:]
# 		#print(prefix_string);
# 	#print prefix_string

# #	args.prefix = prefix_string[0:prefix_string.find('R1')]
# 	args.prefix = ""
# 	if len(prefix_string1) != len(prefix_string2):
# 		print 'ERROR ERROR PROBLEM WITH PREFIX NAMING. FASTQ FILES DIFF LENGTHS'
# 	else:
# 		for i in xrange(len(prefix_string1)):
# 			if prefix_string1[i] != prefix_string2[i]:
# 				break
# 			else:
# 				args.prefix += prefix_string1[i]
# 	print args.prefix
# else:
# 	args.prefix = ''.join(args.prefix)
# '''


# if args.prefix == None:
# 	#print("No prefix given, computing prefix")
# 	file_name = args.fastq1[args.fastq1.rfind('/')+1 : ] 
# 	file_extension_pos = file_name.find('.fastq.gz')
# 	paired_end_indicator_pos = file_name.rfind('_', 0, file_extension_pos)
# 	args.prefix = file_name[0:paired_end_indicator_pos]
# 	print args.prefix
# 	#print("Prefix should have been printed above")
# 	args.prefix += '_'
# else:
# 	args.prefix = ''.join(args.prefix)

# #We now need to sbatch star.sh with some additional batch directives
# #--workdir
# #--output
# #--error
# #All three options are dependent on the user - different users will have different filepaths - must make user indepenedent by modifying filepath according to user

# USER = os.environ['USER']
# PIPE_DIR = os.getcwd()

# workdir = "--workdir=" + args.working_directory
# output = "--output=/home/users/" + USER + "/out/star-aligner.%j.out"
# error = "--error=/home/users/" + USER + "/errout/star-aligner.%j.err"
# mail = "--mail-user=" + USER + "@stanford.edu" #assumes user's email is structured sherlock_username@stanford.edu - add option later in command line parameters

# #subprocess.call(['sbatch', 'star.sh', '-gDir', args.genome_directory, '-prefix', args.prefix, '-f1', args.fastq1, '-f2', args.fastq2, '-wd', '/scratch/users/tbencomo/RNA_seq/pipeline-tests'])

# subprocess.call(['sbatch', workdir, output, error, mail, 'star.sh', '-gDir', args.genome_directory, '-prefix', args.prefix, '-f1', args.fastq1, '-f2', args.fastq2, '-pipeDir', PIPE_DIR])


