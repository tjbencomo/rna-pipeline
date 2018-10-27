#!/usr/bin/env python3

'''
Author:Tomas Bencomo
rna-pipeline.py acts as the main driver script for rna-pipeline.
From this script, the user specifies the files they wish to process
and pipeline settings. 

General Overview:
1) STAR Alignment
2) RSEM and HTSeq-Count Transcript Quantification
3) DESeq2 Differential Expression Analysis

Usage:
Options:
-I | --input : Specify input directory or input files. If specifying files, only 1 pair of files can be specified.
-outDir | --outputDirectory : Directory to store pipeline results files. Creates directory if doesn't already exist
-gDir | --genomeDirectory : Optional. Location of reference directory for STAR. If not used, uses directory from config.ini
-rDir | --rsemDirectory : Optional. Location of reference directory for RSEM. If not used, uses directory from config.ini
-countRef | --htseqCountReferenceFile : Optional. GTF file path for HTSEQ-Count. If not used, uses directory from config.ini
-outLog | --outputLogDirectory : Optional. Location to store output log file. If not used, uses directory from config.ini
-errorLog | --errorLogDirectory : Optional. Location to store error log file. If not used, uses directory from config.ini
-DESeq | --includeDESeqAnalysis : Optional. Boolean value whether or not to run differential expression analysis. Defaults to True if directory given as input.
'''

import argparse
import configparser
import io
import os
import difflib

def createParser():
    parser = argparse.ArgumentParser(description='Pipeline to process RNA-Seq data')
    parser.add_argument('-I', '--input', type=str, nargs='*', 
                        help='Specify input directory or input files.If specifying files, only 1 pair of files can be specified.',
                        dest='input')
    parser.add_argument('-gDir', '--genomeDirectory', type=str, nargs=1, 
                        help='Optional. Location of reference directory for STAR. If not used, uses directory from config.ini',
                        dest='genomeDirectory')
    parser.add_argument('-rDir', '--rsemDirectory', type=str, nargs=1,
                        help='Optional. Location of reference directory for RSEM. If not used, uses directory from config.ini',
                        dest='rsemDirectory')
    parser.add_argument('-countRef', '--htseqCountReferenceFile', type=str, nargs=1,
                        help='Optional. GTF file path for HTSEQ-Count. If not used, uses directory from config.ini',
                        dest='htseqReferenceFile')
    parser.add_argument('-outLog', '--outputLogDirectory', type=str, nargs=1,
                        help='Optional. Location to store output log file. If not used, uses directory from config.ini',
                        dest='outputLogDirectory')
    parser.add_argument('-errorLog', '--errorLogDirectory', type=str, nargs=1, 
                        help='Optional. Location to store error log file. If not used, uses directory from config.ini',
                        dest='errorLogDirectory')
    parser.add_argument('-DESeq', '--includeDESeqAnalysis', type=bool, nargs=1,
                        help='Optional. Boolean value whether or not to run differential expression analysis. Defaults to True.',
                        dest='includeDESeq')
    parser.add_argument('-outDir', '--outputDirectory', type=str, nargs=1,
                        help='Directory to store pipeline results files. Creates directory if doesnt exist',
                        dest='outputDirectory')
    
    args = parser.parse_args()

    return args

def parseArgs():
    
    args = createParser()

    # Argument handling - check input and fix if needed
    if len(args.input) == 1:
        pipelineInput = args.input[0]
        inputIsDirectory = True
    elif len(args.input) == 2:
        pipelineInput = args.input
        inputIsDirectory = False
        includeDESeq = False
    else:
        raise ValueError("More than 2 inputs given. rna-pipeline only accepts 1 input (directory) or 2 inputs (pair of FASTQ files)")

    if args.outputDirectory is None:
        raise ValueError("No output directory specified! Specify with the -outDir flag")
    else:
        outputDirectory = args.outputDirectory[0]
    

    configurationSettings = loadConfigData()

    if args.genomeDirectory is None:
        genomeDirectory = configurationSettings['rna-pipeline']['genomeDirectory']
    else:
        genomeDirectory = args.genomeDirectory[0]
    if args.rsemDirectory is None:
        rsemDirectory = configurationSettings['rna-pipeline']['rsemDirectory']
    else:
        rsemDirectory = args.rsemDirectory[0]
    if args.htseqReferenceFile is None:
        htseqReferenceFile = configurationSettings['rna-pipeline']['htseqGFFFFile']
    else:
        htseqReferenceFile = args.htseqReferenceFile[0]
    if args.outputLogDirectory is None:
        outputLogDirectory = configurationSettings['sherlock']['outputLogDirectory']
    else:
        outputLogDirectory = args.outputLogDirectory[0]
    if args.errorLogDirectory is None:
        errorLogDirectory = configurationSettings['sherlock']['errorLogDirectory']
    else:
        errorLogDirectory = args.errorLogDirectory[0]
    if args.includeDESeq is None:
        if inputIsDirectory:
            includeDESeq = True
    else:
        includeDESeq = args.includeDESeq[0]

    if os.path.isdir(outputDirectory) is False:
        os.makedirs(outputDirectory)
    
    arguments = {'pipelineInput' : pipelineInput, 'outputDirectory' : outputDirectory, 
                'genomeDirectory' : genomeDirectory, 'rsemDirectory' : rsemDirectory,
                'rsemDirectory' : rsemDirectory, 'htseqReferenceFile' : htseqReferenceFile,
                'outputLogDirectory' : outputLogDirectory, 'errorLogDirectory' : errorLogDirectory,
                'includeDESeq' : includeDESeq, 'inputIsDirectory' : inputIsDirectory}
    
    return arguments
    
def loadConfigData():
    userDirectory = os.environ['HOME']

    config = configparser.ConfigParser()
    config.read(os.path.join(userDirectory, 'config.ini'))

    return config

def getSamplePairs(directory):
    files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and '.f' in f]
    samples = []
    for f in files:
        files.remove(f)
        match = difflib.get_close_matches(f, files, n=1)[0]
        print(match)
        samples.append((f, match))
        files.remove(match)
    
    return samples
    
def processSample(sample):
    # initiate sbatch calls through python modules for each part of the pipeline
    pass

def launchPipeline(args):
    # single pair vs directory
    # list with tuples for each paired file to process
    
    if args['inputIsDirectory'] is True:
        # find all pairings
        samples = getSamplePairs(args['pipelineInput'])

        # for loop to processSamples() 
    else:
        # just one sample
        sample = [(args['pipelineInput'][0], args['pipelineInput'][1])]

        processSample(samples)

    


def main():
    args = parseArgs()
    launchPipeline(args)

if __name__ == '__main__':
    main()
