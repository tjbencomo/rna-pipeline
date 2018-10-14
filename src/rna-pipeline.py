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
-outDir | --outputDirectory : Directory to store pipeline results files
-gDir | --genomeDirectory : Optional. Location of reference directory for STAR. If not used, uses directory from config.ini
-rDir | --rsemDirectory : Optional. Location of reference directory for RSEM. If not used, uses directory from config.ini
-countRef | --htseqCountReferenceFile : Optional. GTF file path for HTSEQ-Count. If not used, uses directory from config.ini
-outLog | --outputLogDirectory : Optional. Location to store output log file. If not used, uses directory from config.ini
-errorLog | --errorLogDirectory : Optional. Location to store error log file. If not used, uses directory from config.ini
-DESeq | --includeDESeqAnalysis : Optional. Boolean value whether or not to run differential expression analysis. Defaults to True.
'''

import argparse

def parseArgs():
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
                        dest='outputDirectory')
    parser.add_argument('-errorLog', '--errorLogDirectory', type=str, nargs=1, 
                        help='Optional. Location to store error log file. If not used, uses directory from config.ini',
                        dest='errorDirectory')
    parser.add_argument('-DESeq', '--includeDESeqAnalysis', type=bool, nargs=1,
                        help='Optional. Boolean value whether or not to run differential expression analysis. Defaults to True.',
                        dest='includeDESeq')
    
    args = parser.parse_args()

    if len(args.input) == 1:
        pipelineInput = args.input[0]
    elif len(args.input) == 2:
        pipelineInput = args.input
    else:
        raise ValueError("More than 2 inputs given. rna-pipeline only accepts 1 input (directory) or 2 inputs (pair of FASTQ files)")

    if args.genomeDirectory is None:
            


def main():
    parseArgs()

if __name__ == '__main__':
    main()