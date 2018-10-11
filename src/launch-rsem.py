import os
import subprocess


BASE_DIRECTORY = '/scratch/PI/carilee/CACYPB-analysis/'

all_files = os.listdir(BASE_DIRECTORY + 'align-files/')

files = []

for file in all_files:
        if "Aligned.toTranscriptome.out.bam" in file:
                files.append(file)



print len(files)
for file in files:
        print file



#launch rsem now

workdir = BASE_DIRECTORY + 'rsem/'
basedir = BASE_DIRECTORY + 'align-files/'

for file in files:
	subprocess.call(['python', '/home/users/tbencomo/rna-pipeline/rsem-calculate.py', '-wd', workdir, '-I', basedir + file, '-rDir', '/scratch/users/tbencomo/RNA_seq/refs/out/rsem'])
