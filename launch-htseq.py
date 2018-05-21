import os
import subprocess

BASE_DIR = '/scratch/PI/carilee/CACYPB-analysis/' 
all_files = os.listdir(BASE_DIR + 'align-files/')

files = []

for file in all_files:
	if "Aligned.sortedByCoord.out.bam" in file:
		files.append(file)



print len(files)
for file in files:
	print file



#launch htseq now

workdir = BASE_DIR + 'counts/'
basedir = BASE_DIR + '/align-files/'

for file in files:
        subprocess.call(['python', '/home/users/tbencomo/rna-pipeline/htseq-launch.py', '-wd', workdir, '-I', basedir + file])
