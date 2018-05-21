import os
import subprocess



BASE_DIR = '/scratch/PI/carilee/CACYPB-analysis/'

files = os.listdir(BASE_DIR + 'fastq/')


pairs = []


def findPair(file, files):
        pairs = []
        prefix = file[:file.find('_')]
        for elem in files:
                elem_prefix = elem[:elem.find('_')]
                if(elem_prefix == prefix) and elem != file:
                        return (file, elem)
def findPair2(file, files):
        pairs = []
        start = file.find('_')
        prefix = file[:file.find('_', start+1)]
        for elem in files:
                begin = elem.find('_')
                elem_prefix = elem[:elem.find('_', begin+1)]
                if(elem_prefix == prefix) and elem != file:
                        return (file, elem)

def equals(x, y):
        a = set()
        b = set()
        for i in xrange(len(x)):
                a.add(x[i])
                b.add(y[i])
        return a == b


for file in files:
        pair = findPair(file, files)
        shouldAdd = True
        for elem in pairs:
                if equals(elem, pair):
                        shouldAdd = False
                        break
        if shouldAdd:
                pairs.append(pair)

print len(pairs)
for pair in pairs:
	print pair


#launch star now

workdir = BASE_DIR + 'align-files/'
basedir = BASE_DIR + 'fastq/'

for pair in pairs:
	subprocess.call(['python', '/home/users/tbencomo/rna-pipeline/star-aligner.py', '-wd', workdir, '-f1', basedir + pair[0], '-f2', basedir + pair[1], '-gDir', '/scratch/users/tbencomo/RNA_seq/refs/out'])


