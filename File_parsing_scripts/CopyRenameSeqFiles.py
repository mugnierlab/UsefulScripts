from sys import argv
import subprocess
import argparse
import time
import os

###Purpose: This script takes the raw seq output fastq files and renames them to what you make the os.chdir('FASTQ_BBEX17') say. 

filedata = open(argv[1],"r")

if not os.path.exists('Renamed_Files'):
    os.mkdir('Renamed_Files')



header_list = []
firstline = 0
sampledict = {}
filebasenames = []
                
for line in filedata.readlines():
        if firstline == 0:
            firstline = 1
        else:
            line = line.strip()
            linesplit = line.split(',')
            if linesplit[4] in sampledict:
                sampledict[linesplit[4]].append(linesplit[5])
            else:
                sampledict[linesplit[4]] = [linesplit[5]]

os.chdir('FASTQ_BBEX17')
for name in sampledict.iterkeys():
    subprocess.call('cat '+sampledict[name][0]+'_1.fastq '+sampledict[name][1]+'_1.fastq > ../Renamed_Files/'+name+'.fq', shell = True)    




