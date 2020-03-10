from sys import argv
import subprocess
import argparse
import time
import os

# Purpose: This script renames the fastq sequencing files to the sample names using the CSV file that 
# contains the indexes and the sample names (this file is specific to our GRCF seqeucing core and you get this file back with the FASTQs).

# Keep the directory hierarchy the same as when the FASTQ files are downloaded.
# The input file is the .CSV file that contains the indexes and the corresponding sample names. (argv[1])
# The renamed files are placed in a new directory named 'Renamed_Files'.
# note: The sequencing FASTQs need to be in a directory called 'FASTQ' within your working directory, this is the name the core gives us. 

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

os.chdir('FASTQ')
for name in sampledict.keys():
    if sampledict[name] == 1:
        subprocess.call('cat '+sampledict[name][0]+'_1.fastq ', shell = True)
    else:   
        subprocess.call('cat '+sampledict[name][0]+'_1.fastq '+sampledict[name][1]+'_1.fastq > ../Renamed_Files/'+name+'.fq', shell = True)    




