from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio import SeqIO 
from Bio import Seq
import subprocess
import argparse
import time
import os

def makeFilesList(filesList): #parses file containing information on files/samples
	file = open(filesList,"r")
	header_list = []
	firstline = 0
	sampledict = {}
	filebasenames = []
			
	for line in file.readlines():
		if firstline == 0:
			line = line.strip('\n')
			header_list = line.split('\t')[0:]
			firstline = 1
		else:
			line = line.strip('\n')
			linesplit = line.split('\t')
			sampledict['_'.join(linesplit[0].split('.')[:-1])] = line.split('\t')[0:] 
			filebasenames.append('_'.join(linesplit[0].split('.')[:-1]))
	return (filebasenames, header_list, sampledict)


parser = argparse.ArgumentParser()

# arguments about input files and naming
parser.add_argument('-i',help='text file with .fastq names of sequence files to be put through the pipeline, same for entering any step of pipeline, fq', action="store", dest="i", default='')
parser.add_argument('-header', help='name for output folder, overrides default Y-m-d-H_M format ', action="store", dest='head', default='')

# trimming setting
parser.add_argument('-trim', help='perform quality trimming, default is to skip this step', action='store_false', dest='trim')
parser.add_argument('-g', help='trim_galore, stringency for trim galore', action ="store", dest = "g", default="3") 
parser.add_argument('-trimlen', help='cutadapt/trim_galore, minimum length of read for trimming', action ="store", dest = "trimlen", default="50") 

# trinity settings
parser.add_argument('-minp', help='Trinity, minimum protein length you are filtering for', action ="store", type=int, dest = "minp", default=300) 
parser.add_argument('-mem', help='Trinity, max memory allocation for trinity, G', action ="store", dest = "mem", default="10") 
parser.add_argument('-cpu', help='Trinity/cd-hit-est/MULTo/bowtie, number of processors to use', action="store", dest='cpu', default='2')

# Blast settings
parser.add_argument('-vsgdb', help='BLAST, name of the vsg database to blast against', action ="store", dest = "vsgdb", default="concatAnTattb427")
parser.add_argument('-NoNonVSGdb', help='BLAST, adding this flag means you dont want to BLAST ORFS against the non-VSG database' , dest='NoNonVSGdb', action='store_true')

# cd-hit-est parameters
parser.add_argument('-sit', help='cd-hit-est, sequence identiy threshold - how much the alignment has to match, percentage. value is 0.0 through 1.0 ', action ="store", dest = "sit", default=".98")
parser.add_argument('-t', help='cd-hit-est, number of threads for cd-hit-est to use, 0 default(all CPU will be used)', action ="store", dest = "t", default="0") #check this

# MULTO settings
parser.add_argument('-p', help='MULTo, path to MULTo1.0 folder. default is current directory, please dont use "~/", python doesnt like this in the path', action="store", dest='p', default='') # default assumes MULTo is in your home dir
parser.add_argument('-v', help='number of mismatches allowed for bowtie, default 2', action="store", dest='v', default='2')
parser.add_argument('-reuseMulto', help='MULTo, name of the multo files to be reused, typically same as header, default will make new MULTo files. Final output won\'t provide information on VSG', action="store", dest='rmulto', default='')

arguments = parser.parse_args()

sample_info =  makeFilesList(arguments.i)
filebasenames = sample_info[0]
header_list = sample_info[1]
sampledict = sample_info[2]
min_pro_len = arguments.minp
trim = arguments.trim
seqIdenThresh = arguments.sit
numCPU = arguments.cpu
max_memory_cdhit = 1000*int(arguments.mem)
path = arguments.p
numMismatch = arguments.v
rmulto = arguments.rmulto
trimlen = arguments.trimlen
stringency = arguments.g
currDir = os.getcwd()


#checking for all required software
#if path != '':
#	fullmultopath = path + '/MULTo1.0'
#else:
#	fullmultopath = currDir + '/MULTo1.0'

#try:
	#subprocess.check_call(["bowtie","--version"])
	#subprocess.check_call(["blastn","-version"])
#	subprocess.check_call(["trim_galore","--version"]) 
	#subprocess.check_call(["cutadapt","--version"])
	#subprocess.check_call(["python", str(fullmultopath+'/src/MULTo1.0.py'),"-h"])
	#subprocess.check_call(["python", str(fullmultopath+'/src/rpkmforgenes.py'),"-h"])
#except:
#	raise SystemExit, "Problem opening required software. Is everything in your path?"
#try:	
#	subprocess.check_output(["cd-hit","--version"])
#except subprocess.CalledProcessError as cpe:
#		out = cpe.output
#		if out[9:15] != 'CD-HIT':
#			raise SystemExit, "Problem opening required software. Is cd-hit in your path?"	
#except OSError:
#	raise SystemExit, "Problem opening required software. Is cd-hit in your path?"
#try:	
#	subprocess.check_output(["Trinity","--version"])
#except subprocess.CalledProcessError as cpe:
#		out = cpe.output
#		print out[0:7]
#		if out[0:7] != 'Trinity':
#			raise SystemExit, "Problem opening required software. Is Trinity in your path?"	
#except OSError:
#	raise SystemExit, "Problem opening required software. Is Trinity in your path?"	

header = time.strftime("%Y-%m-%d-%H_%M")  # Year/Month/Day-Hour:Minute , names the folder for output files

if header != "":
	header = str(arguments.head)
	
if not os.path.exists(header):
		os.makedirs(header) # creates the folder
		
if not os.path.exists(header+"/StandardError"):
	os.makedirs(header+"/StandardError")

scoredict = False

if trim == True:
	for file in filebasenames:
		if os.path.exists(str(file+'.fq')): 
			filepath = str(str(file) +'.fq')
		elif os.path.exists(str(file+'.fastq')):
			filepath = str(str(file) +'.fastq')
		stderr_tg = " 2> " + header + "/StandardError/trim_galore-"+file+".txt"
		stderr_ca = " 2> " + header + "/StandardError/cutadapt-"+file+".txt"
		subprocess.call(['trim_galore --stringency '+str(stringency)+' --length '+trimlen+' --dont_gzip --output_dir ' + header + "/ " +str(filepath)+stderr_tg], shell=True) # trim off sequenceing adapters
		subprocess.call(['cutadapt -m '+trimlen+' -b ATTTAGGTGACACTATAG -b CTATAGTGTCACCTAAAT '+ header + "/"+str(file)+'_trimmed.fq > '+ header + "/"+str(file)+'_trimmed2.fq'+stderr_ca], shell=True) # trim off SP6 sequences (from VSG PCR step)
		subprocess.call(['rm '+ header + "/"+str(file)+'_trimmed.fq'], shell=True) # removes intermediate trimmed file 