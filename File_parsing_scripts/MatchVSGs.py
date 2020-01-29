import glob
import os, sys
import subprocess

## Purpose:This script takes the clustered VSG file and matches the representative VSG Trinity name with the cluster name for each cluse. The output is a text file that has [trinity name] [tab] [Cluster name].
## This lets you use the output to rename VSG Trinity names to the cluster names for easier use in R and for matching VSGs and clusters between samples run throught the pipeline seperatly. 


currDir = os.getcwd()

cluster_dict = {}
with open('BBEX17VSGSeq_orf_VSGs_merged.fa.clstr','r') as infile_split:
    with open('concat_RESULTS.txt','r') as outfile:
        with open('cluster_reference_table.txt','w+') as reference:
            for line in infile_split:
                if line.startswith('>'):
                    clean_line = line.split('>')
                    cluster_dict[clean_line[1]]=[]
                    lastline = clean_line[1]
                else:
                    newline = line.split()
                    if len(newline) >1:
                        final = newline[2].split('>')
                        final = final[1].split('...')
                        cluster_dict[lastline].append(final[0])
            
            for key, value in cluster_dict.iteritems():
                for item in value:
                    reference.write(str(item+'\t'+key))
                

infile_split.close()
outfile.close()