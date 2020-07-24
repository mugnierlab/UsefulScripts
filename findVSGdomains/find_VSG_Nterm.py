from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd
import subprocess
import linecache
import argparse
import sys
import os


def iterate_trim(sequence, list_of_seq_variants):
	"""
	Iteratively trims 1 amino acid from the 3' end of a sequence. Stops when 1 amino acid remains.
	Returns a list of trimmed sequence variants.
	"""

	if len(sequence.seq) > 1: # If sequence is longer than the minimum of 1 amino acid.
		single_trim = SeqRecord(
		seq = sequence.seq[:-1], # Remove an amino acid from the 3' end.
		id = list_of_seq_variants[0].id.split("Ntrimmed0")[0] + "Ntrimmed" + str(len(list_of_seq_variants[0]) - len(sequence.seq[:-1])), # Append the total number of amino acids removed to the sequence ID.
		description = "",
		)
		list_of_seq_variants.append(single_trim) # Append trimmed sequence to a list.
		iterate_trim(list_of_seq_variants[-1], list_of_seq_variants) # The last sequence appended to the list of sequence variants is trimmed.
	return list_of_seq_variants

def hmmscan_type(hmmscan_infile, path_hmm):
	"""
	Runs HMMER hmmscan.
	hmmscan searches each sequence variant against the VSG N-terminal Type A and B HMM profile databases.
	"""

	hmmscan_out_base = os.path.splitext(hmmscan_infile)[0] # FASTA file containing trimmed variants of the input sequence(s) w/o extension.

	# TypeA N-terminal Domain HMM Profile
	TypeA_file = hmmscan_out_base + "_TypeA.out" # Name for hmmscan table output file (TypeA).
	subprocess.call(['hmmscan --noali -o '+hmmscan_out_base+'_TypeA_hmmscan.txt --tblout '+TypeA_file+' '+path_hmm+'VSG-N-TypeA.hmm '+hmmscan_infile+''], shell=True)
	# TypeB N-terminal Domain HMM Profile
	TypeB_file = hmmscan_out_base + "_TypeB.out" # Name for hmmscan table output file (TypeB).
	subprocess.call(['hmmscan --noali -o '+hmmscan_out_base+'_TypeB_hmmscan.txt --tblout '+TypeB_file+' '+path_hmm+'VSG-N-TypeB.hmm '+hmmscan_infile+''], shell=True)

	hmmscan_out_files = []

	for file in os.listdir(os.getcwd()):
		if file.endswith("_TypeA.out") or file.endswith("_TypeB.out"):
			hmmscan_out_files.append(file)
	if len(hmmscan_out_files) != 2:
		sys.exit("HMMER hmmscan output files were not created or more output files were found than created.")

def nterm_type_ID(hmmscan_out, all_original_seqIDs):
	"""
	Identifies the most probable N-terminal domain sequence of the input sequence(s).
	Returns a list of the most probable N-terminal domain sequence ID(s).
	"""

	row_values = []
	seq_IDs = []
	seq_scores = []

	with open(hmmscan_out, "r") as hmmscan_out_table:
		for line in hmmscan_out_table.readlines():
			filler = "#" # Lines without sequence data contain "#".
			if filler not in line:
				line = line.strip("\n")
				row_values = line.split()[0:] # Remove white-space delimiters and place values in list.
				seq_IDs.append(row_values[2]) # Append sequence ID (with "_NtrimmedX") to list.
				seq_scores.append(float(row_values[5])) # Append corresponding sequential bit score to list.

	seqID_score_dict = dict(zip(seq_IDs, seq_scores)) # Dictionary - {sequence ID: sequential bit score}.

	most_probable_variant_seqIDs = []

	for original_seqID in all_original_seqIDs: # For each sequence ID in original FASTA file (without "_NtrimmedX").
		if any(original_seqID + "_" in key for key in seqID_score_dict):
			# Make subset dictionary that contains trimmed variants of an input sequence.
			one_seq_dict = {key: value for key, value in seqID_score_dict.items() if original_seqID + "_" in key}
			one_seq_maximum_score = max(one_seq_dict.values()) # Maximum score given to a trimmed variant in subset dictionary.
			variants_with_max_score = []
			for key, value in one_seq_dict.items():
				if value == one_seq_maximum_score:
					variants_with_max_score.append(key) # Append the sequence ID of the maximum scoring variant(s) to list.
					# Make subset dictionary that contains highest scoring trimmed variant(s) of an input sequence.
					max_score_variant_length_dict = {}
					for max_score_variant in variants_with_max_score:
							trimmed = max_score_variant.split("Ntrimmed")[1] # Number of amino acids trimmed to obtain variant.
							max_score_variant_length_dict.update({max_score_variant: trimmed}) # Dictionary - {highest scoring variant sequence ID: number of amino acids trimmed}.
			# The highest scoring variant with the smallest number of trimmed amino acids is the most probable sequence.
			most_probable_variant_seqIDs.append(min(max_score_variant_length_dict, key = max_score_variant_length_dict.get)) # List of sequence ID(s) of the most probable N-terminal domain sequence(s).
	return most_probable_variant_seqIDs

def hmmscan_subtype(hmmscan_infile, path_hmm):
	"""
	Runs HMMER hmmscan.
	hmmscan searches each sequence variant against the VSG N-terminal Type A(1-3) and B(1-2) HMM profile databases.
	"""

	hmmscan_out_base = os.path.splitext(hmmscan_infile)[0] # FASTA file containing the identified N-terminal domain(s) w/o extension.

	# TypeA1-3 and TypeB1-2 N-terminal Domain HMM Profiles
	ABsubtype_file = hmmscan_out_base + "_TypeABsubtype.out" # Name for hmmscan table output file (TypeA1-3/TypeB1-2).
	subprocess.call(['hmmscan --noali -o '+hmmscan_out_base+'_TypeABsubtype_hmmscan.txt --tblout '+ABsubtype_file+' '+path_hmm+'VSG-N-TypeSubtype.hmm '+hmmscan_infile+''], shell=True)

	hmmscan_out_files = []

	for file in os.listdir(os.getcwd()):
		if file.endswith("_TypeABsubtype.out"):
			hmmscan_out_files.append(file)
	if len(hmmscan_out_files) != 1:
		sys.exit("HMMER hmmscan output file was not created or more output files were found than created.")

def nterm_subtype_ID(hmmscan_out, all_original_seqIDs):
	"""
	Identifies the most probable N-terminal domain subtype of the N-terminal TypeA and/or TypeB domain sequence(s).
	Returns a dictionary containing the N-terminal domain sequence ID(s) and the corresponding type-subtype.
	"""

	row_values = []
	seq_profile = []
	seq_IDs = []
	seq_scores = []

	with open(hmmscan_out, "r") as hmmscan_out_table:
		for line in hmmscan_out_table.readlines():
			filler = "#" # Lines without sequence data contain "#".
			if filler not in line:
				line = line.strip("\n")
				row_values = line.split()[0:] # Remove white-space delimiters and place values in list.
				seq_profile.append(str(row_values[0])[-2:]) # Appends hit profile type-subtype to list.
				seq_IDs.append(row_values[2]) # Append N-terminal domain sequence ID (with "_NtrimmedX") to list.
				seq_scores.append(float(row_values[5])) # Appends corresponding sequential bit score to list.

	seq_profile_seqID_score_list = list(zip(seq_profile, seq_IDs, seq_scores))
	seq_profile_seqID_score_df = pd.DataFrame(seq_profile_seqID_score_list, columns = ["seq_profile", "seqID", "seq_score"]) # DataFrame with columns for hit profile type-subtype, N-terminal domain sequence ID, and sequential bit score.

	most_probable_variant_seqIDs = []
	most_probable_variant_type = []

	for original_seqID in all_original_seqIDs: # For each sequence ID in original FASTA file (without "_NtrimmedX").
		if any(original_seqID + "_" in id for id in seq_profile_seqID_score_df["seqID"]):
			# Make subset DataFrame that contains all N-terminal domain sequence type-subtype hits of an input sequence.
			one_seq_df = seq_profile_seqID_score_df.loc[seq_profile_seqID_score_df["seqID"].str.contains(original_seqID + "_")]
			one_seq_df = one_seq_df.drop_duplicates(subset = ["seqID"], keep = "first") # Remove duplicate type-subtype hits, keep first occurence (most probable).
			most_probable_variant_seqIDs.append(one_seq_df.iloc[0]["seqID"]) # Append N-terminal domain sequence ID to list.
			most_probable_variant_type.append(one_seq_df.iloc[0]["seq_profile"]) # Append N-terminal domain type-subtype to list.
		most_probable_dict = dict(zip(most_probable_variant_seqIDs, most_probable_variant_type)) # Dictionary - {N-terminal domain sequence ID: type-subtype}.
	return most_probable_dict

def summary_csv(original_seq_file, nterms_file):
	"""Creates a CSV file containing information about the input VSG sequence(s) and the identified N-terminal domain of said sequence(s)."""

	original_seqID = []
	original_seqlength = []
	nterm_seqID = []
	nterm_seqlength = []
	nterm_type = []
	nterm_typesubtype = []
	nterm_seq = []

	for nterm_sequence in SeqIO.parse(nterms_file, "fasta"):
		for original_sequence in SeqIO.parse(original_seq_file, "fasta"):
			if original_sequence.id + "_" in nterm_sequence.id:
				original_seqID.append(original_sequence.id)

				trimmed = int(nterm_sequence.id.split("Ntrimmed")[1].split("_Type")[0])
				original_seqlength.append(len(nterm_sequence.seq) + trimmed)

				nterm_seqID.append(nterm_sequence.id)
				nterm_seqlength.append(len(nterm_sequence.seq))
				nterm_type.append(nterm_sequence.id[-2:-1])
				nterm_typesubtype.append(nterm_sequence.id[-2:])
				nterm_seq.append(nterm_sequence.seq)

	summary_df_file = os.path.splitext(original_seq_file)[0] + "_NtermSummary.csv" # Name of file that will contain summary DataFrame.
	summary_df_columns = ["original_seqID", "original_seq_length", "nterminal_seqID", "nterm_seq_length", "nterm_type", "nterm_typesubtype", "nterm_sequence"]
	summary_list = list(zip(original_seqID, original_seqlength, nterm_seqID, nterm_seqlength, nterm_type, nterm_typesubtype, nterm_seq))
	summary_df = pd.DataFrame(summary_list, columns = summary_df_columns)
	summary_df.to_csv(summary_df_file, index = False)


# Arguments entered on command line.
parser = argparse.ArgumentParser(description="Identify the N-terminal domain type and subtype of one or more VSG sequences.")
parser.add_argument('file', help="FASTA file containing VSG protein sequences.", action="store")
parser.add_argument('path', help="Path to directory containing VSG N-terminal TypeA/TypeB/TypeSubtype HMM profiles.", action="store")

args = parser.parse_args()
infile = args.file
path_hmm = args.path


# Confirm that HMMER is installed in PATH.
try:
	subprocess.check_output(["hmmscan","-h"])
except subprocess.CalledProcessError:
	sys.exit("HMMER not installed in PATH.")


infile_base = os.path.splitext(infile)[0] # Input filename as string w/o extension.

sys.stdout = open(infile_base + "_Nterm_stdout.txt", 'w') # Only errors will be printed to the UNIX interpreter.

if not os.path.exists(infile):
	sys.exit("%r does not exist." % infile)
else:
	print "%r found.\n" % infile
	if not os.path.exists(path_hmm):
		sys.exit("Path to HMM profiles does not exist.")
	else:
		if not path_hmm.endswith("/"):
			path_hmm += "/"

seq_variants_file = infile_base + "_Nterm_variants.fa" # Name of file that will contain sequence variants.
if os.path.exists(seq_variants_file):
	sys.exit("%r already exists." % seq_variants_file)


all_original_seqIDs = []

for sequence in SeqIO.parse(infile, "fasta"):
	all_original_seqIDs.append(sequence.id) # Append the original sequence ID to a list.
	if sequence.seq.endswith("*"): # If the sequence has a stop codon identifier.
		untrimmed_sequence = SeqRecord(
		seq = sequence.seq.rstrip("*"), # Remove identifier.
		id = sequence.id + "_Ntrimmed0", # Modify the original sequence ID to indicate untrimmed status.
		description = "",
		)
	else:
		untrimmed_sequence = SeqRecord(
		seq = sequence.seq,
		id = sequence.id + "_Ntrimmed0", # Modify the original sequence ID to indicate untrimmed status.
		description = "",
		)
	seq_variants = []
	seq_variants.append(untrimmed_sequence)
	print "Trimming of %r in progress..." % sequence.id
	seq_variants_complete = iterate_trim(seq_variants[0], seq_variants) # List of trimmed variants of the sequence.
	with open(seq_variants_file, "a") as seq_variants_outfile:
		SeqIO.write(seq_variants_complete, seq_variants_outfile, "fasta") # Append trimmed variants of sequence to a file.

print "\nProceeding to HMMER hmmscan analysis that will determine the type of the N-terminal domain(s)..."
hmmscan_type(seq_variants_file, path_hmm)
print "HMMER hmmscan analysis complete.\n"

for file in os.listdir(os.getcwd()):
	if file.endswith("TypeA.out") or file.endswith("TypeB.out"):
		empty = "#" # If there are no hmmscan hits, then the fourth line in hmmscan .out file contains "#".
		fourth_line = str(linecache.getline(file, 4))
		if not empty in fourth_line:
			type_letter = os.path.splitext(str(file))[0][-1] # "A" or "B"
			nterm_seqIDs = nterm_type_ID(file, all_original_seqIDs) # List of sequence ID(s) of the N-terminal domain sequence(s).
			print "Found most probable Type%s N-terminal domain(s) in %r.\n" % (type_letter, file)
			nterms_file = infile_base + "_Nterm_Type" + type_letter + "_FINAL.fa" # Name of file that will contain the N-terminal domain sequence(s).
			for sequence in SeqIO.parse(seq_variants_file, "fasta"):
				for nterm_seqID in nterm_seqIDs:
					if nterm_seqID == sequence.id:
						nterm_seqRecord = SeqRecord(
						seq = sequence.seq,
						id = sequence.id + "_Type" + type_letter, # Add the VSG N-terminal domain Type to sequence ID.
						description = "",
						)
						with open(nterms_file, "a") as nterms_outfile: # Append N-terminal domain sequence to a file.
							SeqIO.write(nterm_seqRecord, nterms_outfile, "fasta")

# Concatenate the FASTA files containing the TypeA and/or TypeB N-terminal domain sequence(s).
seq_nterms_file_wType = infile_base + "_Nterm_wType.fa"
subprocess.call(['cat *_Nterm_Type*_FINAL.fa > '+seq_nterms_file_wType+''], shell=True)
subprocess.call(['rm *_Nterm_Type*_FINAL.fa'], shell=True)

if not os.path.exists(seq_nterms_file_wType):
	sys.exit("%r was not created." % seq_nterms_file_wType)
else:
	seq_nterms_file_woType = infile_base + "_Nterm_woType.fa"
	for sequence in SeqIO.parse(seq_nterms_file_wType, "fasta"):
		nterm_woType = SeqRecord(
		seq = sequence.seq,
		id = sequence.id.split("_Type")[0], # Modify N-terminal domain sequence ID to remove "_TypeX".
		description = "",
		)
		with open(seq_nterms_file_woType, "a") as seq_nterms_file_woType_outfile:
			SeqIO.write(nterm_woType, seq_nterms_file_woType_outfile, "fasta")

print "Proceeding to second HMMER hmmscan analysis that will determine the subtype of the identified N-terminal domain(s)..."
hmmscan_subtype(seq_nterms_file_woType, path_hmm)
print "HMMER hmmscan analysis complete.\n"

for file in os.listdir(os.getcwd()):
	if file.endswith("TypeABsubtype.out"):
		empty = "#" # If there are no hmmscan hits, then the fourth line in hmmscan .out file contains "#".
		fourth_line = str(linecache.getline(file, 4))
		if not empty in fourth_line:
			typed_nterms = nterm_subtype_ID(file, all_original_seqIDs) # Dictionary - {N-terminal domain sequence ID: type-subtype}.
			print "Found most probable type and subtype of the identified N-terminal domain(s) in %r.\n" % file
			nterms_typesubtype_file = infile_base + "_Nterm_wType_FINAL.fa" # Name of file that will contain the N-terminal domain sequence(s) with type-subtype.
			for sequence in SeqIO.parse(seq_variants_file, "fasta"):
				for nterm_seqID, nterm_seqTypeSubtype in typed_nterms.items():
					if nterm_seqID == sequence.id:
						typedsubtyped_nterm_seqRecord = SeqRecord(
						seq = sequence.seq,
						id = sequence.id + "_Type" + nterm_seqTypeSubtype, # Add the VSG N-terminal domain type-subtype to sequence ID.
						description = "",
						)
						with open(nterms_typesubtype_file, "a") as typedsubtyped_nterms_outfile: # Append N-terminal domain sequence (w/ type-subtype) to a file.
							SeqIO.write(typedsubtyped_nterm_seqRecord, typedsubtyped_nterms_outfile, "fasta")
			print "%r contains the most probable N-terminal domain(s) with the corresponding type-subtype." % nterms_typesubtype_file

print "Preparing summary dataframe..."
summary_csv(infile, nterms_typesubtype_file)
print "%r contains information about the identified N-terminal domain sequence(s).\n" % (infile_base+"_NtermSummary.csv")

# Move all HMMER hmmscan input and output files to a subdirectory.
dir_make = infile_base + "_HMMSCAN"
subprocess.call(['mkdir -p '+dir_make+''], shell=True)
subprocess.call(['mv '+seq_variants_file+' '+seq_nterms_file_wType+' *.out *hmmscan.txt '+dir_make+''], shell=True)
subprocess.call(['rm '+seq_nterms_file_woType+''], shell=True)

print "Finished with %r." % infile

