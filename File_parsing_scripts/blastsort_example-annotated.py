from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from sys import argv
from Bio.Blast import NCBIXML



result_handle = open(argv[1])

blast_records = NCBIXML.parse(result_handle)

hit_list = []
exclude_list = []



for blast_record in blast_records:
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			if hsp.expect < 1.0e-10: 
                                if not blast_record.query in hit_list:
                                        hit_list.append(str(blast_record.query))
					alignment_length = float(hsp.query_end - hsp.query_start + 1) #length of the highest-scoring alignment (actual length of the alignment, not of either VSG, but the region in the query that matches. this won't count gaps)
					percent_query_aligned = float((100 * alignment_length) / blast_record.query_letters) #percentage of the query which is covered by the highest-scoring alignment <- this is our y-axis value
					query_length = blast_record.query_letters #length of the query!
					hit_vsg_length = alignment.length #length of the matching VSG for the highest-scoring alignment (I checked a bunch of these in macvector and this length appears to be reported correctly)
					alignment_identity = float((100.0 * hsp.identities) / alignment_length) #percent identity within the highest scoring alignment
					overall_identity = float((100.0 * hsp.identities) /  query_length) #overall identity relative to the length of the query VSG (this is the one to filter for >99% identity)
					overall_identity2 = float((100.0 * hsp.identities) /  hit_vsg_length) #overall identity relative to the matching VSG for the highest-scoring alignment <- this can be low even for a good hit, if the assembled VSG is short. don't filter on this!
					print alignment.title
					print str(alignment.length)
					print str(alignment_length)
					print str(alignment_identity)
					print str(overall_identity)

