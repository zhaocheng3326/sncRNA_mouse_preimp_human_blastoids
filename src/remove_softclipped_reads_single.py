from __future__ import division
import os, pysam, argparse

"""
This script:
	- removes all hardclipped reads, reads with insertions and deletions
	- removes reads with any 5' soft-clipping
	- removes reads with more than 3 bases on 3' via soft-clipping 
 
"""

def remove_clipped_reads(inbam,wofbam,outbam):
	inbamPysamObj = pysam.Samfile(inbam, "rb" )
	wofbamPysamObj = pysam.Samfile(wofbam, "wb", template=inbamPysamObj)
	outbamPysamObj = pysam.Samfile(outbam, "wb", template=inbamPysamObj)
	total_reads = 0; reads_after_clip_removal = 0
	for read in inbamPysamObj:
		total_reads += 1
		read_name = read.qname
		tid = read.rname
		readchr  = inbamPysamObj.getrname(tid)
		cigar = read.cigar  #this is a list of tuples;  example: cigar=[(0, 25), (4, 8)]  

		num_softclipped_5p = 0
		num_softclipped_3p = 0
		num_hardclipped = 0
		num_insertions = 0
		num_deletions = 0

		for i,e in enumerate(cigar):
			if e[0] == 5:
				num_hardclipped += e[1]
			elif e[0] == 1:
				num_insertions += e[1]
			elif e[0] == 2:
				num_deletions += e[1]
			elif i == 0 and e[0] == 4:
				num_softclipped_5p += e[1]
			elif i == len(cigar)-1 and e[0] == 4:
				num_softclipped_3p += e[1]

		if num_softclipped_5p > o.allowed_num_clipped_bases_5p: continue
		elif num_softclipped_3p > o.allowed_num_clipped_bases_3p: continue
		elif num_hardclipped > 0: continue
		elif num_insertions > 0: continue 
		elif num_deletions > 0: continue   

		if len(read.query_sequence) <= o.maximum_smRNA:
			outbamPysamObj.write(read)
		else:
			wofbamPysamObj.write(read)
		reads_after_clip_removal += 1
	wofbamPysamObj.close()
	outbamPysamObj.close()
	pysam.index(outbam, template=inbamPysamObj)
	pysam.index(wofbam, template=inbamPysamObj)

	print("%s percentage of removed reads = %.2f" %(os.path.basename(inbam), 100*(1 -reads_after_clip_removal/total_reads)))


#main function
if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', required=True)
	parser.add_argument('-w', '--withoutfilter', required=True)
	parser.add_argument('-o', '--output', required=True)
	parser.add_argument('-m', '--maximum_smRNA', type=int, default=40, help='maximum read length to determine small RNA')
	parser.add_argument('-f', '--allowed_num_clipped_bases_5p', type=int, default=0, help='both hard and soft clipped, S:4, H:5')
	parser.add_argument('-t', '--allowed_num_clipped_bases_3p', type=int, default=3, help='both hard and soft clipped, S:4, H:5')
	o = parser.parse_args()
	inbam = o.input
	wofbam = o.withoutfilter
	outbam = o.output
	remove_clipped_reads(inbam,wofbam,outbam)


