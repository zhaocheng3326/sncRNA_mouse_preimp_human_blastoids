import sys, pysam
import os, argparse
from pyfasta import Fasta

"""
This script removes reads when 3' has TGG on the genome
Works for smallrna star folder. aka max32, max38

>adapter_three_prime
TGGAATTCTCGGGTGCCAAGG
>polyA
AAAAAAAAAAAAA
"""
def remove_reads_from_precursor(inbam,outbam,gr,minRlen,readlen_cutoff):
	"""
	prepare input/output files
	"""
	inbamPysamObj = pysam.Samfile(inbam, "rb" )
	outbamPysamObj = pysam.Samfile(outbam, "wb", template=inbamPysamObj)

	"""
	create genome fetch object
	"""
	gf = Fasta(gr)

	"""
	remove reads when 3' has TGG on the genome
	"""
	for read in inbamPysamObj:
		read_name = read.qname
		tid = read.rname
		readchr  = inbamPysamObj.getrname(tid)
		readstart = int(read.pos) + 1
		readend = read.aend
		strand = read.flag
		readlen = len(read.seq) #this is the actual read length (41M, means readlen=41)
		read_len = read.qlen  #this only considers matches (8S30M, means read_len=30)
		if readlen <= readlen_cutoff:
			outbamPysamObj.write(read)
			continue
		
#		if strand ==0 :  #read maps to forward strand
		if strand ==0 or strand ==256:  #read maps to forward strand
			upperlimit = minRlen - readlen
			#print(readchr,readend+1,readend+upperlimit)
			bpwindow = gf.sequence({'chr':readchr,'start':readend+1,'stop':readend+upperlimit})
                        #print bpwindow
			if readlen==minRlen-1 and (bpwindow == "T" or bpwindow == "A"): continue #TGGAATTCTCGGGTGCCAAGG
			elif readlen==minRlen-2 and (bpwindow == "TG" or bpwindow == "AA"): continue
			elif readlen==minRlen-3 and (bpwindow == "TGG" or bpwindow == "AAA"): continue
			elif readlen==minRlen-4 and (bpwindow == "TGGA" or bpwindow == "AAAA"): continue
			elif readlen==minRlen-5 and (bpwindow == "TGGAA" or bpwindow == "AAAAA"): continue
			else: outbamPysamObj.write(read)

#		elif strand ==16:  #read maps to reverse strand
		elif strand ==16 or strand ==272:  #read maps to reverse strand
			upperlimit = minRlen - readlen
			bpwindow = gf.sequence({'chr':readchr, 'start':readstart-upperlimit,'stop':readstart-1})
			if readlen==minRlen-1 and (bpwindow == "A" or bpwindow == "T"): continue #TTCCA
			elif readlen==minRlen-2 and (bpwindow == "CA" or bpwindow == "TT"): continue
			elif readlen==minRlen-3 and (bpwindow == "CCA" or bpwindow == "TTT"): continue
			elif readlen==minRlen-4 and (bpwindow == "TCCA" or bpwindow == "TTTT"): continue
			elif readlen==minRlen-5 and (bpwindow == "TTCCA" or bpwindow == "TTTTT"): continue
			else: outbamPysamObj.write(read)

	outbamPysamObj.close()
	#sort and index the final bam file
	#pysam.sort(outbam,outbam.rstrip(".bam")+".sort")
	#pysam.index(outbam.rstrip(".bam")+".sort.bam", template=inbamPysamObj)
	#os.remove(outbam)

#main function
if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--inputbam', required=True)
	parser.add_argument('-o', '--outputbam', required=True)
	parser.add_argument('-g', '--genome_ref', required=True)
	parser.add_argument('-c', '--readlen_cutoff', default=35)
	parser.add_argument('-x', '--minRlen', default=41, type=int)  #minimum read length to define a precursor
	o = parser.parse_args()
	
	remove_reads_from_precursor(o.inputbam,o.outputbam,o.genome_ref,o.minRlen,o.readlen_cutoff)

