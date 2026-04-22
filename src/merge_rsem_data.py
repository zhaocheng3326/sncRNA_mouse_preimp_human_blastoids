#!/usr/bin/env python

"""
"""
# standard python libraries
import os, sys, re
import pandas as pd
from argparse import ArgumentParser
parser = ArgumentParser(description='Merge rpkmforgenes info from several samples')

# [Required input]
parser.add_argument('-i', '--infiles', metavar='infiles', default=[], help="Input rpkm-files, separated by space", required=False,nargs='+')
parser.add_argument('-l', '--samplelist', metavar='samplelist', default=[], help="File with a list of all rpkm-files to merge", required=False)
parser.add_argument('-o', '--outprefix', metavar='outprefix', help='Output prefix', required=True)
parser.add_argument('-f', '--file_extension', metavar='file_extension', help='Remove these strings from the file name to get a sample name', required=False,nargs='+')
parser.add_argument('-s', '--separator', metavar='separator', default=",", help='Separator to use in output files', required=False)

args = parser.parse_args()

outfile_rpkm = args.outprefix + ".rsem_rpkm.csv"
outfile_tpm = args.outprefix + ".rsem_tpm.csv"
outfile_counts = args.outprefix + ".rsem_counts.csv"

if bool(args.samplelist) == bool(args.infiles):
    raise ValueError("Exactly one of --infiles and --samplefile must be set")

if args.samplelist:
    files = []
    for line in open(args.samplelist):
        line = line.strip()
        files.append(line)
else:
    files = args.infiles


# read in data
counts = pd.DataFrame()
tpm = pd.DataFrame()
rpkm = pd.DataFrame()
for file in files:
    pdf = pd.read_csv(file,header=0,sep="\t")
    sname = os.path.split(file)[1]
    for replace in args.file_extension:
        sname = re.sub(replace,'',sname)
    rpkm[sname] = pdf["FPKM"]
    counts[sname] = pdf["expected_count"]
    tpm[sname] = pdf["TPM"]


# make one file with information about the genes, transcripts,
#gene_info = pdf
#gene_info = gene_info.loc[,["gene_id","transcript_id(s)"]]
#gene_info.to_csv(outfile_genes,sep=",")


# fix coulmn and row names and write to file.
rpkm.index = pdf['gene_id']
# add up data from duplicated gene names
rpkm=rpkm.groupby(rpkm.index).sum()
# write to file
rpkm.to_csv(outfile_rpkm,sep=args.separator)

tpm.index = pdf['gene_id']
tpm=tpm.groupby(tpm.index).sum()
tpm.to_csv(outfile_tpm,sep=args.separator)


counts.index = pdf['gene_id']
counts = counts.astype(int)
counts=counts.groupby(counts.index).sum()
counts.to_csv(outfile_counts,sep=args.separator)
