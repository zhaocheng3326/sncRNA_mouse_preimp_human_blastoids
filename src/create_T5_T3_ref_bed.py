import sys

output_T5=open("temp.aligned.sort.dedup.reMap.bam.T5_3bp.bed","w")
output_T3=open("temp.aligned.sort.dedup.reMap.bam.T3_3bp.bed","w")
for line in sys.stdin:
    line=line.strip("\n").split("\t")
    r_id = line[0]
    r_type= line[1]
    r_chr=line[2]
    r_begin = line[3]
    r_end = line[4]
    r_strand = line[6]
    if r_strand=="+":
        if r_type=="T3.0bp":
            o_T5_begin=r_begin; o_T5_end=str(int(r_begin)+3);
            o_T3_begin = str(int(r_end)-3); o_T3_end = r_end;
            o_id="_".join([r_chr,o_T5_begin,o_T3_end,r_strand,r_type,r_id])
        if r_type=="T3.1bp":
            o_T5_begin=r_begin; o_T5_end=str(int(r_begin)+3);
            o_T3_begin = str(int(r_end)-2); o_T3_end = str(int(r_end)+1);
            o_id="_".join([r_chr,o_T5_begin,o_T3_end,r_strand,r_type,r_id])
        if r_type=="T3.2bp":
            o_T5_begin=r_begin; o_T5_end=str(int(r_begin)+3);
            o_T3_begin = str(int(r_end)-1); o_T3_end = str(int(r_end)+2);
            o_id="_".join([r_chr,o_T5_begin,o_T3_end,r_strand,r_type,r_id])
        if r_type=="T3.3bp":
            o_T5_begin=r_begin; o_T5_end=str(int(r_begin)+3);
            o_T3_begin = str(int(r_end)); o_T3_end = str(int(r_end)+3);
            o_id="_".join([r_chr,o_T5_begin,o_T3_end,r_strand,r_type,r_id])
        if r_type=="T5.1bp":
            o_T5_begin=str(int(r_begin)-1); o_T5_end=str(int(r_begin)+2);
            o_T3_begin = str(int(r_end)-3); o_T3_end = r_end;
            o_id="_".join([r_chr,o_T5_begin,o_T3_end,r_strand,r_type,r_id])
        if r_type == "T5.2bp":
            o_T5_begin = str(int(r_begin) - 2);o_T5_end = str(int(r_begin) + 1);
            o_T3_begin = str(int(r_end) - 3); o_T3_end = r_end;
            o_id = "_".join([r_chr, o_T5_begin, o_T3_end, r_strand, r_type, r_id])
        if r_type == "T5.3bp":
            o_T5_begin = str(int(r_begin) - 3);o_T5_end = str(int(r_begin ));
            o_T3_begin = str(int(r_end) - 3); o_T3_end = r_end;
            o_id = "_".join([r_chr, o_T5_begin, o_T3_end, r_strand, r_type, r_id])
        print(r_chr,o_T5_begin,o_T5_end,o_id,r_type,r_strand,sep="\t",file=output_T5)
        print(r_chr, o_T3_begin, o_T3_end, o_id, r_type, r_strand, sep="\t", file=output_T3)

    elif r_strand=="-":
        if r_type=="T3.0bp":
            o_T5_begin=str(int(r_begin)); o_T5_end=str(int(r_begin)+3);
            o_T3_begin = str(int(r_end)-3); o_T3_end = str(int(r_end));
            o_id="_".join([r_chr,o_T5_begin,o_T3_end,r_strand,r_type,r_id])
        if r_type=="T5.1bp":
            o_T5_begin=r_begin; o_T5_end=str(int(r_begin)+3);
            o_T3_begin = str(int(r_end)-2); o_T3_end = str(int(r_end)+1);
            o_id="_".join([r_chr,o_T5_begin,o_T3_end,r_strand,r_type,r_id])
        if r_type=="T5.2bp":
            o_T5_begin=r_begin; o_T5_end=str(int(r_begin)+3);
            o_T3_begin = str(int(r_end)-1); o_T3_end = str(int(r_end)+2);
            o_id="_".join([r_chr,o_T5_begin,o_T3_end,r_strand,r_type,r_id])
        if r_type=="T5.3bp":
            o_T5_begin=r_begin; o_T5_end=str(int(r_begin)+3);
            o_T3_begin = str(int(r_end)); o_T3_end = str(int(r_end)+3);
            o_id="_".join([r_chr,o_T5_begin,o_T3_end,r_strand,r_type,r_id])
        if r_type=="T3.1bp":
            o_T5_begin=str(int(r_begin)-1); o_T5_end=str(int(r_begin)+2);
            o_T3_begin = str(int(r_end)-3); o_T3_end = r_end;
            o_id="_".join([r_chr,o_T5_begin,o_T3_end,r_strand,r_type,r_id])
        if r_type == "T3.2bp":
            o_T5_begin = str(int(r_begin) - 2);o_T5_end = str(int(r_begin) + 1);
            o_T3_begin = str(int(r_end) - 3); o_T3_end = r_end;
            o_id = "_".join([r_chr, o_T5_begin, o_T3_end, r_strand, r_type, r_id])
        if r_type == "T3.3bp":
            o_T5_begin = str(int(r_begin) - 3);o_T5_end = str(int(r_begin ));
            o_T3_begin = str(int(r_end) - 3); o_T3_end = r_end;
            o_id = "_".join([r_chr, o_T5_begin, o_T3_end, r_strand, r_type, r_id])
        print(r_chr,o_T5_begin,o_T5_end,o_id,r_type,r_strand,sep="\t",file=output_T3)
        print(r_chr, o_T3_begin, o_T3_end, o_id, r_type, r_strand, sep="\t", file=output_T5)

output_T5.close()
output_T3.close()

