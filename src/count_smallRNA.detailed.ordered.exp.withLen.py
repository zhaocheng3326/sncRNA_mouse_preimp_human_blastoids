import sys 
#sys.path.insert(0, './')
#import dr_tools, pysam
import argparse, os


parser = argparse.ArgumentParser()
parser.add_argument('-l', '--long_anno', help='long reads order',required=True)
parser.add_argument('-s', '--short_anno', help='short reads order',required=True)
#parser.add_argument('-r', '--reads_length', help='read length distribution',required=True)
parser.add_argument('-e', '--InforDetailed', required=True, help='inforDetailed: 17 columns (small RNA reads with ref bp , read length, annotated pos) ')
parser.add_argument('-o', '--output', required=True, help='output directory')
o = parser.parse_args()


type_order_long={}
with open(o.long_anno,mode="r") as file:
    for line in file:
        line=line.strip("\n").split("\t")
        type_order_long[line[0]]=int(line[1])


type_order_short={}
with open(o.short_anno,mode="r") as file:
    for line in file:
        line=line.strip("\n").split("\t")
        type_order_short[line[0]]=int(line[1])


reads_brief={}
reads_type={}
reads_length={}
reads_pos={}

with open(o.InforDetailed,mode="r") as file:
    for line in file:
        line=line.strip("\n").split("\t")
        if line[0]=="mapped_pos_readID":
            continue
        r_chr=line[0]
        r_begin=line[1]  
        r_end=line[2]
        r_strand=line[3].split("_")[3]
        if not r_strand in ["+","-"]:
            continue
        r_id_shift_type=line[3].split("_")[4]
        r_id=line[3].split("_")[5]+"_"+line[3].split("_")[6]
        r_refbp=line[6]
        r_length=line[7]
        r_fabp=line[8]
        r_combp=""
        r_pos=":".join([r_chr,r_begin,r_end,r_strand])
        for n in range(0,len(r_refbp)):
            if r_refbp[n]==r_fabp[n]:
                r_combp=r_combp+"*"
            else:
                r_combp=r_combp+r_fabp[n]


        m_begin=line[10]
        m_end=line[11]
        m_id=line[15]
        m_type=line[16]
        if r_strand=="+":
            if int(m_begin) +10 < int(r_begin):
                rm_5p="Lv5plt"
            elif int(m_begin)  < int(r_begin):
                rm_5p="Lv5pt"
            elif int(m_begin) == int(r_begin):
                rm_5p = "Lv5p"
            elif int(m_begin)-10 > int(r_begin):
                rm_5p = "Lv5ple"
            else:
                rm_5p = "Lv5pe"
            if int(r_end) +10 < int(m_end):
                rm_3p="Lv3plt"
            elif int(r_end) < int(m_end):
                rm_3p="Lv3pt"
            elif int(r_end) == int(m_end):
                rm_3p="Lv3p"
            elif int(r_end) > int(m_end)+10:
                rm_3p="Lv3ple"
            elif int(r_end) > int(m_end):
                rm_3p = "Lv3pe"
        elif r_strand=="-":
            if int(m_begin) +10 < int(r_begin):
                rm_3p="Lv3plt"
            elif int(m_begin)  < int(r_begin):
                rm_3p="Lv3pt"
            elif int(m_begin) == int(r_begin):
                rm_3p = "Lv3p"
            elif int(m_begin)-10 > int(r_begin):
                rm_3p = "Lv3ple"
            else:
                rm_3p = "Lv3pe"
            if int(r_end) +10 < int(m_end):
                rm_5p="Lv5plt"
            elif int(r_end) < int(m_end):
                rm_5p="Lv5pt"
            elif int(r_end) == int(m_end):
                rm_5p="Lv5p"
            elif int(r_end) > int(m_end)+10:
                rm_5p="Lv5ple"
            elif int(r_end) > int(m_end):
                rm_5p = "Lv5pe"
        r_ft=":".join([m_id,m_type,r_length,r_refbp,r_fabp,r_combp,rm_5p,rm_3p])
        if int(r_length) < 50:
            if m_type not in type_order_short:
                m_type="OTHER"
            if r_id in reads_brief:
                if type_order_short[m_type] > type_order_short[reads_type[r_id]]:
                    continue
                elif type_order_short[m_type] == type_order_short[reads_type[r_id]]:
                    if r_ft not in reads_brief[r_id]:
                        reads_brief[r_id].append(r_ft)
                        reads_pos[r_ft] = r_pos
                elif type_order_short[m_type] < type_order_short[reads_type[r_id]]:
                    reads_brief[r_id]=[]
                    reads_brief[r_id].append(r_ft)
                    reads_pos[r_ft]=r_pos
                    reads_type[r_id]=m_type
            else:
                reads_brief[r_id]=[]
                reads_brief[r_id].append(r_ft)
                reads_pos[r_ft] = r_pos
                reads_type[r_id]=m_type
        else:
            if m_type not in type_order_long:
                m_type="OTHER"
            if r_id in reads_brief:
                if type_order_long[m_type] > type_order_long[reads_type[r_id]]:
                    continue
                elif type_order_long[m_type] == type_order_long[reads_type[r_id]]:
                    if r_ft not in reads_brief[r_id]:
                        reads_brief[r_id].append(r_ft)
                        reads_pos[r_ft] = r_pos
                elif type_order_long[m_type] < type_order_long[reads_type[r_id]]:
                    reads_brief[r_id]=[]
                    reads_brief[r_id].append(r_ft)
                    reads_type[r_id]=m_type
                    reads_pos[r_ft] = r_pos
            else:
                reads_brief[r_id]=[]
                reads_brief[r_id].append(r_ft)
                reads_type[r_id]=m_type
                reads_pos[r_ft] = r_pos


with open(o.output+".od.exp.tmp.txt","w") as output_exp_tmp:
    anno_count={}
    for read in reads_brief:
        for anno in reads_brief[read]:
            if anno not in anno_count:
                anno_count[anno]=1/len(reads_brief[read])
            else:
                anno_count[anno]=1/len(reads_brief[read])+ anno_count[anno]
            anno_pos=reads_pos[anno].split(":")
            anno_id=anno.split(":")[0]
            print(anno_pos[0],anno_pos[1],anno_pos[2],read,"len",anno_pos[3],anno_id,anno,sep="\t",file=output_exp_tmp)
output_exp_tmp.close()

with open(o.output + ".od.exp.txt", "w") as output_exp:
    for anno in anno_count:
        print(anno, round(anno_count[anno], 3), sep="\t", file=output_exp)
output_exp.close()
