#!/bin/bash
set -e
DIR=/home/chenzh/My_project/mouse_smallseq_preimp
TMPD=$DIR/tmp_data
SRC=$DIR/src
DATA=$DIR/data
RE=$DIR/results
DOC=$DIR/doc
BDOC=/home/chenzh/My_project/Extra_smncRNA/big_doc
BIN=$DIR/bin

cd $DIR
## bash smallrna_work_driver.bash E5_24_33
# failed EC10
#E4_8_3
SA=$1

WD=$TMPD/extra_HB_smallseq/${SA}
cutada=$BDOC/sc_smallRNA_annotation/Human/cutadapt_3prime.fa ### From paper
RGB=$BDOC/sc_smallRNA_annotation/Human/RefSeq/clean_chr_refdata/bowtie/bowtie
RGH=$BDOC/sc_smallRNA_annotation/Human/RefSeq/clean_chr_refdata/hisat/genome.trans
REF=$BDOC/sc_smallRNA_annotation/Human/RefSeq/clean_chr_refdata/fasta/genome.fa

BREF=$BDOC/sc_smallRNA_annotation/Human/RefSeq/clean_chr_refdata/bowtie/bowtie
HsaM=$BDOC/sc_smallRNA_annotation/Human/mirBase/hsa.mature.fa
NHsaM=$BDOC/sc_smallRNA_annotation/Human/mirBase/non.hsa.mature.fa
HsaH=$BDOC/sc_smallRNA_annotation/Human/mirBase/hsa.hairpin.fa
Hsa_GTF=$BDOC/sc_smallRNA_annotation/Human/mirBase/hsa.gff3


mkdir -p $WD
cd $WD


N=69 ## 84-8-2-5 ### further require the reads  were overlapped with at least 5bp adapter

#######
#SA=E4_8_3
#SA=E3_1_5
#cd ~/My_project/Extra_smncRNA/tmp_data/raw_data/SecondBatch/E4_12_4/TEMP
#zcat /home/chenzh/My_project/Extra_smncRNA/data/P12709/P12709_1403/02-FASTQ/190305_A00187_0122_AH75G5DRXX/P12709_1403_S392_L002_R1_001.fastq.gz > /home/chenzh/My_project/Extra_smncRNA/tmp_data/raw_data/MergeBatch_101bp/E4_8_3/E4_8_3.fastq
#######

source /home/chenzh/miniconda3/bin/activate small
zcat ../${SA}.fastq.gz > ${SA}.fastq
umi_tools extract --bc-pattern=NNNNNNNN -I ${SA}.fastq -S ${SA}.extract.fastq -L ${SA}.extract.log
cutadapt -a file:${cutada}   -e 0.1 -o 3 -m 18 -M ${N} -u 2 -o ${SA}.extract.cutadp.fastq ${SA}.extract.fastq  >${SA}.extract.cutadp.log  ## trim 2bp from begining, overlap at least 10bp 
#rm ${SA}.fastq ${SA}.extract.fastq
rm  -f ${SA}.extract.fastq ${SA}.fastq


bowtie  -a --best --strata -v 2 -m 50  -S -q -p 2 $RGB ${SA}.extract.cutadp.fastq  --al aligned.T3.0bp.fastq --un unaligned.T3.0bp.reads.fastq 2>bowie.log |python $SRC/filter_sam_mismatch.py |samtools view -Sb -F 4 > aligned.T3.0bp.bam  
seqkit fx2tab aligned.T3.0bp.fastq|sed -e 's/ /_/' |cut -f 1 > aligned.T3.0bp.ID
rm -f aligned.T3.0bp.fastq
#|python $SRC/filter_sam_mismatch.py|samtools view -Sb -F 4 > ${SA}.bowtie.bam 

#' trim 3'bp reads
for ((i=1; $i<4; i=$i+1))
do
   previous=$(($i-1));
   bowtie  -3 $i -a --best --strata -v 1 -m 50  -S -q -p 2 $RGB unaligned.T3.${previous}bp.reads.fastq --al aligned.T3.${i}bp.fastq  --un unaligned.T3.${i}bp.reads.fastq  2>>bowie.log |python $SRC/filter_sam_mismatch.py |samtools view -Sb -F 4 > aligned.T3.${i}bp.bam  
   if [ -e aligned.T3.${i}bp.fastq ];then
   	seqkit fx2tab aligned.T3.${i}bp.fastq|sed -e 's/ /_/' |cut -f 1 > aligned.T3.${i}bp.ID
   fi
   touch aligned.T3.${i}bp.ID
   rm -f unaligned.T3.${previous}bp.reads.fastq aligned.T3.${i}bp.fastq
done

#' trim 5' bp
mv unaligned.T3.3bp.reads.fastq unaligned.T5.0bp.reads.fastq
for ((i=1; $i<4; i=$i+1))
do
   previous=$(($i-1));
   bowtie  -5 $i -a --best --strata -v 1 -m 50  -S -q -p 2 $RGB unaligned.T5.${previous}bp.reads.fastq --al aligned.T5.${i}bp.fastq  --un unaligned.T5.${i}bp.reads.fastq  2>>bowie.log |python $SRC/filter_sam_mismatch.py |samtools view -Sb -F 4 > aligned.T5.${i}bp.bam  
   if [ -e aligned.T5.${i}bp.fastq ];then
   	seqkit fx2tab aligned.T5.${i}bp.fastq|sed -e 's/ /_/' |cut -f 1 > aligned.T5.${i}bp.ID
   fi
   touch aligned.T5.${i}bp.ID
   rm -f unaligned.T5.${previous}bp.reads.fastq aligned.T5.${i}bp.fastq
done


samtools merge aligned.temp.bam aligned.T3.0bp.bam aligned.T3.1bp.bam aligned.T3.2bp.bam aligned.T3.3bp.bam aligned.T5.1bp.bam aligned.T5.2bp.bam aligned.T5.3bp.bam
rm -f aligned.T3.0bp.bam aligned.T3.1bp.bam aligned.T3.2bp.bam aligned.T3.3bp.bam aligned.T5.1bp.bam aligned.T5.2bp.bam aligned.T5.3bp.bam

samtools view -h aligned.temp.bam |awk -F "\t" 'length($10) > 17  || $1 ~ /^@/' | samtools view -bS - > aligned.bam
samtools sort aligned.bam  -o aligned.sort.bam
samtools index aligned.sort.bam
rm -f unaligned.T5.3bp.reads.fastq aligned.bam aligned.temp.bam
#' UMI deduplication
umi_tools dedup --method directional --output-stats dedup.log1  --read-length -I aligned.sort.bam -S aligned.sort.dedup.temp.bam --random-seed 123 > /dev/null

samtools sort aligned.sort.dedup.temp.bam -o aligned.sort.dedup.bam
picard SamToFastq INPUT=aligned.sort.dedup.bam FASTQ=aligned.sort.dedup.bam.fastq.temp  2>> log.txt
seqkit fx2tab aligned.sort.dedup.bam.fastq.temp|sort |uniq |seqkit tab2fx -w 0 > aligned.sort.dedup.bam.fastq #### need save
bowtie  -a --best --strata -v 2 -m 50  -S -q -p 2 $RGB aligned.sort.dedup.bam.fastq 2>bowie.dedup.log |python $SRC/filter_sam_mismatch.py|samtools view -Sb -F 4 > aligned.sort.dedup.reMap.bam  ### recover the multiple mapping results
samtools sort aligned.sort.dedup.reMap.bam -o aligned.sort.dedup.reMap.sort.bam
rm -f aligned.sort.bam aligned.sort.dedup.temp.bam aligned.sort.dedup.bam.fastq.temp dedup.log1_edit_distance.tsv dedup.log1_per_umi_per_position.tsv dedup.log1_per_umi.tsv aligned.sort.bam.bai aligned.sort.dedup.bam

samtools view aligned.sort.dedup.reMap.bam |cut -f 1|sort |uniq > aligned.sort.dedup.reMap.bam.fastq.ID
seqkit grep -f <( cat aligned.sort.dedup.reMap.bam.fastq.ID) ${SA}.extract.cutadp.fastq > aligned.sort.dedup.reMap.bam.original.fastq


rm -f ${SA}.extract.cutadp.fastq aligned.sort.dedup.reMap.bam.fastq.ID

cut -f 1 aligned.sort.dedup.reMap.bam.original.fastq  -d " "|seqkit fx2tab  |awk -F "\t" '{ OFS="\t";print $1,length($2),substr($2,0,3)substr($2,length($2)-2,3)}' |sort -k2,2n -t $'\t'>  aligned.sort.dedup.reMap.bam.original.fastq.T5T3.3bp.infor ### IMP


bam2bed <aligned.sort.dedup.reMap.bam |cut -f 1-8 > temp.aligned.sort.dedup.reMap.bam.bed


echo -n "" > aligned.fastq.ID
for i in T3.0bp T3.1bp T3.2bp T3.3bp T5.1bp T5.2bp T5.3bp;do sed -e 's/_1:N:0:/\t/' aligned.${i}.ID| awk -F "\t" -v a=${i} '{print $1"\t"a}' >> aligned.fastq.ID;done
rm -rf aligned.T3.[0-3]bp.ID aligned.T5.[1-3]bp.ID 
join -1 1 -2 4 -t$'\t' <(sort -k1,1 aligned.fastq.ID ) <(sort -k4,4 temp.aligned.sort.dedup.reMap.bam.bed)|python $SRC/create_T5_T3_ref_bed.py  
#temp.aligned.sort.dedup.reMap.bam.T5_3bp.bed temp.aligned.sort.dedup.reMap.bam.T3_3bp.bed


join -1 2 -2 1 -t$'\t' <(join -1 1 -2 1 -t$'\t' <(awk -F "\t" '{if ($2 >0) print $0}' temp.aligned.sort.dedup.reMap.bam.T5_3bp.bed | fastaFromBed -fi ${REF} -s -name -tab  -bed - |sort -k1,1 ) <(awk -F "\t" '{if ($2 >0) print $0}' temp.aligned.sort.dedup.reMap.bam.T3_3bp.bed | fastaFromBed -fi ${REF} -s -name -tab  -bed - |sort -k1,1 ) |sed -e 's/(+)//'|sed -e 's/(-)//'|awk -F "\t" '{print $1"\t"$1"\t"$2$3}'|sed -e  's/_T[3,5].[0-3]bp_/\t/'|sort -k2,2)  <(sort -k1,1 aligned.sort.dedup.reMap.bam.original.fastq.T5T3.3bp.infor) |cut -f 2-6 |sed -e 's/_/\t/'|sed -e 's/_/\t/'|sed -e 's/_/\t/'|awk -F "\t" '{OFS="\t";print $1,$2,$3,$5,"len",$4,$6,$7,$8}' > temp

echo -n > temp.aligned.sort.dedup.reMap.bam.ov.small.detail ### IMP
for i in miRNA rRNA snoRNA snRNA tRNA piRNA os-piRNA #### mapping order
do
	intersectBed -a temp -b <(awk -v a=${i} -F "\t" '{if ($8==a) print $0}' $BDOC/sc_smallRNA_annotation/Human/small.anno.bed) -s -wa -wb >> temp.aligned.sort.dedup.reMap.bam.ov.small.detail
	intersectBed -a temp -b <(awk -v a=${i} -F "\t" '{if ($8==a) print $0}' $BDOC/sc_smallRNA_annotation/Human/small.anno.bed) -s -v > temp.remain
	cat temp.remain > temp
done



########## feature counting
#### Get counts
python $SRC/count_smallRNA.detailed.ordered.exp.withLen.py -l $BDOC/sc_smallRNA_annotation/Human/reads.report.order -s $BDOC/sc_smallRNA_annotation/Human/reads.report.order  -e temp.aligned.sort.dedup.reMap.bam.ov.small.detail -o aligned.detailed



