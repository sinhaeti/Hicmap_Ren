#!/bin/bash

# test if softwares installed
command -v bwa >/dev/null 2>&1 || { echo >&2 "I require bwa but it's not installed.  Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "I require samtools but it's not installed.  Aborting."; exit 1; }

# pass paramters
usage(){
cat << EOF
usage: ${0##*/} [-h] [-f FASTQ1] [-r FASTQ2] [-n PREFIX] [-g BWA_GENOME] [-c CUT_SITES] [-m MIN_INSERT_SIZE]

Map and processing Hi-C reads
(1) map FASTQ1/FASTQ2 using BWA indepedently;
(2) filter reads with MAPQ < 10;
(3) filter reads fell >500bp far from restriction enzyme cutter sites (strand sensitive);
(4) sort reads by read names;
(5) pair up two ends and filter invalid hic pairs;
(6) Remove PCR duplication using Picard - markDuplicates;

Example:
	bash bin/hicmap.sh -t 20 -m 8G -f data/JL_H4_R1.fastq.bz2 -r data/JL_H4_R2.fastq.bz2 -p scripts/MarkDuplicates.jar -n JL_H4 -g /oasis/tscc/scratch/r3fang/data/Mus_musculus/UCSC/mm9/Sequence/BWAIndex/genome.fa -c data/mm9.MboI.500bp -d 1000

Options:    
	-h, --help			show this help message and exit.
	-t  THREADS			threads [1].
	-m  MAX_MEM			max memory usage [4G].
	-f  FASTQ1			first mate of pair-end sequencing data [.fq/.fastq/.gz/.bz2].
	-r  FASTQ2			second mate of pair-end sequencing data [.fq/.fastq/.gz/.bz2].
	-p  MARK_DUPLICATE  		path to picard MarkDuplicates.jar
	-n  NAME			prefix of output files.
	-g  BWA_GENOME			BWA indexed reference genome.
	-c  CUT_ENZ			restriction cutting enzyme files. 
	-d  MIN_INSERT_SIZE		min insert size for valid "DIFFERENT-STRAND" pairs.
EOF
} 

THREADS=1
MAX_MEM="4G"

while getopts ":t:m:f:r:p:n:g:c:d:" opt;
do
	case "$opt" in
		t) THREADS=$OPTARG;;
		m) MAX_MEM=$OPTARG;;
		f) FASTQ1=$OPTARG;;
		r) FASTQ2=$OPTARG;;
		p) MARK_DUPLICATE=$OPTARG;;
		n) PREFIX=$OPTARG;;
		g) GENOME=$OPTARG;;
		c) CUT_ENZ=$OPTARG;;
		d) MIN_INSERT_SIZE=$OPTARG;;
		\?) usage
			exit 1
			;;
	esac
done

if [ $# -lt 10 ] ; then
   usage
   echo "error: too few arguments"
   exit 1
fi

re='^[0-9]+$'
if ! [[ $THREADS =~ $re ]] ; then
   echo "error: '$THREADS' Not a number" >&2; 
   exit 1
fi

# check if input file exists
if [ ! -f $FASTQ1 ]; then
	usage
    echo "error: '$FASTQ1' not exists.";
	exit 1
fi

if [ ! -f $FASTQ2 ]; then
	usage
    echo "error: '$FASTQ2' not exists.";
	exit 1
fi

if [ ! -f $GENOME ]; then
	usage
    echo "error: '$GENOME' not exists.";
	exit 1
fi

if [ ! -f $CUT_ENZ.neg.merged.bed ]; then
	usage
    echo "error: '$CUT_ENZ.neg.merged.bed' not exists.";
	exit 1
fi

if [ ! -f $CUT_ENZ.pos.merged.bed ]; then
	usage
    echo "error: '$CUT_ENZ.pos.merged.bed' not exists.";
	exit 1
fi

# check if input type
re='^[0-9]+$'
if ! [[ $MIN_INSERT_SIZE =~ $re ]] ; then
   echo "error: '$MIN_INSERT_SIZE' Not a number" >&2; 
   exit 1
fi

mkdir $PREFIX\_tmp

echo "Step1. map reads then filter non-uniquely and secondary alignment" 
if [ ${FASTQ1: -4} == ".bz2" ];
then
	bwa mem -t $THREADS $GENOME <(bzip2 -dc $FASTQ1) | samtools view -F 2048 -q 10 -bS - > $PREFIX\_tmp/$PREFIX\_R1.uniq.bam
else
	bwa mem -t $THREADS $GENOME $FASTQ1 | samtools view -F 2048 -q 10 -bS - > $PREFIX\_tmp/$PREFIX\_R1.uniq.bam
fi

if [ ${FASTQ2: -4} == ".bz2" ];
then
	bwa mem -t $THREADS $GENOME <(bzip2 -dc $FASTQ2) | samtools view -F 2048 -q 10 -bS - > $PREFIX\_tmp/$PREFIX\_R2.uniq.bam
else
	bwa mem -t $THREADS $GENOME $FASTQ2 | samtools view -F 2048 -q 10 -bS - > $PREFIX\_tmp/$PREFIX\_R2.uniq.bam
fi

samtools flagstat $PREFIX\_tmp/$PREFIX\_R1.uniq.bam > $PREFIX\_tmp/$PREFIX\_R1.uniq.bam.flagstat &
samtools flagstat $PREFIX\_tmp/$PREFIX\_R2.uniq.bam > $PREFIX\_tmp/$PREFIX\_R2.uniq.bam.flagstat &

echo "Step2. filter reads that are  > 500bp far from restriction cutter sites" 
# positive strand
samtools view -b -F 16 -L $CUT_ENZ.pos.merged.bed $PREFIX\_tmp/$PREFIX\_R1.uniq.bam > $PREFIX\_tmp/$PREFIX\_R1.uniq.filtered.pos.bam 
samtools view -b -F 16 -L $CUT_ENZ.pos.merged.bed $PREFIX\_tmp/$PREFIX\_R2.uniq.bam > $PREFIX\_tmp/$PREFIX\_R2.uniq.filtered.pos.bam 

# negative strand
samtools view -b -f 16 -L $CUT_ENZ.neg.merged.bed $PREFIX\_tmp/$PREFIX\_R1.uniq.bam > $PREFIX\_tmp/$PREFIX\_R1.uniq.filtered.neg.bam 
samtools view -b -f 16 -L $CUT_ENZ.neg.merged.bed $PREFIX\_tmp/$PREFIX\_R2.uniq.bam > $PREFIX\_tmp/$PREFIX\_R2.uniq.filtered.neg.bam 

# merge both strands
samtools cat -o $PREFIX\_tmp/$PREFIX\_R1.uniq.filtered.bam $PREFIX\_tmp/$PREFIX\_R1.uniq.filtered.pos.bam $PREFIX\_tmp/$PREFIX\_R1.uniq.filtered.neg.bam 
samtools cat -o $PREFIX\_tmp/$PREFIX\_R2.uniq.filtered.bam $PREFIX\_tmp/$PREFIX\_R2.uniq.filtered.pos.bam $PREFIX\_tmp/$PREFIX\_R2.uniq.filtered.neg.bam 

echo "Step3. sort reads by read names" 
samtools sort -n $PREFIX\_R1.uniq.filtered.bam $PREFIX\_tmp/$PREFIX\_R1.uniq.filtered.sorted 
samtools sort -n $PREFIX\_R2.uniq.filtered.bam $PREFIX\_tmp/$PREFIX\_R2.uniq.filtered.sorted 

echo "Step4. pair up two ends" 
pair2mates -m $MIN_INSERT_SIZE -o $PREFIX\_tmp/$PREFIX.filtered.paired.bam $PREFIX\_tmp/$PREFIX\_R1.filtered.sorted.bam $PREFIX\_tmp/$PREFIX\_R2.filtered.sorted.bam

echo "Step5. sort based on genomic coordinates" 
samtools sort -m $MAX_MEM $PREFIX\_tmp/$PREFIX.filtered.paired.bam $PREFIX\_tmp/$PREFIX.filtered.paired.sorted

echo "Step6. Filter PCR duplication"
java -Xmx10g -jar $MARK_DUPLICATE INPUT=$PREFIX\_tmp/$PREFIX.filtered.paired.sorted.bam OUTPUT=$PREFIX.filtered.paired.sorted.nodup.bam ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE=metrics.$PREFIX.txt TMP_DIR=$PREFIX\_tmp

#$SAMTOOLS flagstat $PREFIX.merged.sorted.paired.sorted.nodup.bam > $PREFIX.merged.sorted.paired.sorted.nodup.bam.flagstat &
#
#echo "Total Input:" >> hicmap.log
#bzip2 -dc $FASTQ1 | wc -l >> hicmap.log
#echo "Uniquely mapped for R1:" >> hicmap.log
#$SAMTOOLS view $PREFIX\_R1.bam | wc -l >> hicmap.log
#echo "Uniquely mapped for R2:" >> hicmap.log
#$SAMTOOLS view $PREFIX\_R2.bam | wc -l >> hicmap.log
#echo "Closed to cutter sites for R1:" >> hicmap.log
#$SAMTOOLS view $PREFIX\_R1.filtered.bam | wc -l >> hicmap.log
#echo "Closed to cutter sites for R2:" >> hicmap.log
#$SAMTOOLS view $PREFIX\_R2.filtered.bam | wc -l >> hicmap.log
#echo "Total usable read pairs:" >> hicmap.log
#$SAMTOOLS index $PREFIX.merged.sorted.paired.sorted.nodup.bam 
#$SAMTOOLS view $PREFIX.merged.sorted.paired.sorted.nodup.bam | wc -l >> hicmap.log 