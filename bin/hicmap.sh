#!/bin/bash

# test if softwares installed
command -v bwa >/dev/null 2>&1 || { echo >&2 "I require bwa but it's not installed.  Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "I require samtools but it's not installed.  Aborting."; exit 1; }

# pass paramters
usage(){
cat << EOF
usage: ${0##*/} [-h] [-f FASTQ1] [-r FASTQ2] [-n PREFIX] [-g BWA_GENOME] [-c CUT_SITES] [-m MIN_INSERT_SIZE]

Map and processing Hi-C reads
(1) map two ends using BWA;
(2) filter reads with MAPQ < 10;
(3) filter reads fell 500bp outside restriction enzyme cutter sites;
(4) sort reads by read names;
(5) pair up two ends: 
	i)   if pair is "SAME-STRAND" reads on "+" strand, keep reads only if they are 
	     within 500bp upstream of closest cutter sites;
   	ii)  if pair is "SAME-STRAND" reads on "_" strand, keep reads only if they are 
   	     within 500bp downstream of closest cutter sites;
	iii) if pair is "DIFFERENT-STRAND", keep it if insert size > 10000 (default) 
(6) Remove PCR duplication;

Example:
	bash bin/hicmap.sh -t 20 -f data/JL_H4_R1.fastq.bz2 -r data/JL_H4_R2.fastq.bz2 -n JL_H4 -g /mnt/thumper/home/r3fang/data/Mus_musculus/UCSC/mm9/Sequence/BWAIndex/genome.fa -c data/mm9.MboI.500bp -m 1000

Options:    
	-h, --help      	show this help message and exit.
	-t THREADS			threads [1].
	-f FASTQ1       	first mate of pair-end sequencing data.
	-r FASTQ1       	second mate of pair-end sequencing data.
	-n NAME            	prefix of output files.
	-g BWA_GENOME   	BWA indexed reference genome.
	-c CUT_ENZ      	restriction cutting enzyme files. 
	-m MIN_INSERT_SIZE	min insert size for valid "DIFFERENT-STRAND" pairs.
EOF
} 

while getopts ":t:f:r:n:g:c:m:" opt;
do
	case "$opt" in
		t) THREADS=$OPTARG;;
		f) FASTQ1=$OPTARG;;
		r) FASTQ2=$OPTARG;;
		n) PREFIX=$OPTARG;;
		g) GENOME=$OPTARG;;
		c) CUT_ENZ=$OPTARG;;
		m) MIN_INSERT_SIZE=$OPTARG;;
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


#echo "Step1. Mapping reads and filter non-uniquely and secondary alignment" 
#if [ ${FASTQ1: -4} == ".bz2" ];
#then
#	bwa mem -t $THREADS $GENOME <(bzip2 -dc $FASTQ1) | samtools view -F 2048 -q 10 -bS - >$PREFIX\_R1.uniq.bam
#else
#	bwa mem -t $THREADS $GENOME $FASTQ1 | samtools view -F 2048 -q 10 -bS - > $PREFIX\_R1.uniq.bam
#fi
#
#if [ ${FASTQ2: -4} == ".bz2" ];
#then
#	bwa mem -t $THREADS $GENOME <(bzip2 -dc $FASTQ2) | samtools view -F 2048 -q 10 -bS - > $PREFIX\_R2.uniq.bam
#else
#	bwa mem -t $THREADS $GENOME $FASTQ2 | samtools view -F 2048 -q 10 -bS - > $PREFIX\_R2.uniq.bam
#fi
#
#samtools flagstat $PREFIX\_R1.uniq.bam > $PREFIX\_R1.uniq.bam.flagstat &
#samtools flagstat $PREFIX\_R2.uniq.bam > $PREFIX\_R2.uniq.bam.flagstat &

echo "Step2. Filter reads that are far from restriction cutter sites" 
mkdir $PREFIX\_tmp
# positive strand
samtools view -b -F 16 -L $CUT_ENZ.pos.merged.bed $PREFIX\_R1.uniq.bam > $PREFIX\_tmp/$PREFIX\_R1.uniq.filtered.pos.bam 
samtools view -b -F 16 -L $CUT_ENZ.pos.merged.bed $PREFIX\_R2.uniq.bam > $PREFIX\_tmp/$PREFIX\_R2.uniq.filtered.pos.bam 

# negative strand
samtools view -b -f 16 -L $CUT_ENZ.neg.merged.bed $PREFIX\_R1.uniq.bam > $PREFIX\_tmp/$PREFIX\_R1.uniq.filtered.neg.bam 
samtools view -b -f 16 -L $CUT_ENZ.neg.merged.bed $PREFIX\_R2.uniq.bam > $PREFIX\_tmp/$PREFIX\_R2.uniq.filtered.neg.bam 

# merge both strands
samtools cat -o $PREFIX\_R1.uniq.filtered.bam $PREFIX\_tmp/$PREFIX\_R1.uniq.filtered.pos.bam $PREFIX\_tmp/$PREFIX\_R1.uniq.filtered.neg.bam 
samtools cat -o $PREFIX\_R2.uniq.filtered.bam $PREFIX\_tmp/$PREFIX\_R2.uniq.filtered.pos.bam $PREFIX\_tmp/$PREFIX\_R2.uniq.filtered.neg.bam 

#$SAMTOOLS view -b -L $CUT_INTERVAL $PREFIX\_R1.bam > $PREFIX\_R1.filtered.bam
#sleep 10m
#$SAMTOOLS view -b -L $CUT_INTERVAL $PREFIX\_R2.bam > $PREFIX\_R2.filtered.bam
#
#$SAMTOOLS flagstat $PREFIX\_R1.filtered.bam > $PREFIX\_R1.filtered.bam.flagstat &
#$SAMTOOLS flagstat $PREFIX\_R2.filtered.bam > $PREFIX\_R2.filtered.bam.flagstat &
#
#echo "Step3. Merge two mates"  >> hicmap.log
#$SAMTOOLS cat -o $PREFIX.merged.bam $PREFIX\_R1.filtered.bam $PREFIX\_R2.filtered.bam
#
#echo "Step4. Sort by read names"  >> hicmap.log
#$SAMTOOLS sort -m 10G -n $PREFIX.merged.bam $PREFIX.merged.sorted
#rm $PREFIX.merged.bam 
#
#echo "Step5. Pair two mates up and filter read pairs with small distance"  >> hicmap.log
#$SAMTOOLS view $PREFIX.merged.sorted.bam | python $PAIR2MATES $MINLOOPSIZE - > $PREFIX.merged.sorted.paired.sam
#
#echo "Step6. Convert to bam file and sort by genomic coordiantes"  >> hicmap.log
#$SAMTOOLS view -bS -ht $FAI $PREFIX.merged.sorted.paired.sam | \
#	$SAMTOOLS sort -m 10G - $PREFIX.merged.sorted.paired.sorted
#
#$SAMTOOLS flagstat $PREFIX.merged.sorted.paired.sorted.bam > $PREFIX.merged.sorted.paired.sorted.bam.flagstat &
#
#echo "Step7. Filter PCR duplication"  >> hicmap.log
#java -Xmx10g -jar $MARK_DUPLICATE INPUT=$PREFIX.merged.sorted.paired.sorted.bam OUTPUT=$PREFIX.merged.sorted.paired.sorted.nodup.bam ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE=metrics.$PREFIX.txt TMP_DIR=$PREFIX\_tmp
#
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