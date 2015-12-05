#!/bin/bash

# test if softwares installed
command -v bwa >/dev/null 2>&1 || { echo >&2 "I require bwa but it's not installed.  Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "I require samtools but it's not installed.  Aborting."; exit 1; }

# pass paramters
usage(){
cat << EOF
usage: ${0##*/} [-h] [-f FASTQ1] [-r FASTQ2] [-g BWA_GENOME] [-c CUT_SITES] [-m MIN_INSERT_SIZE]

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

Options:    
	-h, --help      	show this help message and exit.
	-f FASTQ1       	First mate of pair-end sequencing data.
	-r FASTQ1       	Second mate of pair-end sequencing data.
	-g BWA_GENOME   	BWA indexed reference genome.
	-c CUT_ENZ      	Restriction cutting enzyme files. 
	-m MIN_INSERT_SIZE	Min insert size for valid "DIFFERENT-STRAND" pairs.
EOF
} 

while getopts ":f:r:g:c:m:" opt;
do
	case "$opt" in
		f) FASTQ1=$OPTARG;;
		r) FASTQ2=$OPTARG;;
		g) GENOME=$OPTARG;;
		c) CUT_ENZ=$OPTARG;;
		m) MIN_INSERT_SIZE=$OPTARG;;
		\?) usage
			exit 1
			;;
	esac
done

if [ $# -lt 12 ] ; then
   usage
   echo "error: too few arguments, at least 3 arguments"
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

if [ ! -f $CUT_ENZ ]; then
	usage
    echo "error: '$CUT_ENZ' not exists.";
	exit 1
fi

# check if input type
re='^[0-9]+$'
if ! [[ $MIN_INSERT_SIZE =~ $re ]] ; then
   echo "error: '$MIN_INSERT_SIZE' Not a number" >&2; 
   exit 1
fi

#MINLOOPSIZE=10000
#
#cd $DIR
#
#if [ "$CUT" == "Hind3" ]
#then
#	CUT_INTERVAL=$HIND3
#fi
#if [ "$CUT" == "MboI" ]
#then
#	CUT_INTERVAL=$MBOI
#fi
#
#echo "DIR="$DIR > hicmap.log
#echo "PREFIX="$PREFIX >> hicmap.log
#echo "MINLOOPSIZE="$MINLOOPSIZE >> hicmap.log
#echo "MAPPER="$BWA >> hicmap.log
#echo "REFERENCE="$MM9_BWA >> hicmap.log
#echo "CUT="$CUT_INTERVAL >> hicmap.log
#echo >> hicmap.log
#
#echo "Step1. Mapping reads and filter non-uniquely and secondary alignment" >> hicmap.log
#$BWA mem -t $THREADS $MM9_BWA <(bzip2 -dc $FASTQ1) | $SAMTOOLS view -F 2048 -q 10 -bS - > $PREFIX\_R1.bam
#$BWA mem -t $THREADS $MM9_BWA <(bzip2 -dc $FASTQ2) | $SAMTOOLS view -F 2048 -q 10 -bS - > $PREFIX\_R2.bam
#
#$SAMTOOLS flagstat $PREFIX\_R1.bam > $PREFIX\_R1.bam.flagstat &
#$SAMTOOLS flagstat $PREFIX\_R2.bam > $PREFIX\_R2.bam.flagstat &
#
#echo "Step2. Filter reads that are far from restriction cutter sites"  >> hicmap.log
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
#
#
#