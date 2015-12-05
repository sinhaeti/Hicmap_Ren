## Hicmap_Ren
**hicmap** is a fast pipeline for analyzing Hi-C data which includes alignment, reads filtering.

##Depedency
- samtools 1.2+
- bwa
- Picards

##Introduction

```
$./bin/hicmap

usage: hicmap [-h] [-f FASTQ1] [-r FASTQ2] [-n PREFIX] [-g BWA_GENOME] [-c CUT_SITES] [-m MIN_INSERT_SIZE]

Map and processing Hi-C reads
(1) map FASTQ1/FASTQ2 using BWA indepedently;
(2) filter reads with MAPQ < 10;
(3) filter reads fell >500bp far from restriction enzyme cutter sites (strand sensitive);
(4) sort reads by read names;
(5) pair up two ends and filter invalid hic pairs;
(6) Remove PCR duplication using Picard - markDuplicates;

Example:
	bash bin/hicmap.sh -t 20 -m 8G -f data/JL_H4_R1.fastq.bz2 -r data/JL_H4_R2.fastq.bz2 -n JL_H4 -g /mnt/thumper/home/r3fang/data/Mus_musculus/UCSC/mm9/Sequence/BWAIndex/genome.fa -c data/mm9.MboI.500bp -d 1000

Options:    
	-h, --help			show this help message and exit.
	-t  THREADS			threads [1].
	-m  MAX_MEM			max memory usage [4G].
	-f  FASTQ1			first mate of pair-end sequencing data.
	-r  FASTQ2			second mate of pair-end sequencing data.
	-n  NAME			prefix of output files.
	-g  BWA_GENOME			BWA indexed reference genome.
	-c  CUT_ENZ			restriction cutting enzyme files. 
	-d  MIN_INSERT_SIZE		min insert size for valid "DIFFERENT-STRAND" pairs.
```

##Version     
12.05-15

##Author     
Rongxin Fang    
r3fang@ucsd.edu