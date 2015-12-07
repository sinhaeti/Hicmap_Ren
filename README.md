##Get Started     
```
$ git clone https://github.com/r3fang/Hicmap_Ren.git
$ cd Hicmap_Ren
$ bash install.sh
$ PATH=$PATH:/path/to/Hicmap_Ren/bin
$ export PATH
$ hicmap -t 2 -m 3G -d 1000 -f data/JL_H4_R1.fastq.bz2 -r  -n JL_H4 data/JL_H4_R2.fastq.bz2 -p bin/MarkDuplicates.jar -g BWAIndex/genome.fa -c data/mm9.MboI.500bp 
```

##Depedency
- samtools 1.2+
- bwa
- picards (MarkDuplicates.jar)

##Introduction

**Hicmap** is a simple but efficient and fast pipeline for analyzing Hi-C data which includes alignment, reads filtering and heatmap generation.

```
$ bash /bin/hicmap

usage: hicmap.sh [-h] [-t THREADS] [-m MAX_MEM] [-f FASTQ1] [-r FASTQ2] [-p MarkDuplicates.jar] [-n PREFIX] [-g BWA_GENOME] [-c CUT_SITES] [-m MIN_INSERT_SIZE]

Map and processing Hi-C reads
(1) indepedently map FASTQ1/FASTQ2 using BWA;
(2) filter reads with MAPQ < 10;
(3) filter reads far from restriction enzyme cutter sites (500bp upstream for + strand/500 downstream for - strand);
(4) pair up two ends and filter invalid hic pairs;
(5) Remove PCR duplication using Picard - MarkDuplicates;

Example:
	 hicmap -t 2 -m 3G -d 1000 -f data/JL_H4_R1.fastq.bz2 -r  -n JL_H4 data/JL_H4_R2.fastq.bz2 -p bin/MarkDuplicates.jar -g BWAIndex/genome.fa -c data/mm9.MboI.500bp 

Options:    
	-h, --help			show this help message and exit.
	-t  THREADS			threads [1].
	-m  MAX_MEM			max memory usage [4G].
	-f  FASTQ1			first mate of pair-end sequencing data [.fq/.fastq/.gz/.bz2].
	-r  FASTQ2			second mate of pair-end sequencing data [.fq/.fastq/.gz/.bz2].
	-p  MARK_DUPLICATE  path to picard MarkDuplicates.jar
	-n  NAME			prefix of output files.
	-g  BWA_GENOME		BWA indexed reference genome.
	-c  CUT_ENZ			restriction cutting enzyme files. 
	-d  MIN_INSERT_SIZE	min insert size for valid "DIFFERENT-STRAND" pairs.
```

##Version     
12.05-15

##Author     
Rongxin Fang    
r3fang@ucsd.edu
