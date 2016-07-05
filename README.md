##Get Started     
```bash
$ git clone --depth=1 https://github.com/r3fang/Hicmap_Ren.git
$ cd Hicmap_Ren
$ bash install.sh
$ # download demo HiC data (Dixon 2016)
$ wget http://enhancer.sdsc.edu/r3fang/Hicmap_Ren/demo_R1.fastq.bz2
$ wget http://enhancer.sdsc.edu/r3fang/Hicmap_Ren/demo_R2.fastq.bz2
$ wget http://enhancer.sdsc.edu/r3fang/Hicmap_Ren/cutter_sites.tar.gz 
$ tar -xvzf cutter_sites.tar.gz
$ mkdir BWAIndex/; wget http://enhancer.sdsc.edu/r3fang/Hicmap_Ren/mm9.fa -P BWAIndex
$ bwa index BWAIndex/mm9.fa
$ hicmap -f demo_R1.fastq.bz2 -r demo_R2.fastq.bz2 -n demo -p Picard/MarkDuplicates.jar \
  -g BWAIndex/mm9.fa -c cutter_sites/mm9.Hind3.1000bp.bed
```

##Depedency
- [bwa](https://github.com/lh3/bwa)
- [samtools 1.2+](http://www.htslib.org/doc/samtools.html)
- [Python 2.7](https://www.python.org/download/releases/2.7/)

##Introduction

**Hicmap** is a simple but efficient and fast pipeline for analyzing Hi-C data which includes alignment, reads filtering and heatmap generation.

```
$ bash bin/hicmap

Program: hicmap (Hi-C analysis pipeline)
Version: v07.04.2016
Contact: Rongxin Fang <r3fang@ucsd.edu>
Lab: http://bioinformatics-renlab.ucsd.edu/rentrac/

usage: hicmap [-h] [-t THREADS] [-m MAX_MEM] [-f FASTQ1] [-r FASTQ2] [-p MarkDuplicates.jar] [-n PREFIX] [-g BWA_GENOME] [-c CUT_SITES] [-i INSERT_SIZE]

Example:
hicmap -f demo_R1.fastq.bz2 -r demo_R2.fastq.bz2 -p Picard/MarkDuplicates.jar -n demo -g genome.fa -c mm9.MboI.bed

Options:

-- required:
    -f  FASTQ1           first mate of pair-end sequencing data  [.fq|.gz|.bz2].
    -r  FASTQ2           second mate of pair-end sequencing data [.fq|.gz|.bz2].
    -p  MARK_DUPLICATE   path to picard MarkDuplicates.jar [Picard/MarkDuplicates.jar].
	-n  NAME             prefix of output files.
	-g  BWA_GENOME       BWA indexed reference genome.
	-c  CUT_ENZ          restriction cutter sites files.
-- optional:
	-i  INSERT_SIZE      pairs with distance smaller that MIN_INSERT_SIZE will be filtered [10,000].
	-t  INT              threads to be used for BWA mapping [2].
	-m  MAX_MEM          max memory will be used for .bam sorting and MarkDuplicates.jar [2G].

```

##Version     
v07.04.2016

##Author     
Rongxin Fang    
r3fang@ucsd.edu
