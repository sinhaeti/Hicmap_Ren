##Get Started     
```bash
$ git clone --recursive --depth=1 https://github.com/r3fang/Hicmap_Ren.git
$ cd Hicmap_Ren
$ bash install.sh
$ # download demo HiC data (Dixon 2016)
$ wget http://enhancer.sdsc.edu/r3fang/Hicmap_Ren/demo_R1.fastq.bz2
$ wget http://enhancer.sdsc.edu/r3fang/Hicmap_Ren/demo_R2.fastq.bz2
$ wget http://enhancer.sdsc.edu/r3fang/Hicmap_Ren/cutter_sites.tar.gz 
$ tar -xvzf cutter_sites.tar.gz
$ mkdir BWAIndex/; wget http://enhancer.sdsc.edu/r3fang/Hicmap_Ren/mm9.fa -P BWAIndex
$ bwa index BWAIndex/mm9.fa
$ hicmap -t 2 -m 3G -f demo_R1.fastq.bz2 -r demo_R2.fastq.bz2 -n demo \
  -p Picard/MarkDuplicates.jar -g BWAIndex/mm9.fa -c cutter_sites/mm9.Hind3.1000bp.bed -d 1000 
```

##Depedency
- [bwa](https://github.com/lh3/bwa)
- [samtools 1.2+](http://www.htslib.org/doc/samtools.html)
- [Python 2.7](https://www.python.org/download/releases/2.7/)
- [numpy](http://www.numpy.org/)
- [pysam](https://github.com/pysam-developers/pysam)

##Introduction

**Hicmap** is a simple but efficient and fast pipeline for analyzing Hi-C data which includes alignment, reads filtering and heatmap generation.

```
$ bash bin/hicmap

Program: hicmap (Hi-C analysis pipeline)
Version: 12.05-r15
Contact: Rongxin Fang <r3fang@ucsd.edu>
Lab: http://bioinformatics-renlab.ucsd.edu/rentrac/

Map and processing Hi-C reads
(1) indepedently map FASTQ1/FASTQ2 using BWA;
(2) filter reads with MAPQ < 10;
(3) filter reads far from restriction enzyme cutter sites (500bp upstream for + strand/500 downstream for - strand);
(4) pair up two ends and filter invalid hic pairs;
(5) remove PCR duplication using Picard - MarkDuplicates;

usage: hicmap [-h] [-t THREADS] [-m MAX_MEM] [-f FASTQ1] [-r FASTQ2] [-p MarkDuplicates.jar] [-n PREFIX] [-g BWA_GENOME] [-c CUT_SITES] [-m MIN_INSERT_SIZE]

Example:
hicmap -t 20 -m 8G -f data/JL_H4_R1.fastq.bz2 -r data/JL_H4_R2.fastq.bz2 -p bin/MarkDuplicates.jar -n JL_H4 -g BWAIndex/genome.fa -c data/mm9.MboI.500bp -d 1000

Options:    
	-h, --help			 show this help message and exit.
	-t  THREADS			 threads [1].
	-m  MAX_MEM			 max memory usage [4G].
	-f  FASTQ1			 first mate of pair-end sequencing data [.fq/.fastq/.gz/.bz2].
	-r  FASTQ2			 second mate of pair-end sequencing data [.fq/.fastq/.gz/.bz2].
	-p  MARK_DUPLICATE   path to picard MarkDuplicates.jar [bin/MarkDuplicates.jar].
	-n  NAME			 prefix of output files.
	-g  BWA_GENOME		 BWA indexed reference genome.
	-c  CUT_ENZ			 restriction cutting enzyme files. 
	-d  MIN_INSERT_SIZE	 min insert size for valid "DIFFERENT-STRAND" pairs.
```

## Workflow

![workflow](https://github.com/r3fang/Hicmap_Ren/blob/master/img/workflow.jpg)

## FAQ

 1. **What is Hi-C?**     
 2. **Why BWA is chosen as default aligner?**
 3. **What other aligners could be used for aligning Hi-C reads?**
 4. **What does each parameter mean?**
 5. **How does hicmap deal with chimeric reads?**
 6. **What is valid pairs for Hi-C?** 

##Version     
12.05-15

##Author     
Rongxin Fang    
r3fang@ucsd.edu
