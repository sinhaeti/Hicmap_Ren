#!/bin/bash
#PBS -q hotel 
#PBS -N Wrapper
#PBS -l nodes=1:ppn=4
#PBS -l walltime=100:00:00
#PBS -o Wrapper.out
#PBS -e Wrapper.err
#PBS -V
#PBS -M r3fang@ucsd.edu
#PBS -m ae
#PBS -A ren-group

# tools 
MARK_DUPLICATE="/oasis/tscc/scratch/r3fang/software/picard-tools-1.49/MarkDuplicates.jar"

MM9_BWA="/oasis/tscc.old/scratch/r3fang/data/Mus_musculus/UCSC/mm9/Sequence/BWAIndex/genome.fa"

MM9_HIND3="/oasis/tscc/scratch/r3fang/github/Hicmap_Ren/data/mm9.Hind3.500bp"
MM9_MBOI="/oasis/tscc/scratch/r3fang/github/Hicmap_Ren/data/mm9.MboI.500bp"

cd /oasis/tscc/scratch/r3fang/github/Hicmap_Ren
hicmap -t 4 -m 2G -f data/demo_R1.fastq.bz2 -r data/demo_R2.fastq.bz2 -n demo_bwa -p $MARK_DUPLICATE -g $MM9_BWA -c $MM9_MBOI -d 10000

