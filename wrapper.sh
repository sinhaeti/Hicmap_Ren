#!/bin/bash
#PBS -q hotel 
#PBS -N hicmap
#PBS -l nodes=1:ppn=5
#PBS -l walltime=100:00:00
#PBS -o hicmap.out
#PBS -e hicmap.err
#PBS -V
#PBS -M r3fang@ucsd.edu
#PBS -m ae
#PBS -A ren-group


MARK_DUPLICATE="/oasis/tscc/scratch/r3fang/github/Hicmap_Ren/Picard/MarkDuplicates.jar"

MM9_BWA="/oasis/tscc/scratch/r3fang/data/Mus_musculus/UCSC/mm9/Sequence/BWAIndex/genome.fa"
HG19_BWA="/oasis/tscc/scratch/r3fang/data/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"

HG19_HIND3="/oasis/tscc/scratch/r3fang/github/Hicmap_Ren/cutter_sites/hg19.Hind3.1000bp.bed"
MM9_HIND3="/oasis/tscc/scratch/r3fang/github/Hicmap_Ren/cutter_sites/mm9.Hind3.1000bp.bed"
HG19_MBOI="/oasis/tscc/scratch/r3fang/github/Hicmap_Ren/cutter_sites/hg19.MboI.1000bp.bed"
MM9_MBOI="/oasis/tscc/scratch/r3fang/github/Hicmap_Ren/cutter_sites/mm9.MboI.1000bp.bed"

cd /to/your/current/directory/
hicmap -t 5 -m 2G -f demo_R1.fastq.bz2 -r demo_R2.fastq.bz2 -n demo_ -p $MARK_DUPLICATE -g $MM9_BWA -c $MM9_MBOI -d 10000


