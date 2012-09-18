#!/bin/sh

### Create ../bin/MosaikBuild
# cd ../src
# make
###

# Convert the reference sequence to our binary format
../bin/MosaikBuild -fr reference/e.coli.fa -oa reference/e.coli.dat
# You may need to create the jump database for large genome (> 100 million basepair)
#../bin/MosaikJump -ia reference/e.coli.dat -hs 15 -out reference/e.coli.15

# Convert the reads to our binary format
../bin/MosaikBuild -q fastq/mate1.fq -q2 fastq/mate2.fq -out fastq/read.mkb -st illumina
