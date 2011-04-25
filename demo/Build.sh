#!/bin/sh

# convert the reference sequence to our binary format
../bin/MosaikBuild -fr reference/c.elegans_chr2.fasta -oa reference/c.elegans_chr2.dat

# convert the reads to our binary format
../bin/MosaikBuild -q fastq/c_elegans_chr2_test.fastq -out sequence_archives/c_elegans_chr2_test.mkb -st illumina
