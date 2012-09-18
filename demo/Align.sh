#!/bin/sh

### create ../bin/MosaikAligner
# cd ../src
# make
###

# align the reads
ANN_PATH=../src/networkFile
../bin/MosaikAligner -in fastq/read.mkb -out fastq/read.mka -ia reference/e.coli.dat -annpe $ANN_PATH/2.1.26.pe.100.0065.ann -annse $ANN_PATH/2.1.26.se.100.005.ann
# you may need to use the jump database for large genome (> 100 million basepair)
#../bin/MosaikAligner -in fastq/read.mkb -out fastq/read.mka -ia reference/e.coli.dat -annpe $ANN_PATH/2.1.26.pe.100.0065.ann -annse $ANN_PATH/2.1.26.se.100.005.ann -j reference/e.coli.15

