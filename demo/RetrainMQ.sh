#!/bin/bash

### define your bamtools
# https://github.com/pezmaster31/bamtools
BAMTOOLS=/share/software/bamtools/bin/bamtools
###

### create ../bin/MosaikBuild and ../bin/MosaikAligner
# cd ../src
# make
###

### create executable files
RETRAIN_PATH=../src/networkFile/retrainCode
# cd $RETRAIN_PATH
# make
###

# MosaikBuild
sh Build.sh

# MosaikAlign
ANN_PREFIX=../src/networkFile/2.1.26.
../bin/MosaikAligner -in fastq/read.mkb -out fastq/read.mka -ia reference/e.coli.dat -annpe $ANN_PREFIX\pe.100.0065.ann -annse $ANN_PREFIX\se.100.005.ann -zn

# sort bam by read names and convert it to
BAM_PREFIX=fastq/read.mka
$BAMTOOLS sort -byname -in $BAM_PREFIX.bam \
  | $BAMTOOLS convert -format sam -noheader -out $BAM_PREFIX.sorted.sam

# attach ZC
$RETRAIN_PATH/attachXC/xc_pe fastq/gold.sam $BAM_PREFIX.sorted.sam > $BAM_PREFIX.xc.sam

# train neural network
FRAGMENT=1000
$RETRAIN_PATH/trainNetwork/train_mq -i $BAM_PREFIX.xc.sam -o $BAM_PREFIX.pe.ann -p -f $FRAGMENT -e 0.0015
