#!/bin/sh

### create ../bin/MosaikAligner
# cd ../src
# make
###

# align the reads
../bin/MosaikAligner -in sequence_archives/c_elegans_chr2_test.mkb -out sequence_archives/c_elegans_chr2_test.mka -ia reference/c.elegans_chr2.dat

