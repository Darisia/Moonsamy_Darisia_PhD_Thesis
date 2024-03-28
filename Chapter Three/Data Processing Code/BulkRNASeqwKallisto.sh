#!/bin/bash

# This bash script contains the code to pseudoalign the bulk RNA sequences to the reference transcriptome



## PRJNA694583 - THP-1 stimulated with IFN-alpha 
# H1 (THP1 - control 1) - 
kallisto quant -i ./kallisto_index.idx \
    -o ./6.0.BulkRNASeq/24011103/24011105 \
    -t 2 \
    ../Raw/PRJNA694583/SRR13520435_1.fastq.gz ../Raw/PRJNA694583/SRR13520435_2.fastq.gz

# H2 (THP1 - control 2) - 
kallisto quant -i ./kallisto_index.idx \
    -o ./6.0.BulkRNASeq/24011103/24011106 \
    -t 2 \
    ../Raw/PRJNA694583/SRR13520436_1.fastq.gz ../Raw/PRJNA694583/SRR13520436_2.fastq.gz

# H3 (THP1 - control 3) - 
kallisto quant -i ./kallisto_index.idx \
    -o ./6.0.BulkRNASeq/24011103/24011107 \
    -t 2 \
    ../Raw/PRJNA694583/SRR13520437_1.fastq.gz ../Raw/PRJNA694583/SRR13520437_2.fastq.gz

# IFNA-1 (THP-1)
kallisto quant -i ./kallisto_index.idx \
    -o ./6.0.BulkRNASeq/24011103/24011108 \
    -t 2 \
    ../Raw/PRJNA694583/SRR13520438_1.fastq.gz ../Raw/PRJNA694583/SRR13520438_2.fastq.gz

# IFNA-2 (THP-1)
kallisto quant -i ./kallisto_index.idx \
    -o ./6.0.BulkRNASeq/24011103/24011109 \
    -t 2 \
    ../Raw/PRJNA694583/SRR13520439_1.fastq.gz ../Raw/PRJNA694583/SRR13520439_2.fastq.gz

# IFNA-3 (THP-1)
kallisto quant -i ./kallisto_index.idx \
    -o ./6.0.BulkRNASeq/24011103/24011110 \
    -t 2 \
    ../Raw/PRJNA694583/SRR13520440_1.fastq.gz ../Raw/PRJNA694583/SRR13520440_2.fastq.gz


## PRJNA922988 - STAT2 Knockout RNA-seq
kallisto quant -i ./kallisto_index.idx -o ./23112102/23112103/ --single -l 75 -s 1 ./23112102/SRR22897427.fastq.gz
kallisto quant -i ./kallisto_index.idx -o ./23112102/23112104/ --single -l 75 -s 1 ./23112102/SRR22897426.fastq.gz
kallisto quant -i ./kallisto_index.idx -o ./23112102/23112105/ --single -l 75 -s 1 ./23112102/SRR22897423.fastq.gz
kallisto quant -i ./kallisto_index.idx -o ./23112102/23112106/ --single -l 75 -s 1 ./23112102/SRR22897422.fastq.gz

## PRJNA922988 - IRF9 Knockout RNA-seq
kallisto quant -i ./kallisto_index.idx -o ./23112107/23112108/ --single -l 75 -s 1 ./23112107/SRR22897467.fastq.gz
kallisto quant -i ./kallisto_index.idx -o ./23112107/23112109/ --single -l 75 -s 1 ./23112107/SRR22897466.fastq.gz
kallisto quant -i ./kallisto_index.idx -o ./23112107/23112110/ --single -l 75 -s 1 ./23112107/SRR22897457.fastq.gz
kallisto quant -i ./kallisto_index.idx -o ./23112107/23112111/ --single -l 75 -s 1 ./23112107/SRR22897456.fastq.gz


# Note - all bulk RNA-seq pseudoalignments were run as above (they were run exactly as described in the chapter three methods)