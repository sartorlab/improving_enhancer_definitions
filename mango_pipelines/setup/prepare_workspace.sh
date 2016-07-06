#!/bin/bash

ROOT_DIR=""

cd $ROOT_DIR

# prepare workspace for running mango
mkdir ref
mkdir data
mkdir bowtie_index
mkdir mango_out
mkdir data

# download hg19
cd ref
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*
cd ..

# download the bowtie2 reference rather than building ourselves...
cd bowtie_index
wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19.ebwt.zip
unzip bowtie_index/hg19.ebwt.zip
cd ..

# download a bedtools genome file
cd ref
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
cd ..
