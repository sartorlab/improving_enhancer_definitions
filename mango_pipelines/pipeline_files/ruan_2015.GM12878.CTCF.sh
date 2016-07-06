#!/bin/bash

# SRX1210750
# http://www.ebi.ac.uk/ena/data/view/SRP063492

# Fetch the data
#wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR231/006/SRR2312566/SRR2312566_1.fastq.gz
#wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR231/006/SRR2312566/SRR2312566_2.fastq.gz

#gunzip SRR2312566_1.fastq.gz
#gunzip SRR2312566_2.fastq.gz

mkdir GM12878.CTCF.new_out

Rscript /users/porchard/mango/mango/mango.R --fastq1 SRR2312566_1.fastq --fastq2 SRR2312566_2.fastq --prefix ruan.CTCF.GM12878.new --chrominclude chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 --keepempty TRUE --shortreads FALSE --maxlength 1000 --linkerA CGCGATATCTTATCTGACT --linkerB GTCAGATAAGATATCGCGT --bowtieref bowtie_index/hg19 --bedtoolsgenome ref/hg19.chrom.sizes --outdir GM12878.CTCF.new_out --FDR 0.05 &> GM12878.CTCF.new.out &
