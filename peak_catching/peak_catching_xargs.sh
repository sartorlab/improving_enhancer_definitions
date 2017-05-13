#!/bin/bash
set -e
set -u
set -o pipefail

input=$1
ldef=`echo ${input} | cut -f 1 -d ' '`
peak=`echo ${input} | cut -f 2 -d ' '`
baseldef=`basename ${ldef} '.ldef'`
basepeak=`basename ${peak} '.narrowPeak'`
outfile=${basepeak}_${baseldef}.counts

# -a is sorted midpoints of the peaks in the .narrowPeak
# -b is the locus definition (skipping the header column)
# The pipe is to grab only the unique peaks (cols 1-3 should be ordered)
NUM_CAUGHT=`bedtools intersect -wa -wb \
    -a <(awk -v OFS='\t' '{print $1, int(($3 + $2)/2), int(($3 + $2)/2)}' ${peak} | sort -T . -k1,1 -k2,2n) \
    -b <(awk -v OFS='\t' 'NR > 1 {print $0}' ${ldef}) \
| cut -f 1-3 | uniq | wc -l`
echo ${basepeak} ${baseldef} ${NUM_CAUGHT} > ${outfile}
