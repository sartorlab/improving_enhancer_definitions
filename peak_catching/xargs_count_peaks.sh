#!/bin/bash
set -e
set -u
set -o pipefail

peak=$1
basepeak=`basename ${peak} '.narrowPeak'`
outfile=${basepeak}.peak_counts

NUM_PEAKS=`wc -l ${peak} | cut -f 1 -d ' '`
echo ${basepeak} ${NUM_PEAKS} > ${outfile}
