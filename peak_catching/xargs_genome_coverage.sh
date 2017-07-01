#!/bin/bash
set -e
set -u
set -o pipefail

ldef=$1
baseldef=`basename ${ldef} '.ldef.gz'`
outfile=${baseldef}.coverage

COVERAGE=`awk -v OFS='\t' 'NR > 1 {print $1, $2, $3}' <(gunzip -c ${ldef}) | sort -T . -k1,1 -k2,2n | uniq | awk -v OFS='\t' '{SUM += $3 - $2} END {print SUM}'`
echo ${baseldef} ${COVERAGE} > ${outfile}
