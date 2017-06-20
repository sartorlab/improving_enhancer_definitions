#!/bin/bash
set -e
set -u
set -o pipefail

rm genome_coverage.txt

for file in `ls ~/latte/share/improving_enhancer_definitions/locusdefs/*.ldef`
do
	echo `basename ${file}` `awk 'NR > 1 {SUM += $3 - $2} END {print SUM}' ${file}` >> genome_coverage.txt
done

echo `basename ${file}` `awk -v OFS='\t' 'NR > 1 {print $1, $2, $3}' ${file} | sort -T . -k1,1 -k2,2n | uniq | awk -v OFS='\t' '{SUM += $3 - $2} END {print SUM}'` >> genome_coverage.txt
