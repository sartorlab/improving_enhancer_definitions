#!/bin/bash
set -e
set -u
set -o pipefail

for file in `ls ~/latte/share/improved_enhancer_definitions/locusdefs/*.ldef`
do
	echo `basename ${file}` `awk 'NR > 1 {SUM += $3 - $2} END {print SUM}' ${file}` >> genome_coverage.txt
done
