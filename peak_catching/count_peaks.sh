#!/bin/bash
set -e
set -u
set -o pipefail

for file in `ls ~/espresso/share/ENCODE_chipseq/*.narrowPeak`
do
	echo `basename ${file}` `wc -l ${file} | cut -f 1 -d ' '` >> num_peaks.txt
done