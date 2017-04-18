#!/bin/bash
set -e
set -u
set -o pipefail

FANTOM_DATA=../data/fantom/enhancer_tss_associations.bed

# Grab first three columns
# Grab the first element of ; delimited 4th column, replace : with -, and then - with tabs
paste <(cut -f 1-3 ${FANTOM_DATA}) <(cut -f 4 ${FANTOM_DATA} | cut -f 1 -d ';' | sed -e 's/:/-/g' | awk '{gsub("-","\t",$0); print}') > fantom.interactions
