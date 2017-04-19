#!/bin/bash
set -e
set -u
set -o pipefail

THURMAN_DATA=../data/thurman_dnase/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz

gunzip -c ${THURMAN_DATA} | cut -f 1-3,5-7 > thurman.final_interactions
