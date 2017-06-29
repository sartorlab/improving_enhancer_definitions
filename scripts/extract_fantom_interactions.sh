#!/bin/bash
set -e
set -u
set -o pipefail

FANTOM_DATA=../data/fantom/enhancer_tss_associations.bed

Rscript ../scripts/extract_fantom_interactions.R --fantom_interactions ${FANTOM_DATA}
