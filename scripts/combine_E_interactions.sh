#!/bin/bash
set -e
set -u
set -o pipefail

cat `ls *.interactions.with_motif` | sort -T . -k1,1 -k2,2n > E.final_interactions
