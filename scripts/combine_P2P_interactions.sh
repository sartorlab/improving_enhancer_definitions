#!/bin/bash
set -e
set -u
set -o pipefail

cat `ls *.interactions` | sort -T . -k1,1 -k2,2n | cut -f 1-6 > P2P.final_interactions
