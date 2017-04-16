#!/bin/bash

wget -nv http://enhancer.binf.ku.dk/presets/permissive_enhancers.bed
wget -nv http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed

awk -v OFS='\t' 'NR > 1 {print $0}' permissive_enhancers.bed > tmp.bed
mv tmp.bed permissive_enhancers.bed

awk -v OFS='\t' 'NR > 2 {print $0}' enhancer_tss_associations.bed > tmp.bed
mv tmp.bed enhancer_tss_associations.bed
