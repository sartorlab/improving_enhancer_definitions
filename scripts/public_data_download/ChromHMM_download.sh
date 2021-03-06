#!/bin/bash

# ChromHMM
wget -nv http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmGm12878HMM.bed.gz
wget -nv http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmH1hescHMM.bed.gz
wget -nv http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHepg2HMM.bed.gz
wget -nv http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHmecHMM.bed.gz
wget -nv http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHsmmHMM.bed.gz
wget -nv http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHuvecHMM.bed.gz
wget -nv http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmK562HMM.bed.gz
wget -nv http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmNhekHMM.bed.gz
wget -nv http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmNhlfHMM.bed.gz

cat wgEncodeBroadHmmGm12878HMM.bed.gz wgEncodeBroadHmmH1hescHMM.bed.gz wgEncodeBroadHmmHepg2HMM.bed.gz wgEncodeBroadHmmHmecHMM.bed.gz wgEncodeBroadHmmHsmmHMM.bed.gz wgEncodeBroadHmmHuvecHMM.bed.gz wgEncodeBroadHmmK562HMM.bed.gz wgEncodeBroadHmmNhekHMM.bed.gz wgEncodeBroadHmmNhlfHMM.bed.gz > chromhmm_combined.bed.gz

rm wgEncodeBroadHmmGm12878HMM.bed.gz wgEncodeBroadHmmH1hescHMM.bed.gz wgEncodeBroadHmmHepg2HMM.bed.gz wgEncodeBroadHmmHmecHMM.bed.gz wgEncodeBroadHmmHsmmHMM.bed.gz wgEncodeBroadHmmHuvecHMM.bed.gz wgEncodeBroadHmmK562HMM.bed.gz wgEncodeBroadHmmNhekHMM.bed.gz wgEncodeBroadHmmNhlfHMM.bed.gz
