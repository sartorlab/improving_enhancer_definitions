#!/bin/bash

curl http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/individual/MA0139.1.pfm > MA0139.1.pfm

awk 'NR > 1 {print $0}' MA0139.1.pfm > tmp.pfm
sed -ie 's/\[//g' tmp.pfm
sed -ie 's/\]//g' tmp.pfm
tr -s " " < tmp.pfm > ctcf_motif.pfm
