### For data downloads:

download_data : dnase chromhmm thurman fantom ldef motif

dnase : ./data/ENCODE_dnase/wgEncodeAwgDnaseMasterSites.bed.gz
./data/ENCODE_dnase/wgEncodeAwgDnaseMasterSites.bed.gz :
	mkdir -p ./data/ENCODE_dnase && cd ./data/ENCODE_dnase && bash ../../scripts/public_data_download/DNase_seq_download.sh

# List them all, but for the make rule, use the last one as a way to make sure you got the previous ones
# Not ideal, but what are you going to do?
CHROMHMM_FILES := wgEncodeBroadHmmGm12878HMM.bed.gz wgEncodeBroadHmmH1hescHMM.bed.gz wgEncodeBroadHmmHepg2HMM.bed.gz wgEncodeBroadHmmHmecHMM.bed.gz wgEncodeBroadHmmHsmmHMM.bed.gz wgEncodeBroadHmmHuvecHMM.bed.gz wgEncodeBroadHmmK562HMM.bed.gz wgEncodeBroadHmmNhekHMM.bed.gz wgEncodeBroadHmmNhlfHMM.bed.gz

chromhmm : ./data/ENCODE_chromhmm/chromhmm_combined.bed.gz
./data/ENCODE_chromhmm/chromhmm_combined.bed.gz :
	mkdir -p ./data/ENCODE_chromhmm && cd ./data/ENCODE_chromhmm && bash ../../scripts/public_data_download/ChromHMM_download.sh

thurman : ./data/thurman_dnase/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz
./data/thurman_dnase/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz :
	mkdir -p ./data/thurman_dnase && cd ./data/thurman_dnase && bash ../../scripts/public_data_download/thurman_download.sh

fantom : ./data/fantom/enhancer_tss_associations.bed
./data/fantom/enhancer_tss_associations.bed :
	mkdir -p ./data/fantom && cd ./data/fantom && bash ../../scripts/public_data_download/fantom_download.sh

ldef : ./data/genes/chipenrich_5kb_locusdef.txt
./data/genes/chipenrich_5kb_locusdef.txt :
	mkdir -p ./data/genes && cd ./data/genes && Rscript ../../scripts/extract_5kb_ldef.R

# Use what Heming had downloaded from what was supposed to be the same location
# NOTE: Her version does not match the website http://jaspar.genereg.net/?ID=MA0139.1&rm=present&collection=CORE
# exactly for some reason, and downloading it fresh does not work because the column
# sums are not all the same. NO CLUE!
# motif : ./data/motif/MA0139.1.pfm
# ./data/motif/MA0139.1.pfm :
# 	mkdir -p ./data/motif/ && cd ./data/motif && bash ../../scripts/public_data_download/CTCF_motif_download.sh
