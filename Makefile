### For data downloads:

dnase : ./ENCODE_dnase/wgEncodeAwgDnaseMasterSites.bed.gz
./ENCODE_dnase/wgEncodeAwgDnaseMasterSites.bed.gz :
	@mkdir ENCODE_dnase && cd ENCODE_dnase && bash ../scripts/public_data_download/DNase_seq_download.sh

# List them all, but for the make rule, use the last one as a way to make sure you got the previous ones
# Not ideal, but what are you going to do?
CHROMHMM_FILES := wgEncodeBroadHmmGm12878HMM.bed.gz wgEncodeBroadHmmH1hescHMM.bed.gz wgEncodeBroadHmmHepg2HMM.bed.gz wgEncodeBroadHmmHmecHMM.bed.gz wgEncodeBroadHmmHsmmHMM.bed.gz wgEncodeBroadHmmHuvecHMM.bed.gz wgEncodeBroadHmmK562HMM.bed.gz wgEncodeBroadHmmNhekHMM.bed.gz wgEncodeBroadHmmNhlfHMM.bed.gz

chromhmm : ./ENCODE_chromhmm/wgEncodeBroadHmmNhlfHMM.bed.gz
./ENCODE_chromhmm/wgEncodeBroadHmmNhlfHMM.bed.gz :
	@mkdir ENCODE_chromhmm && cd ENCODE_chromhmm && bash ../scripts/public_data_download/ChromHMM_download.sh

thurman : ./thurman_dnase/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz
./thurman_dnase/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz :
	@mkdir thurman_dnase && cd thurman_dnase && bash ../scripts/public_data_download/thurman_download.sh

fantom : ./fantom/enhancer_tss_associations.bed
./fantom/enhancer_tss_associations.bed :
	@mkdir fantom && cd fantom && bash ../scripts/public_data_download/fantom_download.sh
