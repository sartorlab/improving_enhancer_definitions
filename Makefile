### For data downloads:

chipseq:
	@bash ../scripts/public_data_download/ENCODE_chipseq.sh

dnase:
	@mkdir ENCODE_dnase && cd ENCODE_dnase && bash ../scripts/public_data_download/DNase_seq_download.sh

chromhmm:
	@mkdir ENCODE_chromhmm && cd ENCODE_chromhmm && bash ../scripts/public_data_download/ChromHMM_download.sh

