peak_catching :
	bash commands

count_peaks : num_peaks.txt
num_peaks.txt :
	bash count_peaks.sh

genome_coverage : genome_coverage.txt
genome_coverage.txt : 
	bash genome_coverage.sh

summarize_catching:
	Rscript peak_catching_with_ldefs.R
