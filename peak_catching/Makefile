.PHONY : combinations
combinations : combinations.txt
combinations.txt :
	Rscript create_combinations.R

peak_catching : combinations
	bash commands

count_peaks :
	bash count_peaks.sh

genome_coverage :
	bash genome_coverage.sh

summarize_catching:
	Rscript peak_catching.R
