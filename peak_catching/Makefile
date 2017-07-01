.PHONY : chipenrich_ldefs
chipenrich_ldefs : nearest_tss.ldef.gz 5kb.ldef.gz 5kb_outside.ldef.gz
5kb.ldef.gz :
	Rscript extract_chipenrich_locusdefs.R
5kb_outside.ldef.gz : 5kb.ldef.gz
nearest_tss.ldef.gz : 5kb.ldef.gz

.PHONY : combinations
combinations : combinations.txt
combinations.txt :
	Rscript create_combinations.R

peak_catching : combinations
	bash commands

summarize_catching:
	Rscript summarize_peak_catching.R
