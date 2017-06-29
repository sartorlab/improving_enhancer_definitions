suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--fantom_interactions"), action = "store", type = "character", help = "[Required] Path to the enhancer tss association file from FANTOM5 (http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed)")
  )

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)

fantom_raw = read.table(opts$fantom_interactions, header = F, sep='\t', as.is=T)

block_size1 = as.integer(sapply(fantom_raw$V11, function(bs) {unlist(strsplit(bs, ','))[1]}, USE.NAMES = F))
block_size2 = as.integer(sapply(fantom_raw$V11, function(bs) {unlist(strsplit(bs, ','))[2]}, USE.NAMES = F))

fantom_interactions = data.frame(
	'chr1' = fantom_raw$V1,
	'start1' = fantom_raw$V2,
	'end1' = fantom_raw$V2 + block_size1,
	'chr2' = fantom_raw$V1,
	'start2' = fantom_raw$V3 - block_size2,
	'end2' = fantom_raw$V3,
	stringsAsFactors = F)

write.table(fantom_interactions, file = 'fantom.final_interactions', sep='\t', row.names = F, col.names = F, quote = F)
