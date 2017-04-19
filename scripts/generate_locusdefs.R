chromhmm_pool = c(NA, 'chromhmm')
dnase_pool = c(NA, 'dnase')
thurman_pool = c(NA, 'thurman')
fantom_pool = c(NA, 'fantom')

expand_pool = c(0, 1000)

P2P_pool = c(NA, 'P2P')
E_pool = c(NA, 'E')


combinations = expand.grid(chromhmm_pool, dnase_pool, thurman_pool, fantom_pool, expand_pool, P2P_pool, E_pool, thurman_pool, fantom_pool, stringsAsFactors=F)
# Get rid of 4 NAs in 1-4 and 6-9
combinations = combinations[!apply(combinations, 1, function(row){all(is.na(row[1:4]))}),]
combinations = combinations[!apply(combinations, 1, function(row){all(is.na(row[6:9]))}),]

# Make codes
enhancer_codes = apply(combinations[,1:4], 1, function(row){paste(row[!is.na(row)], collapse='_')})
extension_codes = combinations[,5]

enhancer_base = unique(paste(enhancer_codes, extension_codes, sep='.'))

# Go by enhancer base to minimize the number of file reads
for(base in enhancer_base) {
	# Read files
	message('Reading P2P...')
	base_p2p = read.table(file = sprintf('%s/%s.P2P.pairs', '../pair_lists', base), sep='\t', header = T, as.is = T)
	message('Reading E...')
	base_e = read.table(file = sprintf('%s/%s.E.pairs', '../pair_lists', base), sep='\t', header = T, as.is = T)
	message('Reading thurman...')
	base_thurman = read.table(file = sprintf('%s/%s.thurman.pairs', '../pair_lists', base), sep='\t', header = T, as.is = T)
	message('Reading fantom...')
	base_fantom = read.table(file = sprintf('%s/%s.fantom.pairs', '../pair_lists', base), sep='\t', header = T, as.is = T)

	# Make GRanges
	base_p2p_gr = GenomicRanges::makeGRangesFromDataFrame(df = base_p2p, keep.extra.columns = TRUE)
	base_e_gr = GenomicRanges::makeGRangesFromDataFrame(df = base_e, keep.extra.columns = TRUE)
	base_thurman_gr = GenomicRanges::makeGRangesFromDataFrame(df = base_thurman, keep.extra.columns = TRUE)
	base_fantom_gr = GenomicRanges::makeGRangesFromDataFrame(df = base_fantom, keep.extra.columns = TRUE)

	P2P_pool = c(NA, 'P2P')
	E_pool = c(NA, 'E')
	thurman_pool = c(NA, 'thurman')
	fantom_pool = c(NA, 'fantom')
	interaction_combinations = expand.grid(P2P_pool, E_pool, thurman_pool, fantom_pool, stringsAsFactors=F)
	interaction_combinations = interaction_combinations[!apply(interaction_combinations, 1, function(row){all(is.na(row))}),]

	# Go through each interaction combination and create the locus definition
	for(i in 1:nrow(interaction_combinations)) {
		E_code = interaction_combinations[i,2]
		P2P_code = interaction_combinations[i,1]
		thurman_code = interaction_combinations[i,3]
		fantom_code = interaction_combinations[i,4]

		ldef_gr = GenomicRanges::GRanges()

		interaction_row = interaction_combinations[i,]
		interaction_code = paste(interaction_row[!is.na(interaction_row)], collapse='_')

		ldef_file = paste(base, interaction_code, 'ldef', sep='.')

		if(file.exists(ldef_file)) {
			message(sprintf('%s, skipping...', ldef_file))
			next
		}

	  if (!is.na(P2P_code)) {
	    message('Adding P2P')
	    ldef_gr = c(ldef_gr, base_p2p_gr)
	  }

	  if (!is.na(E_code)) {
	    message('Adding E')
	    ldef_gr = c(ldef_gr, base_e_gr)
	  }

	  if (!is.na(thurman_code)) {
	    message('Adding thurman')
	    ldef_gr = c(ldef_gr, base_fantom_gr)
	  }

	  if (!is.na(fantom_code)) {
	    message('Adding fantom')
	    ldef_gr = c(ldef_gr, base_thurman_gr)
	  }

		message('Reducing...')
		ldef_gr = sort(ldef_gr)
		ldef_grl = IRanges::splitAsList(ldef_gr, ldef_gr$gene_id)
		ldef_grl = GenomicRanges::reduce(ldef_grl)
		ldef_grl_genes = cds_txname_rle = S4Vectors::Rle(names(ldef_grl), S4Vectors::elementNROWS(ldef_grl))

		ldef_gr = unlist(ldef_grl, use.names = FALSE)
		ldef_gr$gene_id = ldef_grl_genes

		ldef_gr = sort(ldef_gr)

		ldef_df = BiocGenerics::as.data.frame(ldef_gr)
		colnames(ldef_df) = c('chr','start','end','width','strand','gene_id')
		ldef_df = ldef_df[, c('chr','start','end','gene_id')]

		message(sprintf('Writing %s...', ldef_file))
		write.table(ldef_df, file = ldef_file, sep = '\t', quote = F, row.names = F, col.names = T)
	}
}
