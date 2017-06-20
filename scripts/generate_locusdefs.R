library(GenomicRanges)

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
		ldef_grl_genes = S4Vectors::Rle(names(ldef_grl), S4Vectors::elementNROWS(ldef_grl))

		ldef_gr = unlist(ldef_grl, use.names = FALSE)
		ldef_gr$gene_id = ldef_grl_genes

		ldef_gr = sort(ldef_gr)

		ldef_df = BiocGenerics::as.data.frame(ldef_gr)
		colnames(ldef_df) = c('chr','start','end','width','strand','gene_id')
		ldef_df = subset(ldef_df, width != 0)
		ldef_df = ldef_df[, c('chr','start','end','gene_id')]

		message(sprintf('Writing %s...', ldef_file))
		write.table(ldef_df, file = ldef_file, sep = '\t', quote = F, row.names = F, col.names = T)
	}

	dist_within_tss = 500000

	# Go through the enhancer bases again and assign the enhancers to the nearest TSS based on the 5kb definition
	data('locusdef.hg19.5kb', package='chipenrich.data')
	fivekb_gr = locusdef.hg19.5kb@granges

	enh_df = read.table(file = sprintf('%s/%s.enhancers.gz', '../enhancer_lists', base), sep='\t', header = F, as.is = T)
	colnames(enh_df) = c('chr','start','end')
	enh_gr = GenomicRanges::makeGRangesFromDataFrame(df = enh_df)

	ldef_file = paste(base, sprintf('nearest_tss_%s', format(dist_within_tss, scientific = FALSE)), 'ldef', sep='.')

	if(file.exists(ldef_file)) {
		message(sprintf('%s, skipping...', ldef_file))
		next
	}

	# These return Hits objects, but the queryHits() should have unique indices only
	nearests = GenomicRanges::nearest(enh_gr, fivekb_gr, ignore.strand = TRUE, select='all')
	distances = GenomicRanges::distanceToNearest(enh_gr, fivekb_gr, ignore.strand = TRUE, select='all')

	keep = which(GenomicRanges::mcols(distances)$distance < dist_within_tss)

	if(length(queryHits(nearests)) != length(unique(queryHits(nearests)))) {
		stop(sprintf('There are ties for nearest 5kb with the enhancer base: %s', base))
	}

	ldef_gr = enh_gr[queryHits(nearests)[keep]]
	GenomicRanges::mcols(ldef_gr)$gene_id = fivekb_gr[subjectHits(nearests)[keep]]$gene_id

	message('Reducing...')
	ldef_gr = sort(ldef_gr)
	ldef_grl = IRanges::splitAsList(ldef_gr, ldef_gr$gene_id)
	ldef_grl = GenomicRanges::reduce(ldef_grl)
	ldef_grl_genes = S4Vectors::Rle(names(ldef_grl), S4Vectors::elementNROWS(ldef_grl))

	ldef_gr = unlist(ldef_grl, use.names = FALSE)
	ldef_gr$gene_id = ldef_grl_genes

	ldef_gr = sort(ldef_gr)

	ldef_df = BiocGenerics::as.data.frame(ldef_gr)
	colnames(ldef_df) = c('chr','start','end','width','strand','gene_id')
	ldef_df = subset(ldef_df, width != 0)
	ldef_df = ldef_df[, c('chr','start','end','gene_id')]

	message(sprintf('Writing %s...', ldef_file))
	write.table(ldef_df, file = ldef_file, sep = '\t', quote = F, row.names = F, col.names = T)
}

# I just want to put this odd occurrence here just in case anyone makes the mistake
# of thinking that nearest() and distanceToNearest() ignores chromosomes.
# > ldef_gr
# GRanges object with 307773 ranges and 1 metadata column:
#            seqnames                 ranges strand |   gene_id
#               <Rle>              <IRanges>  <Rle> | <integer>
#        [1]     chr1         [41081, 43537]      * |    641702
#        [2]     chr1         [45137, 45537]      * |    641702
#        [3]     chr1         [46337, 47537]      * |    641702
#        [4]     chr1         [52737, 53737]      * |     79501
#        [5]     chr1         [54737, 55137]      * |     79501
#        ...      ...                    ...    ... .       ...
#   [307769]     chrX [155246606, 155247006]      * |      3581
#   [307770]     chrX [155249606, 155251406]      * |      3581
#   [307771]     chrX [155254606, 155256806]      * |      3581
#   [307772]     chrX [155257406, 155258206]      * |      3581
#   [307773]     chrX [155259206, 155259606]      * |      3581
#   -------
#   seqinfo: 23 sequences from an unspecified genome; no seqlengths
# > enh_gr
# GRanges object with 341715 ranges and 0 metadata columns:
#            seqnames                 ranges strand
#               <Rle>              <IRanges>  <Rle>
#        [1]     chr1         [41081, 43537]      *
#        [2]     chr1         [45137, 45537]      *
#        [3]     chr1         [46337, 47537]      *
#        [4]     chr1         [52737, 53737]      *
#        [5]     chr1         [54737, 55137]      *
#        ...      ...                    ...    ...
#   [341711]     chrX [155246606, 155247006]      *
#   [341712]     chrX [155249606, 155251406]      *
#   [341713]     chrX [155254606, 155256806]      *
#   [341714]     chrX [155257406, 155258206]      *
#   [341715]     chrX [155259206, 155259606]      *
#   -------
#   seqinfo: 23 sequences from an unspecified genome; no seqlengths
# > fivekb_gr
# GRanges object with 29436 ranges and 2 metadata columns:
#           seqnames               ranges strand |   gene_id      symbol
#              <Rle>            <IRanges>  <Rle> | <integer> <character>
#       [1]     chr1     [  6874,  14319]      * | 100287102     DDX11L1
#       [2]     chr1     [ 14320,  33021]      * |    653635      WASH7P
#       [3]     chr1     [ 33022,  41080]      * |    641702     FAM138F
#       [4]     chr1     [ 64091,  74090]      * |     79501       OR4F5
#       [5]     chr1     [362659, 372658]      * |    729759      OR4F29
#       ...      ...                  ...    ... .       ...         ...
#   [29432]     chrY [27763264, 27773263]      * |      9085        CDY1
#   [29433]     chrY [27869637, 27879636]      * |    114760       TTTY3
#   [29434]     chrY [59095457, 59105456]      * |     10251       SPRY3
#   [29435]     chrY [59208949, 59218948]      * |      6845       VAMP7
#   [29436]     chrY [59325252, 59335251]      * |      3581        IL9R
#   -------
#   seqinfo: 25 sequences (1 circular) from hg19 genome
# > fivekb_gr[fivekb_gr$gene_id == 3581]
# GRanges object with 2 ranges and 2 metadata columns:
#       seqnames                 ranges strand |   gene_id      symbol
#          <Rle>              <IRanges>  <Rle> | <integer> <character>
#   [1]     chrX [155222246, 155232245]      * |      3581        IL9R
#   [2]     chrY [ 59325252,  59335251]      * |      3581        IL9R
#   -------
#   seqinfo: 25 sequences (1 circular) from hg19 genome
