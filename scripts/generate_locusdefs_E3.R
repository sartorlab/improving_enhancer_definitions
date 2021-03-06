library(GenomicRanges)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Have this handy
data('locusdef.hg19.5kb_outside', package='chipenrich.data')
fivekb_outside_gr = locusdef.hg19.5kb_outside@granges

data('locusdef.hg19.5kb', package='chipenrich.data')
fivekb_gr = locusdef.hg19.5kb@granges

chromhmm_pool = c(NA, 'chromhmm')
dnase_pool = c(NA, 'dnase')
thurman_pool = c(NA, 'thurman')
fantom_pool = c(NA, 'fantom')

expand_pool = c(0, 1000)

P2P_pool = c(NA, 'P2P')
E3_pool = c(NA, 'E3')

combinations_E3 = expand.grid(chromhmm_pool, dnase_pool, thurman_pool, fantom_pool, expand_pool, P2P_pool, E3_pool, thurman_pool, fantom_pool, stringsAsFactors=F)

# Get rid of 4 NAs in 1-4 and 6-9
combinations_E3 = combinations_E3[!apply(combinations_E3, 1, function(row){all(is.na(row[1:4]))}),]
combinations_E3 = combinations_E3[!apply(combinations_E3, 1, function(row){all(is.na(row[6:9]))}),]

# Get rid of link methods without E3

combinations_E3 <- combinations_E3[apply(combinations_E3[,6:9],1,function(x){any(x %in% "E3")}),]

# get rid of link methods containing "thurman" or "fantom"
# combinations_E3 = combinations_E3[apply(combinations_E3[,6:9],1,function(x){!(any(x %in% c("thurman","fantom"))) & any(x %in% "E3")}),]

# Make codes
enhancer_codes = apply(combinations_E3[,1:4], 1, function(row){paste(row[!is.na(row)], collapse='_')})
extension_codes = combinations_E3[,5]

enhancer_base = unique(paste(enhancer_codes, extension_codes, sep='.'))

ldef_qc_file = '/home/qinting/latte/share/improving_enhancer_definitions/locusdefs/E2_E3/QC_ldefs.txt'
if(!file.exists(ldef_qc_file)) {
    header = data.frame(
        'ldef_type' = 'ldef_type',
        'total_unique_regions' = 'total_unique_regions',
        'genome_coverage' = 'genome_coverage',
        'genes_per_enh' = 'genes_per_enh',
        'enhs_per_gene' = 'enhs_per_gene',
		stringsAsFactors = F)
    write.table(header, file = ldef_qc_file, sep = '\t', quote = F, col.names = F, row.names = F, append = T)
}

# Go by enhancer base to minimize the number of file reads
for(base in enhancer_base) {
	####
	# Read files
	message('Reading P2P...')
	base_p2p = read.table(file = sprintf('%s/%s.P2P.pairs', '/home/qinting/latte/share/improving_enhancer_definitions/pair_lists', base), sep='\t', header = T, as.is = T)
	message('Reading E...')
	base_e = read.table(file = sprintf('%s/%s.E3.pairs', '/home/qinting/latte/share/improving_enhancer_definitions/pair_lists', base), sep='\t', header = T, as.is = T)
	message('Reading thurman...')
	base_thurman = read.table(file = sprintf('%s/%s.thurman.pairs', '/home/qinting/latte/share/improving_enhancer_definitions/pair_lists', base), sep='\t', header = T, as.is = T)
	message('Reading fantom...')
	base_fantom = read.table(file = sprintf('%s/%s.fantom.pairs', '/home/qinting/latte/share/improving_enhancer_definitions/pair_lists', base), sep='\t', header = T, as.is = T)
	# message('Reading enhancer base...')
	enh_df = read.table(file = sprintf('%s/%s.enhancers.gz', '/home/qinting/latte/share/improving_enhancer_definitions/enhancer_lists', base), sep='\t', header = F, as.is = T)
	colnames(enh_df) = c('chr','start','end')

	####
	# Make GRanges of all the bases
	base_p2p_gr = GenomicRanges::makeGRangesFromDataFrame(df = base_p2p, keep.extra.columns = TRUE)
	base_e_gr = GenomicRanges::makeGRangesFromDataFrame(df = base_e, keep.extra.columns = TRUE)
	base_thurman_gr = GenomicRanges::makeGRangesFromDataFrame(df = base_thurman, keep.extra.columns = TRUE)
	base_fantom_gr = GenomicRanges::makeGRangesFromDataFrame(df = base_fantom, keep.extra.columns = TRUE)
	base_enh_gr = GenomicRanges::makeGRangesFromDataFrame(df = enh_df)

	P2P_pool = c(NA, 'P2P')
	E3_pool = c(NA, 'E3')
	thurman_pool = c(NA, 'thurman')
	fantom_pool = c(NA, 'fantom')
	nearest_pool = c(NA, 'nearest')
  nearest_All_pool = c(NA, 'nearest_All')
	interaction_combinations = expand.grid(P2P_pool, E3_pool, thurman_pool, fantom_pool, nearest_pool, nearest_All_pool, stringsAsFactors=F)
  # get rid of interactions without E3
  interaction_combinations = interaction_combinations[apply(interaction_combinations, 1, function(row){!(all(is.na(row))) & any(row %in% "E3")}),]
  # get rid of interactions with both nearest and nearest_All
  interaction_combinations = interaction_combinations[!(interaction_combinations[,5] %in% "nearest" & interaction_combinations[,6] %in% "nearest_All"),]

	####
	# Establish needfuls for assigning to nearest TSS
  dist_within_tss = 500000

	####
	# Go through each interaction combination and create the locus definition
	for(i in 1:nrow(interaction_combinations)) {
		E_code = interaction_combinations[i,2]
		P2P_code = interaction_combinations[i,1]
		thurman_code = interaction_combinations[i,3]
		fantom_code = interaction_combinations[i,4]
		nearest_code = interaction_combinations[i,5]
    nearest_All_code = interaction_combinations[i,6]

		ldef_gr = GenomicRanges::GRanges()

		interaction_row = interaction_combinations[i,]
		interaction_code = paste(interaction_row[!is.na(interaction_row)], collapse='_')

		# ldef_file = paste(base, interaction_code, 'ldef', 'gz', sep='.')
    ldef_file = paste('/home/qinting/latte/share/improving_enhancer_definitions/locusdefs/E2_E3/',paste(base, interaction_code, 'ldef', 'gz', sep='.'),sep='')

		qc_prefix = gsub('.ldef.gz', '', ldef_file)

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
    if (!is.na(nearest_code)) {
			message('Adding nearest')
			if(length(ldef_gr) != 0) {
				message('Adding nearest for missing enhancers only')
				# Then
				overlaps = findOverlaps(ldef_gr, base_enh_gr)
				missing_enh_idx = setdiff(seq_along(base_enh_gr), subjectHits(overlaps))

				if(length(findOverlaps(ldef_gr, base_enh_gr[missing_enh_idx])) > 0) {
					stop((sprintf('Some "missing" enhancers still intersect %s', base)))
				}

				base_enh_gr = base_enh_gr[missing_enh_idx]
			} else {
				message('Adding nearest for all enhancers')
			}

			# Then assign ALL of the enhancers to the nearest TSS (5kb)
			nearests = GenomicRanges::nearest(base_enh_gr, fivekb_gr, ignore.strand = TRUE, select='all')
			distances = GenomicRanges::distanceToNearest(base_enh_gr, fivekb_gr, ignore.strand = TRUE, select='all')

			keep = which(GenomicRanges::mcols(distances)$distance < dist_within_tss)

			if(length(queryHits(nearests)) != length(unique(queryHits(nearests)))) {
				warning(sprintf('There are ties for nearest 5kb with the enhancer base: %s', base))
			}

			tmp_gr = base_enh_gr[queryHits(nearests)[keep]]
			GenomicRanges::mcols(tmp_gr)$gene_id = fivekb_gr[subjectHits(nearests)[keep]]$gene_id

			ldef_gr = c(ldef_gr, tmp_gr)
		}
    if (!is.na(nearest_All_code)) {
			message('Adding all peaks outsied enhancer ldefs but 5kb_outside of a gene to that gene')
      ### find overlapped ranges between 5kb_outside ldef and enhancer ldef
      ovlp_outside5kb_vs_ldef = data.frame(findOverlaps(fivekb_outside_gr,ldef_gr))
      ovlpIndex.list = split(ovlp_outside5kb_vs_ldef$subjectHits, ovlp_outside5kb_vs_ldef$queryHits)
      ### for each overlapped 5kb_outside ldef, substract the overlapped enhancer ldef and keep the remaining region linked the nearest tss gene (using 5kb ldef)
      subsetLdef = mapply(function(x,y){tmp=GenomicRanges::setdiff(fivekb_outside_gr[x], ldef_gr[y]); mcols(tmp)=data.frame(gene_id= mcols(fivekb_outside_gr[x])$gene_id); tmp}, as.integer(names(ovlpIndex.list)), ovlpIndex.list)
      ### concatenate all subset ldef
      subsetLdef = do.call(c, subsetLdef)
      ### add all other 5kb_outside ldefs without any enhancer ldef
      tmp = fivekb_outside_gr[-unique(ovlp_outside5kb_vs_ldef[,1])]
      mcols(tmp)$symbol <- NULL
      subsetLdef = c(subsetLdef,tmp)
      ### add subsetLdef back to enhancer ldef
      ldef_gr = c(ldef_gr, subsetLdef)
		}

		message('Reducing...')
    seqlevels(ldef_gr) = sort(seqlevels(ldef_gr))
    ldef_gr = sort(ldef_gr)
		ldef_grl = IRanges::splitAsList(ldef_gr, ldef_gr$gene_id)
		ldef_grl = GenomicRanges::reduce(ldef_grl)
		ldef_grl_genes = S4Vectors::Rle(names(ldef_grl), S4Vectors::elementNROWS(ldef_grl))

		ldef_gr = unlist(ldef_grl, use.names = FALSE)
		ldef_gr$gene_id = ldef_grl_genes

		ldef_gr = sort(ldef_gr)

		ldef_df = BiocGenerics::as.data.frame(ldef_gr)
		colnames(ldef_df) = c('chr','start','end','width','strand','gene_id')
		ldef_df = subset(ldef_df, width > 5)
		ldef_df = ldef_df[, c('chr','start','end','gene_id')]

    # Get Entrez ID to gene symbol mappings for custom locus definitions
    # currently all enhancer ldefs are hg19
    # egSYMBOL = org.Hs.eg.db::org.Hs.egSYMBOL
    # ### Build Entrez ID to gene symbol mapping
    # mapped_genes = AnnotationDbi::mappedkeys(egSYMBOL)
    # eg2symbol = as.data.frame(egSYMBOL[mapped_genes])
    # eg2symbol$gene_id = as.integer(eg2symbol$gene_id)
    #
    # # map the gene Entrez ID in ldef to symbol
    # tmp = data.frame(mcols(ldef_gr))
    # tmp$symbol = eg2symbol[match(tmp$gene_id, eg2symbol[,1]),2]
    # mcols(ldef_gr) = tmp

		message(sprintf('Writing %s...', ldef_file))
		write.table(ldef_df, file = gzfile(ldef_file), sep = '\t', quote = F, row.names = F, col.names = T)

		#################################
	    # QC
		total_unique_regions = length(unique(ldef_gr))
		genome_coverage = sum(width(unique(ldef_gr)))

		enh_hash = paste(ldef_df$chr, ldef_df$start, ldef_df$end, sep='')

		genes_per_enh = as.integer(table(enh_hash))
		enhs_per_gene = as.integer(table(as.integer(ldef_gr$gene_id)))

		pdf(file = sprintf('%s_QC_gene_enh_per.pdf', qc_prefix), width = 12, height = 6)
			par(mfrow = c(1,2))
			hist(genes_per_enh, main = 'Genes per enhancer', xlab='')
			hist(enhs_per_gene, main = 'Enhancers per gene', xlab='')
		dev.off()

	    qc_text = data.frame(
			'ldef_type' = qc_prefix,
			'total_unique_regions' = total_unique_regions,
			'genome_coverage' = genome_coverage,
			'genes_per_enh' = mean(genes_per_enh),
			'enhs_per_gene' = mean(enhs_per_gene),
			stringsAsFactors = F)
	    write.table(qc_text, file = ldef_qc_file, sep = '\t', quote = F, col.names = F, row.names = F, append = T)
		#################################
	}
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
# > base_enh_gr
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
