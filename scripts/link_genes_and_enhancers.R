suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--gene_loci"), action = "store", type = "character", help = "[Required] Path to a file listing the genomic intervals that should be assigned to each gene.  Format is one line per locus, tab-separated, chromosome, start, end, gene_id (no header)"),
    make_option(c("--enhancers"), action = "store", type = "character", help = "[Required] Path to a file listing the enhancers.  Format is one line per locus, tab-separated, chromosome, start, end, enhancer_id (no header)"),
    make_option(c("--interactions"), action = "store", type = "character", help = "[Required] Path to a file of interactions. Format is one interaction per line, tab separated, chromosome_1, start_1, end_1, chromosome_2, start_2, end_2"),
    make_option(c("--method"), action = "store", type = "character", help = "[Required] The method used to assign genes to enhancers; must be either 'point_to_point' or 'encompassing'"),
    make_option(c("--out"), action = "store", type = "character", help = "[Required] File to which to write the resulting pairs"),
    make_option(c("--max_genes_per_region"), action = "store", type = "numeric", help = "[Optional] Ignored unless --method = 'encompassing'.  If this is the case, then gene-enhancer pairs are not assigned within regions (e.g. TADs) that encompass more than max_genes_per_region (default = no threshold)")
  )

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T,
                              description = "\nDescription: Accepts a list of gene loci, enhancers, and an interaction dataset (e.g., ChIA-PET or Hi-C).  Prints a list of gene-enhancer pairs. Generates this list using one of two methods: either by finding the interactions that link one gene to one enhancer (each end of the interaction overlaps with one of the genes/enhancers; 'point_to_point') or by declaring all gene/enhancer pairs encompassed WITHIN an interaction region to be interacting ('encompassing'; useful when e.g. enhancers/genes are expected to be within the loops rather than at the interaction edges, such as would be expected within a TAD).")
opts <- parse_args(option_parser)


#####  VERIFY THE INPUT  #####
if (length(opts$method) == 0) {
  stop("--method is a required argument")
}
if (length(opts$gene_loci) == 0) {
  stop("--gene_loc is a required argument")
}
if (length(opts$enhancers) == 0) {
  stop("--enhancers is a required argument")
}
if (length(opts$interactions) == 0) {
  stop("--interactions is a required argument")
}
if (length(opts$out) == 0) {
  stop("--out is a required argument")
}
if (opts$method != "point_to_point" & opts$method != "encompassing") {
  stop("'--method' argument must be one of 'point_to_point' or 'encompassing'")
}
if (length(opts$max_genes_per_region) > 0 & opts$method == "point_to_point") {
  warning("Ignoring --max_genes_per_region argument (doesn't apply when '--method point_to_point' is used)")
}

#####  LOAD THE REST OF THE LIBRARIES  #####
suppressPackageStartupMessages(library("GenomicRanges"))


#####  MAIN  #####
message(sprintf('Building %s...', opts$out))

###
message('Reading genes...')
genes = read.table(opts$gene_loci, sep='\t', header = F, as.is = T)
colnames(genes) = c('chr','start','end','gene_id','symbol')

genes_gr = GenomicRanges::makeGRangesFromDataFrame(df = genes, keep.extra.columns = T)

###
message(sprintf('Reading %s...', opts$enhancers))
enhancers = read.table(opts$enhancers, sep = '\t', header = F, as.is = T)
colnames(enhancers) = c('chr','start','end')

enhancers_gr = GenomicRanges::makeGRangesFromDataFrame(df = enhancers)

###
message(sprintf('Reading %s...', opts$interactions))
interactions = read.table(opts$interactions, sep='\t', header = F, as.is = T)
colnames(interactions) = c('chrom1','start1','end1','chrom2','start2','end2')
interactions = interactions[interactions$chrom1 == interactions$chrom2,]

if (opts$method == "point_to_point") {

  message('Building P2P pairs...')
  ###
  interactions1_gr = GenomicRanges::makeGRangesFromDataFrame(df = interactions, seqnames.field='chrom1', start.field='start1', end.field='end1')
  interactions2_gr = GenomicRanges::makeGRangesFromDataFrame(df = interactions, seqnames.field='chrom2', start.field='start2', end.field='end2')

  ###
  genes_interactions1 = as.data.frame(GenomicRanges::findOverlaps(genes_gr, interactions1_gr))
  genes_interactions2 = as.data.frame(GenomicRanges::findOverlaps(genes_gr, interactions2_gr))

  enhancers_interactions1 = as.data.frame(GenomicRanges::findOverlaps(enhancers_gr, interactions1_gr))
  enhancers_interactions2 = as.data.frame(GenomicRanges::findOverlaps(enhancers_gr, interactions2_gr))

  ###
  # The indices for the subjectHits are the same. Require enhancer / gene or gene / enhancer
  # overlap in interactions1_gr and interactions2_gr, respectively
  idx1 = merge(genes_interactions1, enhancers_interactions2, by = 'subjectHits', suffixes= c('.genes','.enh'))
  idx2 = merge(genes_interactions2, enhancers_interactions1, by = 'subjectHits', suffixes= c('.genes','.enh'))

  ###
  loci1_gr = enhancers_gr[idx1$queryHits.enh]
  loci1_gr$gene_id = genes_gr$gene_id[idx1$queryHits.genes]
  loci1_gr$symbol = genes_gr$symbol[idx1$queryHits.genes]

  loci2_gr = enhancers_gr[idx2$queryHits.enh]
  loci2_gr$gene_id = genes_gr$gene_id[idx2$queryHits.genes]
  loci2_gr$symbol = genes_gr$symbol[idx2$queryHits.genes]

  loci_gr = sort(c(loci1_gr, loci2_gr))
  loci_gr = unique(loci_gr)

  loci_df = as.data.frame(loci_gr)
  colnames(loci_df) = c('chr','start','end','width','strand','gene_id','symbol')

  loci_df = loci_df[,c('chr','start','end','gene_id','symbol')]

  write.table(x = loci_df, file = opts$out, append = F, quote = F, sep = "\t", row.names = F, col.names = T)

} else if (opts$method == "encompassing") {

  message('Building E pairs...')
  # Create the loops from start to start
  loops_gr = GenomicRanges::GRanges(
  	seqnames = interactions$chrom1,
  	ranges = IRanges::IRanges(start = interactions$start1, end = interactions$start2)
  )

  ###
  genes_in_loops = as.data.frame(GenomicRanges::findOverlaps(genes_gr, loops_gr, type = 'within'))
  enhancers_in_loops = as.data.frame(GenomicRanges::findOverlaps(enhancers_gr, loops_gr, type = 'within'))

  ###
  # Each subjectHit index represents a loop. Use table to count how many genes fall in each loop
  num_genes_in_loops = table(genes_in_loops$subjectHits)

  valid_loop_idx = names(num_genes_in_loops[num_genes_in_loops <= opts$max_genes_per_region])

  ###
  genes_in_loops = subset(genes_in_loops, subjectHits %in% valid_loop_idx)
  enhancers_in_loops = subset(enhancers_in_loops, subjectHits %in% valid_loop_idx)

  ###

  idx = merge(genes_in_loops, enhancers_in_loops, by = 'subjectHits', suffixes = c('.genes', '.enh'))

  ###
  loci_gr = enhancers_gr[idx$queryHits.enh]
  loci_gr$gene_id = genes_gr$gene_id[idx$queryHits.genes]
  loci_gr$symbol = genes_gr$symbol[idx$queryHits.genes]

  loci_gr = sort(loci_gr)
  loci_gr = unique(loci_gr)

  loci_df = as.data.frame(loci_gr)
  colnames(loci_df) = c('chr','start','end','width','strand','gene_id','symbol')

  loci_df = loci_df[,c('chr','start','end','gene_id','symbol')]

  write.table(x = loci_df, file = opts$out, append = F, quote = F, sep = "\t", row.names = F, col.names = T)
}
