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
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))


#####  MAIN  #####
qc_prefix = gsub('.pairs', '', opts$out)
p2p_qc_file = 'QC_P2P_interactions.txt'
if(!file.exists(p2p_qc_file)) {
    header = data.frame(
        'pair_type' = 'pair_type',
        'total_interactions' = 'total_interactions',
        'gene_enh_pairs' = 'gene_enh_pairs',
        'gene_gene_pairs' = 'gene_gene_pairs',
        'enh_enh_pairs' = 'enh_enh_pairs', stringsAsFactors = F)
    write.table(header, file = p2p_qc_file, sep = '\t', quote = F, col.names = F, row.names = F, append = T)
}

e_qc_file = 'QC_E_interactions.txt'
if(!file.exists(e_qc_file)) {
    header = data.frame(
        'pair_type' = 'pair_type',
        'total_loops' = 'total_loops',
        'valid_loops' = 'valid_loops',
        'loops_per_gene' = 'loops_per_gene',
        'genes_per_loop' = 'genes_per_loop',
        'loops_per_enh' = 'loops_per_enh',
        'enhs_per_loop' = 'enhs_per_loop',
        'enhs_per_gene' = 'enhs_per_gene',
        'genes_per_enh' = 'genes_per_enh',
        stringsAsFactors = F)
    write.table(header, file = e_qc_file, sep = '\t', quote = F, col.names = F, row.names = F, append = T)
}


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

  gene_gene_idx = merge(genes_interactions1, genes_interactions2, by = 'subjectHits', suffixes = c('.1','.2'))
  enh_enh_idx = merge(enhancers_interactions1, enhancers_interactions2, by = 'subjectHits', suffixes = c('.1','.2'))

  ###
  loci1_gr = enhancers_gr[idx1$queryHits.enh]
  loci1_gr$gene_id = genes_gr$gene_id[idx1$queryHits.genes]

  loci2_gr = enhancers_gr[idx2$queryHits.enh]
  loci2_gr$gene_id = genes_gr$gene_id[idx2$queryHits.genes]

  loci_gr = sort(c(loci1_gr, loci2_gr))

  # Split by gene_id before doing unique
  loci_grl = split(loci_gr, loci_gr$gene_id)
  loci_grl = unique(loci_grl)
  loci_gr_geneid = S4Vectors::Rle(names(loci_grl), S4Vectors::elementNROWS(loci_grl))

  loci_gr = unlist(loci_grl, use.names = FALSE)
  loci_gr$gene_id = loci_gr_geneid

  loci_gr = sort(loci_gr)

  loci_df = as.data.frame(loci_gr)
  colnames(loci_df) = c('chr','start','end','width','strand','gene_id')

  loci_df = loci_df[,c('chr','start','end','gene_id')]

  write.table(x = loci_df, file = opts$out, append = F, quote = F, sep = "\t", row.names = F, col.names = T)

    #################################
    # QC
    qc_text = data.frame(
        'pair_type' = qc_prefix,
        'total_interactions' = nrow(interactions),
        'gene_enh_pairs' = nrow(idx1) + nrow(idx2),
        'gene_gene_pairs' = nrow(gene_gene_idx),
        'enh_enh_pairs' = nrow(enh_enh_idx),
        stringsAsFactors = F)
    write.table(qc_text, file = p2p_qc_file, sep = '\t', quote = F, col.names = F, row.names = F, append = T)

    # Distances
    # Distance from gene to enhancers
    unique_gene_enh_idx = unique(rbind(idx1[,2:3], idx2[,2:3]))

    # Interaction distances for gene gene
    unique_gene_gene_idx = unique(gene_gene_idx[,1])

    # Interaction distances for enh enh
    unique_enh_enh_idx = unique(enh_enh_idx[,1])

    distance_list = list(
        'gene_enh' = distance(genes_gr[unique_gene_enh_idx$queryHits.genes], enhancers_gr[unique_gene_enh_idx$queryHits.enh]),
        'gene_gene' = distance(interactions1_gr[unique_gene_gene_idx], interactions2_gr[unique_gene_gene_idx]),
        'enh_enh' = distance(interactions1_gr[unique_enh_enh_idx], interactions2_gr[unique_enh_enh_idx]),
        'nearest' = mcols(distanceToNearest(enhancers_gr, genes_gr))$distance)
    combined = reshape2::melt(distance_list)
    colnames(combined) = c('distance', 'type')

    nonzeros = which(combined$distance != 0)
    combined$distance[nonzeros] = log10(combined$distance[nonzeros])

    plot_distance = ggplot(data = combined, aes(x = distance)) +
        geom_histogram(binwidth = 0.25) +
        facet_grid(. ~ type) +
        theme_bw() + ggtitle(sprintf('%s connection distances', qc_prefix)) + xlab('log10_distance')
    ggsave(plot_distance, filename = sprintf('QC_%s_connection_distances.pdf', qc_prefix), height = 5, width = 20)

    # Interaction widths
    interaction_widths = width(c(interactions1_gr, interactions2_gr))
    #################################

} else if (opts$method == "encompassing") {

  message('Building E pairs...')
  # Create the loops from midpoint to midpoint
  # Interactions can be quite wide
  loops_gr = GenomicRanges::GRanges(
  	seqnames = interactions$chrom1,
  	ranges = IRanges::IRanges(start = floor((interactions$start1 + interactions$end1)/2), end = floor((interactions$start2 + interactions$end2)/2))
  )

  ###
  genes_in_loops = as.data.frame(GenomicRanges::findOverlaps(genes_gr, loops_gr, type = 'within'))
  enhancers_in_loops = as.data.frame(GenomicRanges::findOverlaps(enhancers_gr, loops_gr, type = 'within'))

  ###
  # Each subjectHit index represents a loop. Use table to count how many genes fall in each loop
  num_genes_in_loops = table(genes_in_loops$subjectHits)
  num_enh_in_loops = table(enhancers_in_loops$subjectHits)

  ###
  valid_loop_idx = as.integer(names(num_genes_in_loops[num_genes_in_loops <= opts$max_genes_per_region]))
  valid_genes_in_loops = subset(genes_in_loops, subjectHits %in% valid_loop_idx)
  valid_enhancers_in_loops = subset(enhancers_in_loops, subjectHits %in% valid_loop_idx)

  ###
  idx = merge(valid_genes_in_loops, valid_enhancers_in_loops, by = 'subjectHits', suffixes = c('.genes', '.enh'))

  ###
  loci_gr = enhancers_gr[idx$queryHits.enh]
  loci_gr$gene_id = genes_gr$gene_id[idx$queryHits.genes]

  # Split by gene_id before doing unique
  loci_grl = split(loci_gr, loci_gr$gene_id)
  loci_grl = unique(loci_grl)
  loci_gr_geneid = S4Vectors::Rle(names(loci_grl), S4Vectors::elementNROWS(loci_grl))

  loci_gr = unlist(loci_grl, use.names = FALSE)
  loci_gr$gene_id = loci_gr_geneid

  loci_gr = sort(loci_gr)

  loci_df = as.data.frame(loci_gr)
  colnames(loci_df) = c('chr','start','end','width','strand','gene_id')

  loci_df = loci_df[,c('chr','start','end','gene_id')]

  write.table(x = loci_df, file = opts$out, append = F, quote = F, sep = "\t", row.names = F, col.names = T)

    #################################
    # QC
    num_total_loops = length(loops_gr)
    num_valid_loops = length(valid_loop_idx)

    unique_loop_gene = unique(idx[,c('subjectHits','queryHits.genes')])
    unique_loop_enh = unique(idx[,c('subjectHits','queryHits.enh')])
    unique_gene_enh = unique(idx[,c('queryHits.genes','queryHits.enh')])

    # Pers
    per_list = list(
        'loops_per_gene' = table(unique_loop_gene$queryHits.genes),
        'genes_per_loop' = table(unique_loop_gene$subjectHits),
        'loops_per_enh' = table(unique_loop_enh$queryHits.enh),
        'enhs_per_loop' = table(unique_loop_enh$subjectHits),
        'enhs_per_gene' = table(unique_gene_enh$queryHits.genes),
        'genes_per_enh' = table(unique_gene_enh$queryHits.enh))
    per = reshape2::melt(lapply(per_list, as.integer))
    colnames(per) = c('per', 'type')

    pdf(file = sprintf('QC_%s.pdf', qc_prefix), width = 12, height = 18)
    par(mfrow = c(3,2))
    for(t in unique(per$type)) {
        hist(subset(per, type == t)$per, xlab = t, main = sprintf('%s %s', qc_prefix, t))
    }
    dev.off()

    if(!file.exists(sprintf('QC_loop_widths.pdf', qc_prefix))) {
        # Loop widths
        loop_widths_list = list(
            'all' = log10(width(loops_gr)),
            'valid' = log10(width(loops_gr[valid_loop_idx])))
        loop_widths = reshape2::melt(loop_widths_list)
        colnames(loop_widths) = c('log10_width', 'type')

        plot_widths = ggplot(data = loop_widths, aes(x = log10_width)) +
            geom_histogram(binwidth = 0.25) +
            facet_grid(. ~ type) +
            theme_bw() + ggtitle(sprintf('Loop widths', qc_prefix))
        ggsave(plot_widths, filename = sprintf('QC_loop_widths.pdf', qc_prefix), height = 6, width = 8)
    }

    # Text
    qc_text = data.frame(
        'pair_type' = qc_prefix,
        'total_loops' = num_total_loops,
        'valid_loops' = num_valid_loops,
        t(sapply(per_list, function(p){mean(as.integer(p))})),
        stringsAsFactors = F)
    write.table(qc_text, file = e_qc_file, sep = '\t', quote = F, col.names = F, row.names = F, append = T)
    #################################
}
