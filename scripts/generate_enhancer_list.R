suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(c("--genes"), action = "store", type = "character", help = "[Required] Path to a file listing the genomic intervals that should be assigned to each gene.  Format is one line per locus, tab-separated, chromosome, start, end, gene_id (no header)"),
  make_option(c("--chromhmm"), action = "store", type = "character", help = "[Optional] File of concatenated chromHMM tracks"),
  make_option(c("--dnase"), action = "store", type = "character", help = "[Optional] Path to a file listing the enhancers.  Format is one line per locus, tab-separated, chromosome, start, end, enhancer_id (no header)"),
  make_option(c("--thurman"), action = "store", type = "character", help = "[Optional] Path to a file listing enhancers (and interactions) from Thurman et al. paper."),
  make_option(c("--fantom"), action = "store", type = "character", help = "[Optional] Path to a file listing enhancers from the FANTOM5 consortium."),
  make_option(c("--extension"), action = "store", type = "numeric", help = "[Required] Number of base pairs to extend enhancers to"),
  make_option(c("--tissue_threshold_dnase"), action = "store", type = "numeric", default = 1, help = "[Optional] An integer.  If the DNase list is used to create the enhancer list, then only DNase sites supported by this number of experiments are included.  There are ~125 experiments in the ENCODE DNase list (Default: 1).")
)

option_parser = OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T,
                              description = "\nDescription: Accepts a list of gene loci, enhancers, and an interaction dataset (e.g., ChIA-PET or Hi-C).  Prints a list of gene-enhancer pairs. Generates this list using one of two methods: either by finding the interactions that link one gene to one enhancer (each end of the interaction overlaps with one of the genes/enhancers; 'point_to_point') or by declaring all gene/enhancer pairs encompassed WITHIN an interaction region to be interacting ('encompassing'; useful when e.g. enhancers/genes are expected to be within the loops rather than at the interaction edges, such as would be expected within a TAD).")
opts = parse_args(option_parser)

#####  VERIFY THE INPUT  #####
if (length(opts$genes) == 0) {
  stop("--genes is a required argument")
}
if( (length(opts$dnase) + length(opts$chromhmm) + length(opts$fantom) + length(opts$thurman)) != 4) {
  stop("All of --dnase, --chromhmm, --fantom, or --thurman must be given")
}
if (length(opts$extension) == 0) {
  stop("--extension is a required argument")
}

#####  LOAD THE REST OF THE LIBRARIES  #####
suppressPackageStartupMessages(library("GenomeInfoDb"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))

hg19_seqinfo = GenomeInfoDb::Seqinfo(genome = 'hg19')

#####  FUNCTION TO EXPAND AND RESECT GRANGES  #####
expand_and_resect <- function(granges, min_size) {
  # extends each grange to the min_size (if not there already); if this would cause two neighboring ranges to overlap, it extends them only until they bookend

  # make sure everything overlapping already has been merged
  # this WILL NOT merge ranges that bump up against each other
  granges <- reduce(granges, min.gapwidth = 0)

  # resize the ranges that are smaller than min_size; keep the ranges that are bigger than min_size as is.
  was_resized <- rep(F, length(granges))
  was_resized[(end(granges) - start(granges)) < min_size] <- T

  granges[(end(granges) - start(granges)) < min_size] <- resize(granges[(end(granges) - start(granges)) < min_size], width = min_size, fix = "center")

  # find the cases where, after resizing, two ranges overlap
  overlaps <- findOverlaps(granges, granges)
  overlaps <- overlaps[(queryHits(overlaps)+1) == subjectHits(overlaps),]


  resect_first <- queryHits(overlaps)[was_resized[queryHits(overlaps)] & !was_resized[subjectHits(overlaps)]]
  resect_second <- queryHits(overlaps)[!was_resized[queryHits(overlaps)] & was_resized[subjectHits(overlaps)]]
  resect_both <- queryHits(overlaps)[was_resized[queryHits(overlaps)] & was_resized[subjectHits(overlaps)]]

  # resect
  end(granges)[resect_first] <- start(granges)[resect_first + 1]
  start(granges)[resect_second + 1] <- end(granges)[resect_second]
  overlap_amount <- end(granges)[resect_both] - start(granges)[resect_both+1]
  resection_amount <- ceiling(overlap_amount / 2)
  end(granges)[resect_both] <- end(granges)[resect_both] - resection_amount
  start(granges)[resect_both+1] <- start(granges)[resect_both+1] + resection_amount


  # NOTE: need to know the chromosome sizes so we can truncate at the chromosome starts/ends.  For now, just truncate at the chromosome starts...
  start(granges)[start(granges) <= 0] <- 1

  return(granges)
}

expand_and_resect2 = function(granges, min_size) {
    granges = reduce(granges, min.gapwidth = 0)

    resize = rep(F, length(granges))
    resize[width(granges) <= min_size] = T

    granges[resize] = resize(granges[resize], width = min_size, fix = 'center')

    granges = reduce(granges, min.gapwidth = 0)

    return(granges)
}

#####
#####  MAIN  #####
#####

#####  READ GENE LOCI  #####
genes = read.table(opts$genes, header = F, as.is = T, sep = '\t')
colnames(genes) = c("chromosome", "start", "end", "gene_id", "symbol")
genes_gr = makeGRangesFromDataFrame(df = genes, seqnames.field = "chromosome", start.field = "start", end.field = "end")

#####  DETERMINE COMBINATIONS  #####
chromhmm_pool = c(NA, 'chromhmm')
dnase_pool = c(NA, 'dnase')
thurman_pool = c(NA, 'thurman')
fantom_pool = c(NA, 'fantom')

combinations = expand.grid(chromhmm_pool, dnase_pool, thurman_pool, fantom_pool, stringsAsFactors=F)
# Get rid of all NA
combinations = combinations[!apply(combinations, 1, function(row){all(is.na(row))}),]

#####  LOAD ENHANCER COMPONENTS  #####
# load in the chromhmm data for all the tissues
message('Reading chromhmm...')
chromhmm = read.table(opts$chromhmm, header = F, as.is = T, sep = '\t')
colnames(chromhmm) = c("chromosome", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb")
chromhmm = chromhmm[grep("Enhancer", chromhmm$name),]
chromhmm = chromhmm[,c("chromosome", "start", "end")]

# load in the dnase data
message('Reading dnase...')
dnase = read.table(opts$dnase, header = F, as.is = T, sep = '\t')
colnames(dnase) = c("chromosome", "start", "end", "name", "score", "floatScore", "sourceCount", "sourceIds")
dnase = dnase[dnase$sourceCount >= opts$tissue_threshold_dnase, ]
dnase = dnase[, c(1:3)]

# load in the fantom data
message('Reading fantom...')
fantom = read.table(opts$fantom, header = F, as.is = T, sep = '\t')
fantom = fantom[, c(1:3)]
colnames(fantom) = c('chromosome', 'start', 'end')

# load in the thurman data
message('Reading thurman...')
thurman = read.table(opts$thurman, header = F, as.is = T, sep = '\t')
colnames(thurman) = c('promoter.chromosome', 'promoter.start', 'promoter.end', 'promoter.symbol', 'distal.chromosome', 'distal.start', 'distal.end', 'cor')
thurman = thurman[, c('distal.chromosome', 'distal.start', 'distal.end')]
colnames(thurman) = c('chromosome', 'start', 'end')

if (opts$extension == 1) {
    extension = 0
} else {
    extension = opts$extension
}

#################################
# QC
qc_summary = data.frame()

widths_list = lapply(list('chromhmm' = chromhmm, 'dnase' = dnase, 'fantom' = fantom , 'thurman' = thurman), function(e){log2(e$end - e$start)})
widths_df = melt(widths_list)
colnames(widths_df) = c('log2_width', 'enhancer_type')

gg_density = ggplot(widths_df, aes(log2_width, ..density.., fill = enhancer_type)) +
    geom_histogram(binwidth = 0.5) +
    geom_vline(xintercept = log2(1000)) +
    ggtitle('Enhancer widths by type (bins = 0.5)')

gg_counts = ggplot(widths_df, aes(log2_width, fill = enhancer_type)) +
    geom_histogram(binwidth = 0.5) +
    geom_vline(xintercept = log2(1000)) +
    ggtitle('Enhancer widths by type (bins = 0.5)')

ggsave(gg_density, filename = 'QC_enhancer_densities.pdf', width = 8, height = 6)
ggsave(gg_counts, filename = 'QC_enhancer_counts.pdf', width = 8, height = 6)
#################################

#####  BUILD THE COMBINATIONS  #####
for(i in 1:nrow(combinations)) {
    chromhmm_code = combinations[i,1]
    dnase_code = combinations[i,2]
    thurman_code = combinations[i,3]
    fantom_code = combinations[i,4]

    out_code = c()
    if(!is.na(chromhmm_code)) {
        out_code = c(out_code, chromhmm_code)
    }
    if(!is.na(dnase_code)) {
        out_code = c(out_code, dnase_code)
    }
    if(!is.na(thurman_code)) {
        out_code = c(out_code, thurman_code)
    }
    if(!is.na(fantom_code)) {
        out_code = c(out_code, fantom_code)
    }

    qc_prefix = paste(paste(out_code, collapse = '_'), extension, sep = '.')
    out_file = paste(paste(out_code, collapse = '_'), extension, 'enhancers', 'gz', sep = '.')

    message(sprintf('On %s', out_file))

    if(file.exists(out_file)) {
        message('File exists, skipping...')
        next
    }

    #################################
    # Build enhancer set
    # Add the correct enhancer pieces to a data.frame
    enhancers_df = data.frame()
    if (!is.na(chromhmm_code)) {
        message('Adding chromhmm')
        enhancers_df = rbind(enhancers_df, chromhmm)
    }
    if (!is.na(dnase_code)) {
        message('Adding dnase')
        enhancers_df = rbind(enhancers_df, dnase[,c("chromosome", "start", "end")])
    }
    if (!is.na(thurman_code)) {
        message('Adding thurman')
        enhancers_df = rbind(enhancers_df, thurman)
    }
    if (!is.na(fantom_code)) {
        message('Adding fantom')
        enhancers_df = rbind(enhancers_df, fantom)
    }

    # Create a GRanges object
    enhancers_gr = makeGRangesFromDataFrame(enhancers_df, ignore.strand = TRUE)

    # Expand and resect the enhancers
    # Note, this has a reduce() step
    enhancers_merged = expand_and_resect2(enhancers_gr, opts$extension)

    # Remove any part of genes_gr from the enhancers_merged
    # There is another reduce() step here
    enhancers_final = GenomicRanges::setdiff(enhancers_merged, genes_gr)

    # Convert back to data.frame
    enhancers_df = as.data.frame(enhancers_final)
    enhancers_df = enhancers_df[order(enhancers_df$seqnames, enhancers_df$start), ]

    message(sprintf('Writing %s', out_file))
    write.table(enhancers_df[,c("seqnames", "start", "end")], file = gzfile(out_file), quote = F, sep = "\t", row.names = F, col.names = F)
    #################################

    #################################
    # QC
    # Information about numbers of enhancer regions at different stages, and
    # the number of overlaps at different stages
    qc_pre = length(enhancers_gr)
    qc_post = length(enhancers_merged)
    qc_final = length(enhancers_final)
    pre_overlaps_5kb_overall = length(findOverlaps(enhancers_merged, genes_gr))
    pre_overlaps_5kb_within = length(findOverlaps(enhancers_merged, genes_gr, type = 'within'))
    pre_overlaps_5kb_trim = pre_overlaps_5kb_overall - pre_overlaps_5kb_within
    post_overlaps_5kb_overall = length(findOverlaps(enhancers_final, genes_gr))
    post_overlaps_5kb_within = length(findOverlaps(enhancers_final, genes_gr, type = 'within'))
    post_overlaps_5kb_trim = post_overlaps_5kb_overall - post_overlaps_5kb_within

    # Information about the widths of regions and distances to nearest regions
    pre_widths_df = data.frame(type = rep.int('pre', length(enhancers_gr)), log2_width = log2(width(enhancers_gr)), stringsAsFactors = F)
    pre_nearest_distances = as.data.frame(distanceToNearest(enhancers_gr))

    post_widths_df = data.frame(type = rep.int('post', length(enhancers_merged)), log2_width = log2(width(enhancers_merged)), stringsAsFactors = F)
    post_nearest_distances = as.data.frame(distanceToNearest(enhancers_merged))

    # Information about the final regions used
    final_widths_df = data.frame(log2_width = log2(width(enhancers_final)))
    final_nearest_distances_df = as.data.frame(distanceToNearest(enhancers_final))

    # Visualizatinos combining the information pre and post resecttion/expansion
    pre_widths_label = 'Pre-expand/resect widths'
    post_widths_label = 'Post-expand/resect widths'
    fill_manual = c('gray', NA)
    names(fill_manual) = c(pre_widths_label, post_widths_label)
    compare_widths_gg = ggplot(data = pre_widths_df, aes(x = log2_width)) +
      geom_histogram(binwidth = 0.5, aes(fill = pre_widths_label)) +
      geom_histogram(data = post_widths_df, binwidth = 0.5, aes(fill = post_widths_label, color = 'red')) +
      scale_fill_manual(values = fill_manual) +
      geom_vline(xintercept = log2(1000)) +
      ggtitle(sprintf('%s widths (bins = 0.5)', qc_prefix)) +
      xlab('log2_width') +
      guides(color = FALSE) +
      theme_bw() + theme(legend.title=element_blank(), legend.position="bottom", legend.key = element_rect(color = c('red','white')))
    ggsave(compare_widths_gg, filename = sprintf('QC_%s_compare_widths.pdf', qc_prefix), width = 8, height = 6)

    pre_distances_label = 'Pre-expand/resect distances'
    post_distances_label = 'Post-expand/resect distances'
    fill_manual = c('gray', NA)
    names(fill_manual) = c(pre_distances_label, post_distances_label)
    compare_distances_gg = ggplot(data = pre_nearest_distances, aes(x = distance)) +
      geom_histogram(binwidth = 100, aes(fill = pre_distances_label)) +
      geom_histogram(data = post_nearest_distances, binwidth = 100, aes(fill = post_distances_label, col = 'red')) +
      scale_fill_manual(values = fill_manual) +
      scale_x_continuous(limits = c(-100, 1500)) +
      geom_vline(xintercept = 500) +
      ggtitle(sprintf('%s nearest distances (bins = 100)', qc_prefix)) +
      xlab('Distance to nearest') +
      guides(color = FALSE) +
      theme_bw() + theme(legend.title=element_blank(), legend.position="bottom", legend.key = element_rect(color = c('red','white')))
    ggsave(compare_distances_gg, filename = sprintf('QC_%s_compare_distances.pdf', qc_prefix), width = 8, height = 6)

    # Visualizations for the final regions
    final_widths_gg = ggplot(final_widths_df, aes(log2_width)) +
        geom_histogram(binwidth = 0.5) +
        geom_vline(xintercept = log2(1000)) +
        ggtitle(sprintf('%s widths (bins = 0.5)', qc_prefix)) +
        xlab('log2_width')
    ggsave(final_widths_gg, filename = sprintf('QC_%s_final_widths.pdf', qc_prefix), width = 8, height = 6)

    final_distances_gg = ggplot(final_nearest_distances_df, aes(distance)) +
      geom_histogram(binwidth = 100) +
      geom_vline(xintercept = 500) +
      scale_x_continuous(limits = c(-100, 1500)) +
      ggtitle(sprintf('%s nearest distances (bins = 100)', qc_prefix)) +
      xlab('Distance to nearest')
    ggsave(final_distances_gg, filename = sprintf('QC_%s_final_distances.pdf', qc_prefix), width = 8, height = 6)

    # Combine summary information
    qc_summary = rbind(qc_summary,
        data.frame(
            'enhancers' = qc_prefix,
            'num_pre_process' = qc_pre,
            'num_post_process' = qc_post,
            'num_final' = qc_final,
            'num_pre_5kb_overlaps_overall' = pre_overlaps_5kb_overall,
            'num_pre_5kb_overlaps_within' = pre_overlaps_5kb_within,
            'num_pre_5kb_overlaps_trim' = pre_overlaps_5kb_trim,
            'num_post_5kb_overlaps_overall' = post_overlaps_5kb_overall,
            'num_post_5kb_overlaps_within' = post_overlaps_5kb_within,
            'num_post_5kb_overlaps_trim' = post_overlaps_5kb_trim))
    #################################
}

write.table(qc_summary, file = sprintf('QC_summary_%s.txt', extension), sep='\t', quote = F, row.names = F, col.names = T)

#####  TEST expand_and_resect()  #####

#df <- data.frame(chromosome=rep("chr1", 6), start = c(57, 61, 81, 97, 122, 127), end = c(60, 72, 95, 98, 124, 129))
#df.ranges <- makeGRangesFromDataFrame(df)
#df.ranges
#expand_and_resect(df.ranges, min_size = 10)

#for(i in seq(2, 1000)) {
#  seqname <- paste("chr", as.character(i), sep = '')
#  df <- data.frame(chromosome=rep(seqname, 6), start = c(57, 61, 81, 97, 122, 127), end = c(60, 72, 95, 98, 124, 129))
#  df.ranges.new <- makeGRangesFromDataFrame(df)
#  df.ranges <- c(df.ranges, df.ranges.new)
#}
#expand_and_resect(df.ranges, min_size = 10)

# df <- data.frame(chromosome=rep("chr1", 4), start = c(5, 14, 20, 32), end = c(7, 17, 25, 40))
# df.ranges <- makeGRangesFromDataFrame(df)
# df.ranges
# expand_and_resect(df.ranges, min_size = 1)
#
#
# df <- data.frame(chromosome=rep("chr1", 4), start = c(5, 14, 20, 32), end = c(7, 17, 25, 40))
# df.ranges <- makeGRangesFromDataFrame(df)
# df.ranges
# expand_and_resect(df.ranges, min_size = 5)
#
#
# df <- data.frame(chromosome=rep("chr1", 4), start = c(5, 14, 20, 32), end = c(7, 17, 25, 40))
# df.ranges <- makeGRangesFromDataFrame(df)
# df.ranges
# expand_and_resect(df.ranges, min_size = 2000)
#
#
# df <- data.frame(chromosome=rep("chr1", 4), start = c(5, 6, 20, 32), end = c(6, 20, 25, 40))
# df.ranges <- makeGRangesFromDataFrame(df)
# df.ranges
# expand_and_resect(df.ranges, min_size = 10)

# df <- data.frame(chromosome=rep("chr1", 4), start = c(5, 8, 20, 32), end = c(6, 20, 25, 40))
# df.ranges <- makeGRangesFromDataFrame(df)
# df.ranges
# expand_and_resect(df.ranges, min_size = 10)
