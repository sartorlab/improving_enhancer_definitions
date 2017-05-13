suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(c("--genes"), action = "store", type = "character", help = "[Required] Path to a file listing the genomic intervals that should be assigned to each gene.  Format is one line per locus, tab-separated, chromosome, start, end, gene_id (no header)"),
  make_option(c("--chromhmm"), action = "store", type = "character", help = "[Optional] File of concatenated chromHMM tracks"),
  make_option(c("--dnase"), action = "store", type = "character", help = "[Optional] Path to a file listing the enhancers.  Format is one line per locus, tab-separated, chromosome, start, end, enhancer_id (no header)"),
  make_option(c("--thurman"), action = "store", type = "character", help = "[Optional] Path to a file listing enhancers (and interactions) from Thurman et al. paper."),
  make_option(c("--fantom"), action = "store", type = "character", help = "[Optional] Path to a file listing enhancers from the FANTOM5 consortium."),
  make_option(c("--extension"), action = "store", type = "numeric", help = "[Required] Number of base pairs to extend enhancers to"),  make_option(c("--tissue_threshold_dnase"), action = "store", type = "numeric", default = 1, help = "[Optional] An integer.  If the DNase list is used to create the enhancer list, then only DNase sites supported by this number of experiments are included.  There are ~125 experiments in the ENCODE DNase list (Default: 1).")
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
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("plyr"))

#####  FUNCTION TO EXPAND AND RESECT GRANGES  #####
expand_and_resect <- function(granges, min_size) {
  # extends each grange to the min_size (if not there already); if this would cause two neighboring ranges to overlap, it extends them only until they bookend

  # make sure everything overlapping already has been merged
  granges <- reduce(granges, min.gapwidth = 0)

  granges$left = start(granges)
  granges$right = end(granges)


  # resize the ranges that are smaller than min_size; keep the ranges that are bigger than min_size as is.
  was_resized <- rep(F, length(granges))
  was_resized[(end(granges) - start(granges)) < min_size] <- T

  granges[(end(granges) - start(granges)) < min_size] <- resize(granges[(end(granges) - start(granges)) < min_size], width = min_size, fix = "center")
  resized_starts <- start(granges)
  resized_ends <- end(granges)

  # find the cases where, after resizing, two ranges overlap
  overlaps <- findOverlaps(granges, granges)
  overlaps <- overlaps[(queryHits(overlaps)+1)==subjectHits(overlaps),]


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

  out_file = paste(paste(out_code, collapse = '_'), extension, 'enhancers', 'gz', sep = '.')

  message(sprintf('On %s', out_file))

  if(file.exists(out_file)) {
    message('File exists, skipping...')
    next
  }

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

  enhancers_gr = makeGRangesFromDataFrame(enhancers_df, ignore.strand = TRUE)

  # Expand enhancers, but make adjacent enhnacers within extension just bump up
  # against each other at the midpoint
  enhancers_merged = expand_and_resect(enhancers_gr, opts$extension)

  # Remove any part of genes_gr from the enhancers_gr
  enhancers_merged = GenomicRanges::setdiff(enhancers_merged, genes_gr)

  enhancers_df = as.data.frame(enhancers_merged)
  enhancers_df = enhancers_df[order(enhancers_df$seqnames, enhancers_df$start), ]

  message(sprintf('Writing %s', out_file))
  write.table(enhancers_df[,c("seqnames", "start", "end")], file = gzfile(out_file), quote = F, sep = "\t", row.names = F, col.names = F)
}

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
