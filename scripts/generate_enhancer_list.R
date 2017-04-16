suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(c("--genes"), action = "store", type = "character", help = "[Required] Path to a file listing the genomic intervals that should be assigned to each gene.  Format is one line per locus, tab-separated, chromosome, start, end, gene_id (no header)"),
  make_option(c("--dnase"), action = "store", type = "character", help = "[Optional] Path to a file listing the enhancers.  Format is one line per locus, tab-separated, chromosome, start, end, enhancer_id (no header)"),
  make_option(c("--tissue_threshold_dnase"), action = "store", type = "numeric", default = 1, help = "[Optional] An integer.  If the DNase list is used to create the enhancer list, then only DNase sites supported by this number of experiments are included.  There are ~125 experiments in the ENCODE DNase list (Default: 1)."),
  make_option(c("--chromhmm"), action = "store", type = "character", help = "[Optional] File of concatenated chromHMM tracks"),
  make_option(c("--fantom"), action = "store", type = "character", help = "[Optional] Path to a file listing enhancers from the FANTOM5 consortium."),
  make_option(c("--thurman"), action = "store", type = "character", help = "[Optional] Path to a file listing enhancers (and interactions) from Thurman et al. paper."),
  make_option(c("--extension"), action = "store", type = "numeric", help = "[Required] Number of base pairs to extend enhancers to")
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
source("../scripts/create_enhancer_list.R", chdir = T)

#####  MAIN  #####

#####  READ GENE LOCI  #####
genes = read.table(opts$genes, header = F, as.is = T, sep = '\t')
colnames(genes) = c("chromosome", "start", "end", "gene_id", "symbol")
genes_gr = makeGRangesFromDataFrame(df = genes, seqnames.field = "chromosome", start.field = "start", end.field = "end")

#####  DETERMINE COMBINATIONS  #####
dnase_pool = c(NA, 'dnase')
chromhmm_pool = c(NA, 'chromhmm')
fantom_pool = c(NA, 'fantom')
thurman_pool = c(NA, 'thurman')

combinations = expand.grid(dnase_pool, chromhmm_pool, fantom_pool, thurman_pool, stringsAsFactors=F)
# Get rid of all NA
combinations = combinations[!apply(combinations, 1, function(row){all(is.na(row))}),]
# Get rid of combinations with dnase and thurman
combinations = combinations[!apply(combinations, 1, function(row){'dnase' %in% row && 'thurman' %in% row}), ]

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

#####  BUILD THE COMBINATIONS  #####
for(i in 1:nrow(combinations)) {
  dnase_code = combinations[i,1]
  chromhmm_code = combinations[i,2]
  fantom_code = combinations[i,3]
  thurman_code = combinations[i,4]

  message(sprintf('On %s %s %s %s %s', dnase_code, chromhmm_code, fantom_code, thurman_code, opts$extension))

  enhancers_df = data.frame()
  out_code = c()

  if (!is.na(chromhmm_code)) {
    message('Adding chromhmm')
    enhancers_df = rbind(enhancers_df, chromhmm)
    out_code = c(out_code, chromhmm_code)
  }

  if (!is.na(dnase_code)) {
    message('Adding dnase')
    enhancers_df = rbind(enhancers_df, dnase[,c("chromosome", "start", "end")])
    out_code = c(out_code, dnase_code)
  }

  if (!is.na(thurman_code)) {
    message('Adding thurman')
    enhancers_df = rbind(enhancers_df, thurman)
    out_code = c(out_code, thurman_code)
  }

  if (!is.na(fantom_code)) {
    message('Adding fantom')
    enhancers_df = rbind(enhancers_df, fantom)
    out_code = c(out_code, fantom_code)
  }

  enhancers_gr = makeGRangesFromDataFrame(enhancers_df)

  enhancers_merged = create_enhancer_list(enhancers_gr, opts$extension)

  # ditch things that overlap with genes...
  overlaps = findOverlaps(genes_gr, enhancers_merged)
  kick_out = unique(subjectHits(overlaps))
  enhancers_merged = enhancers_merged[-1 * kick_out]

  enhancers_df = as.data.frame(enhancers_merged)
  enhancers_df = enhancers_df[order(enhancers_df$seqnames, enhancers_df$start), ]

  if (opts$extension == 1) {
    extension = 0
  } else {
    extension = opts$extension
  }

  out_file = paste(paste(out_code, collapse = '_'), extension, 'enhancers', 'gz', sep = '.')

  message(sprintf('Writing %s', out_file))
  write.table(enhancers_df[,c("seqnames", "start", "end")], file = gzfile(out_file), quote = F, sep = "\t", row.names = F, col.names = F)
}
