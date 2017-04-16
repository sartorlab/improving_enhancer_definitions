suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(c("--genes"), action = "store", type = "character", help = "[Required] Path to a file listing the genomic intervals that should be assigned to each gene.  Format is one line per locus, tab-separated, chromosome, start, end, gene_id (no header)"),
  make_option(c("--dnase"), action = "store", type = "character", help = "[Optional] Path to a file listing the enhancers.  Format is one line per locus, tab-separated, chromosome, start, end, enhancer_id (no header)"),
  make_option(c("--tissue_threshold_dnase"), action = "store", type = "numeric", default = 1, help = "[Optional] An integer.  If the DNase list is used to create the enhancer list, then only DNase sites supported by this number of experiments are included.  There are ~125 experiments in the ENCODE DNase list (Default: 1)."),
  make_option(c("--chromhmm_directory"), action = "store", type = "character", help = "[Optional] Path to the directory containing the ENCODE chromHMM files. The file names should end in '*HMM.bed.gz'"),
  make_option(c("--fantom"), action = "store", type = "character", help = "[Optional] Path to a file listing enhancers from the FANTOM5 consortium."),
  make_option(c("--thurman"), action = "store", type = "character", help = "[Optional] Path to a file listing enhancers (and interactions) from Thurman et al. paper."),
  make_option(c("--extension"), action = "store", type = "numeric", help = "[Required] Number of base pairs to extend enhancers to"),
)

option_parser = OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T,
                              description = "\nDescription: Accepts a list of gene loci, enhancers, and an interaction dataset (e.g., ChIA-PET or Hi-C).  Prints a list of gene-enhancer pairs. Generates this list using one of two methods: either by finding the interactions that link one gene to one enhancer (each end of the interaction overlaps with one of the genes/enhancers; 'point_to_point') or by declaring all gene/enhancer pairs encompassed WITHIN an interaction region to be interacting ('encompassing'; useful when e.g. enhancers/genes are expected to be within the loops rather than at the interaction edges, such as would be expected within a TAD).")
opts = parse_args(option_parser)

#####  VERIFY THE INPUT  #####
if (length(opts$genes) == 0) {
  stop("--genes is a required argument")
}
if( (length(opts$dnase) + length(opts$chromhmm_directory) + length(opts$fantom) + length(opts$thurman)) != 4) {
  stop("All of --dnase, --chromhmm_directory, --fantom, or --thurman must be given")
}
if (length(opts$extension) == 0) {
  stop("--extension is a required argument")
}
if (length(opts$out) == 0) {
  stop("--out is a required argument")
}

#####  LOAD THE REST OF THE LIBRARIES  #####
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("hash"))
source("../scripts/create_enhancer_list.R", chdir = T)

#####  MAIN  #####

#####  Read gene loci  #####
genes = read.table(opts$genes, header = F, as.is = T, sep = '\t')
colnames(genes) = c("chromosome", "start", "end", "gene_id", "symbol")
genes_gr = makeGRangesFromDataFrame(df = genes, seqnames.field = "chromosome", start.field = "start", end.field = "end")

#####  Load enhancer location files as given  #####
dnase_pool = c(NA, 'dnase')
chromhmm_pool = c(NA, 'chromhmm')
fantom_pool = c(NA, 'fantom')
thurman_pool = c(NA, 'thurman')

combinations = expand.grid(dnase_pool, chromhmm_pool, fantom_pool, thurman_pool, ext_pool, stringsAsFactors=F)
combinations = combinations[!apply(combinations, 1, function(row){all(is.na(row))}),]
combinations = combinations[!apply(combinations, 1, function(row){'dnase' %in% row && 'thurman' %in% row}), ]

for(i in 1:nrow(combinations)) {
  dnase_code = combinations[i,1]
  chromhmm_code = combinations[i,2]
  fantom_code = combinations[i,3]
  thurman_code = combinations[i,4]

  message(sprintf('On %s %s %s %s %s', dnase_code, chromhmm_code, fantom_code, thurman_code, opts$extension))

  enhancers_df = data.frame()
  out_code = ''

  if (!is.na(chromhmm_code)) {
    # load in the chromhmm data for all the tissues
    chrom_hmm_files = list.files(opts$chromhmm_directory, pattern = "HMM.bed.gz")
    chrom_hmm = ""
    for(file in chrom_hmm_files) {
      cell_type = gsub("wgEncodeBroadHmm(.*)HMM.bed.gz", "\\1", file, perl = T)
      dat = read.table(gzfile(paste(opts$chromhmm_directory, file, sep='/')), header = F, as.is = T, sep = '\t')
      colnames(dat) = c("chromosome", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb")
      dat = dat[,c("chromosome", "start", "end", "name")]
      dat$cell_type = cell_type
      if (is.data.frame(chrom_hmm)) {
        chrom_hmm = rbind(chrom_hmm, dat)
      } else {
        chrom_hmm = dat
      }
    }

    chrom_hmm = chrom_hmm[grep("Enhancer", chrom_hmm$name),]
    chrom_hmm = chrom_hmm[,c("chromosome", "start", "end")]

    enhancers_df = rbind(enhancers_df, chrom_hmm)
    out_code = c(out_code, chromhmm_code)
  }

  if (!is.na(dnase_code)) {
    # load in the dnase data
    dnase = read.table(gzfile(opts$dnase), header = F, as.is = T, sep = '\t')
    colnames(dnase) = c("chromosome", "start", "end", "name", "score", "floatScore", "sourceCount", "sourceIds")

    dnase = dnase[dnase$sourceCount >= opts$tissue_threshold_dnase, ]

    enhancers_df = rbind(enhancers_df, dnase[,c("chromosome", "start", "end")])
    out_code = c(out_code, dnase_code)
  }

  if (!is.null(thurman_code)) {
    thurman = read.table(opts$dnase, header = F, as.is = T, sep = '\t')
    colnames(thurman) = c('promoter.chromosome', 'promoter.start', 'promoter.end', 'promoter.symbol', 'distal.chromosome', 'distal.start', 'distal.end', 'cor')

    thurman = thurman[, c('distal.chromosome', 'distal.start', 'distal.end')]
    colnames(thurman) = c('chromosome', 'start', 'end')

    enhancers_df = rbind(enhancers_df, thurman)
    out_code = c(out_code, thurman_code)
  }

  if (!is.na(fantom_code)) {
    fantom = read.table(opts$fantom, header = F, as.is = T, sep = '\t')
    fantom = fantom[, c(1:3)]
    colnames(fantom) = c('chromosome', 'start', 'end')

    enhancers_df = rbind(enhancers_df, fantom)
    out_code = c(out_code, fantom_code)
  }

  enhancers_gr = makeGRangesFromDataFrame(enhancers_df)
  enhancers_merged = reduce(enhancers_gr)

  enhancers_merged = create_enhancer_list(enhancers_merged, opts$extension)

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

  out_file = paste(out_code, collapse = '_')
  out_file = paste(out_file, extension, 'enhancers', collapse = '.')

  write.table(enhancers_df[,c("seqnames", "start", "end")], file = out_file, quote = F, sep = "\t", row.names = F, col.names = F)
}
