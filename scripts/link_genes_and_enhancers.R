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
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("hash"))


#####  MAIN  #####

# read in the genes loci
genes <- read.table(opts$gene_loci, head = F, as.is = T, sep = '\t')
colnames(genes) <- c("chromosome", "start", "end", "gene_id")

# read in the list of enhancers
enhancers <- read.table(opts$enhancers, head = F, as.is = T, sep = '\t')
colnames(enhancers) <- c("chromosome", "start", "end")
enhancers$enhancer_id <- paste(enhancers$chromosome, enhancers$start, enhancers$end, sep = ':')

# read in the interactions
interactions <- read.table(opts$interactions, head = F, as.is = T, sep = '\t')
colnames(interactions) <- c("chromosome_1", "start_1", "end_1", "chromosome_2", "start_2", "end_2")
interactions$region_id_1 <- paste(interactions$chromosome_1, interactions$start_1, interactions$end_1, sep = ':')
interactions$region_id_2 <- paste(interactions$chromosome_2, interactions$start_2, interactions$end_2, sep = ':')

# only interested in intra-chromosomal interactions
interactions <- interactions[interactions$chromosome_1==interactions$chromosome_2,]


if (opts$method == "point_to_point") {
  # generate the entire list of interaction regions
  first_interaction_regions <- interactions[,c("chromosome_1", "start_1", "end_1", "region_id_1")]
  second_interaction_regions <- interactions[,c("chromosome_2", "start_2", "end_2", "region_id_2")]
  colnames(first_interaction_regions) <- c("chromosome", "start", "end", "region_id")
  colnames(second_interaction_regions) <- c("chromosome", "start", "end", "region_id")
  interaction_regions <- rbind(first_interaction_regions, second_interaction_regions)
  interaction_regions <- unique(interaction_regions)
  interaction_regions <- as.data.frame(interaction_regions) 
  
  # find the interaction regions that overlap with genes
  interaction_regions.range <- makeGRangesFromDataFrame(df = interaction_regions, seqnames.field = "chromosome", start.field = "start", end.field = "end")
  genes.range <- makeGRangesFromDataFrame(df = genes, seqnames.field = "chromosome", start.field = "start", end.field = "end")
  gene_interaction.overlaps <- findOverlaps(genes.range, interaction_regions.range)
  gene_interaction.overlaps <- cbind(genes[queryHits(gene_interaction.overlaps),"gene_id"], interaction_regions[subjectHits(gene_interaction.overlaps),"region_id"])
  gene_interaction.overlaps <- as.data.frame(gene_interaction.overlaps, stringsAsFactors = F)
  colnames(gene_interaction.overlaps) <- c("gene_id", "region_id")
  gene_interaction.overlaps <- unique(gene_interaction.overlaps)
  
  
  # find the interaction edges that overlap with enhancers
  interaction_regions.range <- makeGRangesFromDataFrame(df = interaction_regions, seqnames.field = "chromosome", start.field = "start", end.field = "end")
  enhancers.range <- makeGRangesFromDataFrame(df = enhancers, seqnames.field = "chromosome", start.field = "start", end.field = "end")
  enhancer_interaction.overlaps <- findOverlaps(enhancers.range, interaction_regions.range)
  enhancer_interaction.overlaps <- cbind(enhancers[queryHits(enhancer_interaction.overlaps),"enhancer_id"], interaction_regions[subjectHits(enhancer_interaction.overlaps),"region_id"])
  enhancer_interaction.overlaps <- as.data.frame(enhancer_interaction.overlaps, stringsAsFactors = F)
  colnames(enhancer_interaction.overlaps) <- c("enhancer_id", "region_id")
  enhancer_interaction.overlaps <- unique(enhancer_interaction.overlaps)
  
  # now assign the gene-enhancer interactions based on the above overlaps
  
  # populate hash, gene_interactions[["region_id"]] --> c(genes_overlapping_this_region)
  gene_interactions <- hash(keys = gene_interaction.overlaps$region_id)
  for(i in seq(1, nrow(gene_interaction.overlaps))){
    gene_interactions[[gene_interaction.overlaps$region_id[i]]] <- c(gene_interactions[[gene_interaction.overlaps$region_id[i]]], gene_interaction.overlaps$gene_id[i])
  }
  
  # populate hash, enhancer_interactions[["region_id"]] --> c(enhancers_overlapping_this_region)
  enhancer_interactions <- hash(keys = enhancer_interaction.overlaps$region_id)
  for(i in seq(1, nrow(enhancer_interaction.overlaps))){
    enhancer_interactions[[enhancer_interaction.overlaps$region_id[i]]] <- c(enhancer_interactions[[enhancer_interaction.overlaps$region_id[i]]], enhancer_interaction.overlaps$enhancer_id[i])
  }
  
  gene_in_pair <- c()
  enhancer_in_pair <- c()
  
  for(row in seq(1, nrow(interactions))) {
    first_region_id <- interactions$region_id_1[row]
    second_region_id <- interactions$region_id_2[row]
    if (has.key(first_region_id, gene_interactions) & has.key(second_region_id, enhancer_interactions)) {
      for(gene in gene_interactions[[first_region_id]]){
        for(enhancer in enhancer_interactions[[second_region_id]]){
          gene_in_pair <- c(gene_in_pair, gene)
          enhancer_in_pair <- c(enhancer_in_pair, enhancer)
        }
      }
    }
    
    if (has.key(second_region_id, gene_interactions) & has.key(first_region_id, enhancer_interactions)) {
      for(gene in gene_interactions[[second_region_id]]){
        for(enhancer in enhancer_interactions[[first_region_id]]){
          gene_in_pair <- c(gene_in_pair, gene)
          enhancer_in_pair <- c(enhancer_in_pair, enhancer)
        }
      }
    }
  }
  
  gene_enhancer_pairs <- data.frame(gene=gene_in_pair, enhancer=enhancer_in_pair)
  write.table(x = gene_enhancer_pairs, file = opts$out, append = F, quote = F, sep = "\t", row.names = F, col.names = T)

} else if (opts$method == "encompassing") {
  
  interactions$region_chromosome <- interactions$chromosome_1
  interactions$region_start <- apply(interactions[,c("start_1", "start_2")], 1, function(x){min(x)})
  interactions$region_end <- apply(interactions[,c("end_1", "end_2")], 1, function(x){max(x)})
  interactions$region_id <- paste(interactions$region_chromosome, interactions$region_start, interactions$region_end, sep=":")
  
  # find the genes overlapping each region 
  interaction_regions.range <- makeGRangesFromDataFrame(df = interactions[,c("region_chromosome", "region_start", "region_end")], seqnames.field = "region_chromosome", start.field = "region_start", end.field = "region_end")
  genes.range <- makeGRangesFromDataFrame(df = genes, seqnames.field = "chromosome", start.field = "start", end.field = "end")
  gene_interaction.overlaps <- findOverlaps(genes.range, interaction_regions.range)
  gene_interaction.overlaps <- cbind(genes[queryHits(gene_interaction.overlaps),"gene_id"], interactions[subjectHits(gene_interaction.overlaps),"region_id"])
  gene_interaction.overlaps <- as.data.frame(gene_interaction.overlaps, stringsAsFactors = F)
  colnames(gene_interaction.overlaps) <- c("gene_id", "region_id")
  gene_interaction.overlaps <- unique(gene_interaction.overlaps)
  
  if (length(opts$max_genes_per_region) > 0) {
    # a threshold has been set for the maximum number of genes per region
    number_genes_per_region <- table(gene_interaction.overlaps$region_id)
    exclude <- names(number_genes_per_region[number_genes_per_region>opts$max_genes_per_region])
    gene_interaction.overlaps <- gene_interaction.overlaps[!gene_interaction.overlaps$region_id %in% exclude,]
  }
  
  
  # find the enhancers overlapping each region
  interaction_regions.range <- makeGRangesFromDataFrame(df = interactions[,c("region_chromosome", "region_start", "region_end")], seqnames.field = "region_chromosome", start.field = "region_start", end.field = "region_end")
  enhancers.range <- makeGRangesFromDataFrame(df = enhancers, seqnames.field = "chromosome", start.field = "start", end.field = "end")
  enhancer_interaction.overlaps <- findOverlaps(enhancers.range, interaction_regions.range)
  enhancer_interaction.overlaps <- cbind(enhancers[queryHits(enhancer_interaction.overlaps),"enhancer_id"], interactions[subjectHits(enhancer_interaction.overlaps),"region_id"])
  enhancer_interaction.overlaps <- as.data.frame(enhancer_interaction.overlaps, stringsAsFactors = F)
  colnames(enhancer_interaction.overlaps) <- c("enhancer_id", "region_id")
  enhancer_interaction.overlaps <- unique(enhancer_interaction.overlaps)
  
  
  # now assign the gene-enhancer interactions based on the above overlaps
  gene_enhancer_pairs <- join(gene_interaction.overlaps, enhancer_interaction.overlaps)
  gene_enhancer_pairs <- gene_enhancer_pairs[!is.na(gene_enhancer_pairs$enhancer_id),]  
  gene_enhancer_pairs <- gene_enhancer_pairs[,c("gene_id", "enhancer_id")]
  gene_enhancer_pairs <- unique(gene_enhancer_pairs)
  
  write.table(x = gene_enhancer_pairs, file = opts$out, append = F, quote = F, sep = "\t", row.names = F, col.names = T)
}

