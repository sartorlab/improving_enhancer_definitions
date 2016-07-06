# take care of options here
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("--results"), action = "store", type="character", help = "[Required] results.tab output from chipenrich.  Comma-separated list of files to be compared"),
  make_option(c("--labels"), action = "store", type="character", help = "[Required] labels for the pairs plots (one label for each of the results files passed)"),
  make_option(c("--out"), action = "store", help = "[Required] Name of the outfile")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list=option_list, add_help_option = T, prog = "chipenrich_pairs.R")
opts <- parse_args(option_parser)

#### ADD:
# color according to:
# if both have FDR less than .05, one has less than 0.05, or neither (get coloring scheme from chris)
# Negative/positive based on log odds ratio




#####  VERIFY OPTIONS  #####
if (length(opts$results)==0) {
  stop("Option --results not provided")
}
if (length(opts$labels)==0) {
  stop("Option --labels not provided")
}
if (length(opts$out)==0) {
  stop("Option --out not provided")
}

# For testing...

#opts <- list("chromhmm_dnase.2000.E.wgEncodeAwgTfbsBroadGm12878CtcfUniPk.narrowPeak_results.tab,chromhmm_dnase.2000.P2P.wgEncodeAwgTfbsBroadGm12878CtcfUniPk.narrowPeak_results.tab,current.wgEncodeAwgTfbsBroadGm12878CtcfUniPk.narrowPeak_results.tab",
 #            "DNase.2000.E,DNase.2000.P2P,current")
#names(opts) <- c("results", "labels")


results_files <- strsplit(opts$results, ",")[[1]]
if (length(results_files) < 2) {
  stop("At least two results files must be passed to --results (as a comma-separated list)")
}

labels <- strsplit(opts$labels, ",")[[1]]
if (length(labels) != length(results_files)) {
  stop("The number of labels must match the number of results files")
}


#####  LOAD IN ADDITIONAL PACKAGES  #####
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("tidyr"))

#####  READ IN THE RESULTS  ######
all <- ""

for(i in results_files) {
  tmp <- read.table(i, head = T, as.is = T, sep = "\t", quote = "")
  tmp <- tmp[,c("Geneset.ID", "Description", "P.value", "FDR", "Status")]
  tmp$P.value <- -1 * log10(tmp$P.value)
  tmp <- tidyr::unite(data = tmp, col = category, ... = P.value:Status, sep = "::")
  colnames(tmp) <- c("Geneset.ID", "Description", labels[match(i, results_files)])
  if(is.data.frame(all)) {
    all <- join(all, tmp)
  } else {
    all <- tmp
  }
}

#####  MAKE THE PLOTS  #####
# ditch unnecessary columns
all <- all[,!colnames(all) %in% c("Geneset.ID", "Description")]

# reshape the data.frame
all <- tidyr::gather(all, key = category, value, grep("current", colnames(all), invert = T))
all <- tidyr::separate(data = all, col = value, into = c("p.value", "fdr", "status"), sep = "::")
all <- tidyr::separate(data = all, col = current, into = c("current_p.value", "current_fdr", "current_status"), sep = "::")
all$current_p.value <- as.numeric(all$current_p.value)
all$current_fdr <- as.numeric(all$current_fdr)
all$p.value <- as.numeric(all$p.value)
all$fdr <- as.numeric(all$fdr)

# make a column that will be used for color coding (as follows):
# if at least one is not significant: black (000000)
# if both are enriched: A123F0
# if both are depleted: CD9D23
# if current is enriched and other is depleted: 17D017
# if current is depleted and other is enriched: 45E0D0
colors <- c("#000000", "#A123F0", "#CD9D23", "#17D017", "#45E0D0")
all$effect <- "ns"
all$effect[all$current_fdr <= 0.05 & all$fdr <= 0.05 & all$current_status=="enriched" & all$status=="enriched"] <- "enr/enr"
all$effect[all$current_fdr <= 0.05 & all$fdr <= 0.05 & all$current_status=="depleted" & all$status=="enriched"] <- "dep/enr"
all$effect[all$current_fdr <= 0.05 & all$fdr <= 0.05 & all$current_status=="depleted" & all$status=="depleted"] <- "dep/dep"
all$effect[all$current_fdr <= 0.05 & all$fdr <= 0.05 & all$current_status=="enriched" & all$status=="depleted"] <- "enr/dep"

# need to know the maximum values so that the plot limits can be set accordingly...
max_val <- max(c(max(all$current_p.value), max(all$p.value)))

# change the sign of the p-value depending on whether or not it's enriched or depleted...
all$current_p.value[all$current_status=="depleted"] <- -1 * all$current_p.value[all$current_status=="depleted"]
all$p.value[all$status=="depleted"] <- -1 * all$p.value[all$status=="depleted"]

# plot
png(filename = opts$out, width = 10, height = 7, units = "in", res = 300)
p <- ggplot(all) + geom_point(aes(x = current_p.value, y = p.value, color = effect))
p <- p + facet_wrap( ~ category)
p <- p + geom_abline(color = "red", linetype = "dashed", slope = 1, intercept = 0) 
p <- p + xlim(c(max_val*-1.1, max_val*1.1)) + ylim(c(max_val*-1.1, max_val*1.1))
p <- p + scale_colour_manual(values=c("ns" = colors[1], "enr/enr" = colors[2], "dep/dep" = colors[3], "enr/dep" = colors[4], "dep/enr" = colors[5]))
p
