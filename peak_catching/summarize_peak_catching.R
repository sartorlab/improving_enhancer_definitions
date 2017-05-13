library(dplyr)
library(ggplot2)
# library(plotly)

# Peak catching calculated by Tingting
    tt_peak_catching = read.table('genomeCoverage_peakCatching_enhNum.txt', header=T, sep='\t', stringsAsFactors=F)

# Number of peaks in the experiment
    num_peaks = read.table('num_peaks.txt', header=F, sep=' ', stringsAsFactors=F)
    colnames(num_peaks) = c('tf', 'num_peaks')
    num_peaks$tf = gsub('.narrowPeak','',num_peaks$tf)

# Number of unique peaks assigned with the locus definition
    num_uniq_peaks = read.table('num_uniq_peaks_assigned.txt', header=F, sep=' ', stringsAsFactors=F)
    colnames(num_uniq_peaks) = c('tf', 'ldef', 'num_caught')
    num_uniq_peaks$ldef = sapply(num_uniq_peaks$ldef, function(l){gsub('_','.',l)})

    num_uniq_peaks = subset(num_uniq_peaks, tf != 'wgEncodePilotTfbsMyersA549GrDex')

# Read in chromosomze sizes
    chrom_sizes = read.table('~/espresso/share/genomes/hg19/chromInfo_hg19.txt', sep='\t', header=F, as.is = TRUE)
    genome_size = sum(as.numeric(chrom_sizes[,2]))

# Read in genome coverage
    genome_coverage = read.table('genome_coverage.txt', header=F, sep = ' ', stringsAsFactors=F)
    colnames(genome_coverage) = c('ldef', 'bp_covered')
    genome_coverage$ldef = sapply(genome_coverage$ldef, function(l){gsub('_','.',l)})
    genome_coverage$ldef = gsub('.ldef','',genome_coverage$ldef)
    genome_coverage$prop_covered = genome_coverage$bp_covered / genome_size

# Join num_peaks to num_uniq_peaks on the tf
num_uniq_peaks$num_peaks = num_peaks[match(num_uniq_peaks$tf, num_peaks$tf), 'num_peaks']

# Per experiment, what is the percentage of peaks caught?
num_uniq_peaks$prop_caught = num_uniq_peaks$num_caught / num_uniq_peaks$num_peaks

# Join genome coverage
num_uniq_peaks$prop_coverage = genome_coverage[match(num_uniq_peaks$ldef, genome_coverage$ldef), 'prop_covered']

# Overall peak catching per ldef
peak_catching_ldef = summarize(group_by(num_uniq_peaks, ldef), mean_chipenrich_caught = mean(prop_caught), prop_covered = mean(prop_coverage), n = n())

# A tibble: 6 × 3
#                          ldef mean_chipenrich_caught     n
#                         <chr>                  <dbl> <int>
# 1                chromhmm.0.E             0.09788291   298
# 2         chromhmm.0.E.fantom             0.23259272   298
# 3        chromhmm.0.E.thurman             0.21000799   298
# 4 chromhmm.0.E.thurman.fantom             0.25913581   298
# 5           chromhmm.0.fantom             0.19388853   298
# 6              chromhmm.0.P2P             0.10927716   298

# Compare to what Heming calculated based on findOverlaps().
# NOTE: The difference is findOverlaps() considers any overlap whereas chipenrich requires the midpoint be in the locus
combined = merge(peak_catching_ldef, tt_peak_catching, by.x='ldef', by.y='dataset')
colnames(combined) = c('ldef','mean_chipenrich_caught','genome_coverage_raymond','n_datasets','genome_coverage_heming','mean_findOverlap_caught','numEnh_ge1Gene','genome_coverage_tingting')

# Write tabular output
write.table(combined, file='peak_catching_compare_chipenrich_findoverlaps.txt', sep='\t', row.names=F, col.names=T, quote=F)

caught_by_coverage_raymond = ggplot(data = combined, aes(x=mean_chipenrich_caught, y=genome_coverage_raymond, text = paste('ldef:', ldef))) +
    geom_point(alpha=0.5) + xlab('Peaks Caught (by midpoint)') + ylab('Genome Coverage (Raymond)')
ggsave(filename='midpoint_peak_catching_genome_coverage_raymond.pdf', plot = caught_by_coverage_raymond, width=6, height=6)

caught_by_coverage = ggplot(data = combined, aes(x=mean_chipenrich_caught, y=genome_coverage_tingting, text = paste('ldef:', ldef))) +
    geom_point(alpha=0.5) + xlab('Peaks Caught (by midpoint)') + ylab('Genome Coverage (Ting Ting)')
ggsave(filename='midpoint_peak_catching_genome_coverage.pdf', plot = caught_by_coverage, width=6, height=6)

midpoint_by_findoverlaps = ggplot(data = combined, aes(x=mean_chipenrich_caught, y=mean_findOverlap_caught, text = paste('ldef:', ldef))) +
    geom_point(alpha=0.5) + xlab('Peaks Caught (by midpoint)') + ylab('Peaks Caught (by findOverlaps)') +
    geom_abline(intercept = 0, slope = 1)
ggsave(filename='peak_catching_compare_chipenrich_findoverlaps.pdf', plot = midpoint_by_findoverlaps, width=6, height=6)

coverage_comparison = ggplot(data = combined, aes(x = genome_coverage_tingting, y = genome_coverage_raymond), text = paste('ldef:', ldef)) +
    geom_point(alpha=0.5) + xlab('Genome Coverage (old ldef)') + ylab('Genome Coverage (new ldef)')
ggsave(filename='genome_coverage_comparison.pdf', plot = caught_by_coverage, width=6, height=6)

# ggplotly(caught_by_coverage)
# ggplotly(midpoint_by_findoverlaps)
