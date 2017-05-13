library(dplyr)
library(ggplot2)
# library(plotly)

# Peak catching calculated by Tingting
    tt_peak_catching = read.table('genomeCoverage_peakCatching_enhNum.txt', header=T, sep='\t', stringsAsFactors=F)

# Number of peaks in the experiment
    num_peaks = read.table('num_peaks.txt', header=F, sep=' ', stringsAsFactors=F)
    colnames(num_peaks) = c('tf', 'num_peaks')

# Number of unique peaks assigned with the locus definition
    num_uniq_peaks = read.table('num_uniq_peaks_assigned.txt', header=F, sep=' ', stringsAsFactors=F)
    colnames(num_uniq_peaks) = c('tf', 'ldef', 'num_caught')
    num_uniq_peaks$ldef = sapply(num_uniq_peaks$ldef, function(l){gsub('_','.',l)})

    num_uniq_peaks = subset(num_uniq_peaks, tf != 'wgEncodePilotTfbsMyersA549GrDex')

# Join num_peaks to num_uniq_peaks on the tf
num_uniq_peaks$num_peaks = num_peaks[match(num_uniq_peaks$tf, num_peaks$tf), 'num_peaks']

# Per experiment, what is the percentage of peaks caught?
num_uniq_peaks$prop_caught = num_uniq_peaks$num_caught / num_uniq_peaks$num_peaks

# Overall peak catching per ldef
peak_catching_ldef = summarize(group_by(num_uniq_peaks, ldef), mean_chipenrich_caught = mean(prop_caught), n = n())

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
colnames(combined) = c('ldef','mean_chipenrich_caught','n_datasets','genome_coverage_heming','mean_findOverlap_caught','numEnh_ge1Gene','genome_coverage_tingting')

# Write tabular output
write.table(combined, file='peak_catching_compare_chipenrich_findoverlaps.txt', sep='\t', row.names=F, col.names=T, quote=F)

caught_by_coverage = ggplot(data = combined, aes(x=mean_chipenrich_caught, y=genome_coverage_tingting, text = paste('ldef:', ldef))) +
    geom_point(alpha=0.5) + xlab('Peaks Caught (by midpoint)') + ylab('Genome Coverage')
ggsave(filename='midpoint_peak_catching_genome_coverage.pdf', plot = caught_by_coverage, width=6, height=6)

midpoint_by_findoverlaps = ggplot(data = combined, aes(x=mean_chipenrich_caught, y=mean_findOverlap_caught, text = paste('ldef:', ldef))) +
    geom_point(alpha=0.5) + xlab('Peaks Caught (by midpoint)') + ylab('Peaks Caught (by findOverlaps)') +
    geom_abline(intercept = 0, slope = 1)
ggsave(filename='peak_catching_compare_chipenrich_findoverlaps.pdf', plot = midpoint_by_findoverlaps, width=6, height=6)

# ggplotly(caught_by_coverage)
# ggplotly(midpoint_by_findoverlaps)
