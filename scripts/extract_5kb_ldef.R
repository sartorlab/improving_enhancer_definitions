# NOTE: Should be chipenrich.data 1.99.0 or greater for new locus definitions
library(chipenrich.data)

data('locusdef.hg19.5kb', package = 'chipenrich.data')

df = locusdef.hg19.5kb@dframe

write.table(df, file = 'chipenrich_5kb_locusdef.txt', sep = '\t', quote = F, col.names = F, row.names = F)
