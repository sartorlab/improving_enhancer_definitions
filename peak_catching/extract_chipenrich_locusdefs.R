library(chipenrich.data)

data('locusdef.hg19.5kb', package = 'chipenrich.data')
data('locusdef.hg19.5kb_outside', package = 'chipenrich.data')
data('locusdef.hg19.nearest_tss', package = 'chipenrich.data')

fivekb_df = locusdef.hg19.5kb@dframe
fivekb_outside_df = locusdef.hg19.5kb_outside@dframe
ntss_df = locusdef.hg19.nearest_tss@dframe

write.table(fivekb_df, file = gzfile('5kb.ldef'), sep='\t', col.names = T, row.names = F, quote = F
write.table(fivekb_outside_df, file = gzfile('5kb_outside.ldef'), sep='\t', col.names = T, row.names = F, quote = F
write.table(ntss_df, file = gzfile('nearest_tss.ldef'), sep='\t', col.names = T, row.names = F, quote = F)
