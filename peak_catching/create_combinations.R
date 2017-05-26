chipseq_path = '~/espresso/share/ENCODE_chipseq'
ldef_path = '~/latte/share/improving_enhancer_definitions/locusdefs'

chipseq_files = list.files(path = chipseq_path, pattern='.narrowPeak', full.names = TRUE)
ldef_files = list.files(path = ldef_path, pattern='.ldef', full.names = TRUE)

combinations = expand.grid(ldef_files, chipseq_files, stringsAsFactors=F)

write.table(combinations, file = 'combinations.txt', sep=' ', col.names = F, row.names = F, quote = F)
