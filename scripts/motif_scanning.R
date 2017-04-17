#find motif
library(optparse)
library(GenomicRanges)
library(ggplot2)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

option_list <- list(
  make_option(c("--input"), action = "store", type = "character"),
  make_option(c("--pwm"), action = "store", type = "character")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list)

opts <- parse_args(option_parser)
if (length(opts$input) == 0) {
  stop("--input is a required argument")
}

################################################################################
# Read interactions, create GRanges from 1-3 and 4-6, and extract sequences

seq_data = read.delim(opts$input, header=FALSE, stringsAsFactors=FALSE)
seq_data_length = nrow(seq_data)

seq1_gr = GenomicRanges::GRanges(
  seqnames = seq_data[,1],
  ranges = IRanges::IRanges(start = seq_data[,2], end = seq_data[,3])
)

seq2_gr = GenomicRanges::GRanges(
  seqnames = seq_data[,4],
  ranges = IRanges::IRanges(start = seq_data[,5], end = seq_data[,6])
)

seq1_seqs = BSgenome::getSeq(Hsapiens, names = seq1_gr)
seq2_seqs = BSgenome::getSeq(Hsapiens, names = seq2_gr)

################################################################################
# Build the PWM

pfm <- read.table(opts$pwm, row.names=1, quote="\"", comment.char="")
#pfm <- read.table(opts$pwm, row.names=1, quote="\"", comment.char="")
pfm = as.matrix(pfm)
x = as.integer(pfm)
pfm = matrix(x, nrow=4, byrow=FALSE)
rownames(pfm) = c("A","C","G","T")
pwm = Biostrings::PWM(pfm)

################################################################################
# Determine directionality of motifs

threshold = 0.87

result_convergent = rep.int(0, seq_data_length)
result_divergent = rep.int(0, seq_data_length)
result_samedirection = rep.int(0, seq_data_length)
result_one = rep.int(0, seq_data_length)
result_two = rep.int(0, seq_data_length)

for (i in seq_along(seq1_gr)) {
  sequence_1 = seq1_seqs[[i]]
  sequence_2 = seq2_seqs[[i]]

  hit_1 = matchPWM(pwm, sequence_1, with.score = TRUE)
  hit_2 = matchPWM(reverseComplement(pwm), sequence_2, with.score=TRUE)
  hit_3 = matchPWM(pwm, sequence_2, with.score = TRUE)
  hit_4 = matchPWM(reverseComplement(pwm), sequence_1, with.score=TRUE)

  if (max(mcols(hit_1)$score) > threshold & max(mcols(hit_2)$score) > threshold){
    result_convergent[i] = 1
  }

  if ((max(mcols(hit_1)$score) > threshold & max(mcols(hit_3)$score) > threshold) ||
      (max(mcols(hit_2)$score) > threshold & max(mcols(hit_4)$score) > threshold)){
    result_samedirection[i] = 1
  }

  if (max(mcols(hit_3)$score) > threshold & max(mcols(hit_4)$score) > threshold){
    result_divergent[i] = 1
  }

  if (max(mcols(hit_1)$score) > threshold ||
      max(mcols(hit_2)$score) > threshold ||
      max(mcols(hit_3)$score) > threshold ||
      max(mcols(hit_4)$score) > threshold){
        result_one[i] = 1
  }

  if ((max(mcols(hit_1)$score) > threshold & max(mcols(hit_3)$score) > threshold) ||
      (max(mcols(hit_2)$score) > threshold & max(mcols(hit_4)$score) > threshold) ||
      (max(mcols(hit_1)$score) > threshold & max(mcols(hit_2)$score) > threshold) ||
      (max(mcols(hit_3)$score) > threshold & max(mcols(hit_4)$score) > threshold)){
    result_two[i] = 1
  }
}

# Summarize statistics
rate_convergent = sum(result_convergent)/length(seq_data[,1])
rate_divergent = sum(result_divergent)/length(seq_data[,1])
rate_one = sum(result_one)/length(seq_data[,1])
rate_two = sum(result_two)/length(seq_data[,1])
rate_samedirection = sum(result_samedirection)/length(seq_data[,1])

result = c(rate_convergent, rate_divergent, rate_samedirection, rate_two, rate_one)
result = matrix(result, nrow=1)

# Extract the convergent pairs
index = as.logical(result_convergent)
seq_data = seq_data[index,c(1:6)]

# Rename the rows of the statistics data.frame as the input file name
rownames(result) = basename(opts$input)

write.table(x=result, file=paste("statistics_result",sep=""), append=T, quote = F, sep = "\t", row.names = T, col.names = F)
write.table(x=seq_data, file=paste(opts$input, ".with_motif", sep=""), append=F, quote = F, sep = "\t", row.names = F, col.names = F)
