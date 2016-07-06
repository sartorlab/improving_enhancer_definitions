library("GenomicRanges")

create_enhancer_list <- function(granges, min_size) {
  # extends each grange to the min_size (if not there already); if this would cause two neighboring ranges to overlap, it extends them only until they bookend
  
  # make sure everything overlapping already has been merged
  granges <- reduce(granges, min.gapwidth = 0)
  
  # resize the ranges that are smaller than min_size; keep the ranges that are bigger than min_size as is.
  was_resized <- rep(F, length(granges))
  was_resized[(end(granges) - start(granges)) < min_size] <- T
  
  granges[(end(granges) - start(granges)) < min_size] <- resize(granges[(end(granges) - start(granges)) < min_size], width = min_size, fix = "center")
  resized_starts <- start(granges)
  resized_ends <- end(granges)
  
  # find the cases where, after resizing, two ranges overlap
  overlaps <- findOverlaps(granges, granges)
  overlaps <- overlaps[(queryHits(overlaps)+1)==subjectHits(overlaps),]
  
  
  resect_first <- queryHits(overlaps)[was_resized[queryHits(overlaps)] & !was_resized[subjectHits(overlaps)]]
  resect_second <- queryHits(overlaps)[!was_resized[queryHits(overlaps)] & was_resized[subjectHits(overlaps)]]
  resect_both <- queryHits(overlaps)[was_resized[queryHits(overlaps)] & was_resized[subjectHits(overlaps)]]
  
  # resect
  end(granges)[resect_first] <- start(granges)[resect_first + 1]
  start(granges)[resect_second + 1] <- end(granges)[resect_second]
  overlap_amount <- end(granges)[resect_both] - start(granges)[resect_both+1]
  resection_amount <- ceiling(overlap_amount / 2)
  end(granges)[resect_both] <- end(granges)[resect_both] - resection_amount
  start(granges)[resect_both+1] <- start(granges)[resect_both+1] + resection_amount
  
  
  # NOTE: need to know the chromosome sizes so we can truncate at the chromosome starts/ends.  For now, just truncate at the chromosome starts...
  start(granges)[start(granges) <= 0] <- 1
  
  return(granges)
}


#####  TEST  #####

#df <- data.frame(chromosome=rep("chr1", 6), start = c(57, 61, 81, 97, 122, 127), end = c(60, 72, 95, 98, 124, 129))
#df.ranges <- makeGRangesFromDataFrame(df)
#df.ranges
#create_enhancer_list(df.ranges, min_size = 10)
#create_enhancer_list_2(df.ranges, min_size = 10)

#for(i in seq(2, 1000)) {
#  seqname <- paste("chr", as.character(i), sep = '')
#  df <- data.frame(chromosome=rep(seqname, 6), start = c(57, 61, 81, 97, 122, 127), end = c(60, 72, 95, 98, 124, 129))
#  df.ranges.new <- makeGRangesFromDataFrame(df)
#  df.ranges <- c(df.ranges, df.ranges.new)
#}
#create_enhancer_list(df.ranges, min_size = 10)
#create_enhancer_list_2(df.ranges, min_size = 10)

# df <- data.frame(chromosome=rep("chr1", 4), start = c(5, 14, 20, 32), end = c(7, 17, 25, 40))
# df.ranges <- makeGRangesFromDataFrame(df)
# df.ranges
# create_enhancer_list(df.ranges, min_size = 1)
# 
# 
# df <- data.frame(chromosome=rep("chr1", 4), start = c(5, 14, 20, 32), end = c(7, 17, 25, 40))
# df.ranges <- makeGRangesFromDataFrame(df)
# df.ranges
# create_enhancer_list(df.ranges, min_size = 5)
# 
# 
# df <- data.frame(chromosome=rep("chr1", 4), start = c(5, 14, 20, 32), end = c(7, 17, 25, 40))
# df.ranges <- makeGRangesFromDataFrame(df)
# df.ranges
# create_enhancer_list(df.ranges, min_size = 2000)
# 
# 
# df <- data.frame(chromosome=rep("chr1", 4), start = c(5, 6, 20, 32), end = c(6, 20, 25, 40))
# df.ranges <- makeGRangesFromDataFrame(df)
# df.ranges
# create_enhancer_list(df.ranges, min_size = 10)

# df <- data.frame(chromosome=rep("chr1", 4), start = c(5, 8, 20, 32), end = c(6, 20, 25, 40))
# df.ranges <- makeGRangesFromDataFrame(df)
# df.ranges
# create_enhancer_list(df.ranges, min_size = 10)
