#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

inputdir = args[1]
myfiles <- list.files(inputdir, pattern = "\\.csv$")

annot <- read.csv(args[3], header=T, stringsAsFactors = F)

merged_df <- annot
for (file in myfiles){
  fullfile <- file.path(inputdir, file)
  this_df <- read.csv(fullfile, header=T, stringsAsFactors = F, row.names=1)
  this_df <- this_df[c('Gene', 'log2FoldChange', 'padj')]
  names(this_df) <- c('Gene', paste0(file, '.l2FC'), paste0(file, '.padj'))
  merged_df <- merge(merged_df, this_df, by.x = 'TIGR4.old', by.y='Gene', all.x=F, all.y=F, sort=F)
  merged_df <- unique(merged_df)
}


outfilename <- args[2]
write.csv(merged_df, outfilename)
