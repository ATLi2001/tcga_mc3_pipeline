#tsv file used to get the phylowgs trees
phylowgs_tsv <- read.delim("~/Desktop/Austin/MSTC/Capstone Project/Austin Pipeline MC3/Google Cloud/phylowgs_multi/pan12_multi.tsv",
                           header = FALSE,
                           stringsAsFactors = FALSE)

#pan12 metadata file to match patient to cancer type
pan12_meta <- read.delim("~/Desktop/Austin/MSTC/Capstone Project/Austin Pipeline MC3/TitanCNA/phylowgs_mc3-pan12-meta.pan12.paired_bam_hg38_chr.tsv",
                         stringsAsFactors = FALSE)

#old output is now the new input
phylowgs_tsv$V1 <- phylowgs_tsv$V3
phylowgs_tsv$V1 <- paste(phylowgs_tsv$V1, "trees.zip",sep="/")
phylowgs_tsv$V1[1] <- "--input INPUT"

#second column is the SSM name
phylowgs_tsv$V2[1] <- "SSM"
for(i in 2:nrow(phylowgs_tsv)){
  
  #patient sample barcode
  curr_bar <- substr(phylowgs_tsv$V1[i],48,63)
  
  #identify which cancer type this patient is of
  curr_cancer <- pan12_meta$project_short_name[which(pan12_meta$sample_barcode==curr_bar)]
  #format cancer type as lowercase without tcga prefix
  curr_cancer <- tolower(substring(curr_cancer,6))
  
  #all input file locations have a standard format
  #substr is the patient name
  phylowgs_tsv$V2[i] <- paste0("ssm_data_",curr_cancer,"_",curr_bar, ".txt")
}
#output column
phylowgs_tsv$V3 <- "gs://austin_pipeline_mc3/phylowgs/pan12/summ_stats/*"
phylowgs_tsv$V3[1] <- "--output OUTPUT_DIR"

#number of trees (k) to read column
# k = number of chains (n) * number of MCMC (s)
n <- as.numeric(phylowgs_tsv$V4[-1])
s <- as.numeric(phylowgs_tsv$V8[-1])
k <- n*s

phylowgs_tsv <- phylowgs_tsv[,-c(5:9)]
phylowgs_tsv$V4[1] <- "k"
phylowgs_tsv$V4[2:nrow(phylowgs_tsv)] <- k

write.table(phylowgs_tsv,file="pan12_tree_reader_multi.tsv",sep = "\t",row.names = F,col.names = F,quote = F)
