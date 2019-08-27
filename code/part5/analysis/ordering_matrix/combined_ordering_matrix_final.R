#path to overall ordering matrices
path <- "ordering_matrix_num_final"

#finished ordering matrix file name strings
filelist <- list.files(path = path, pattern = ".ordering_matrix.csv")

#need to identify the cancer types
cancer_type <- rep("",length(filelist))
for(i in 1:length(filelist)){
  #stop before the first _
  st <- gregexpr("[:_:]",filelist[i])[[1]]
  st <- st[1]
  cancer_type[i] <- substr(filelist[i],1,st-1)
}

#add on full path name
filelist <- paste(path,filelist, sep = "/")

#get list of the matrices
datalist <- lapply(filelist, FUN = read.csv, header=TRUE, 
                   stringsAsFactors = FALSE)

#get list of all unique genes across all cancer types
all_genes <- c()
for(i in 1:length(datalist)){
  all_genes <- c(all_genes,colnames(datalist[[i]]))
}
all_genes <- unique(all_genes)

#final matrix will be all genes x all cancer types
#NA = gene not in the cancer's driver genes
final_matrix <- matrix(NA, nrow = length(all_genes),ncol = length(cancer_type),
                       dimnames = list(all_genes,cancer_type))

#loop through matrix by col to so can check by cancer type
for(c in 1:ncol(final_matrix)){
  for(r in 1:nrow(final_matrix)){
    curr_gene <- rownames(final_matrix)[r]
    #if the gene is in the cancer type
    if(curr_gene %in% colnames(datalist[[c]])){
      #identify which col
      index <- which(colnames(datalist[[c]]) == curr_gene)
      #calculate the mean of the column and negate it
      #more negative value = more likely to occur early
      #more positive value = more likely to occur late
      value <- -mean(datalist[[c]][,index], na.rm = TRUE)
      final_matrix[r,c] <- value
    }
  }
}

write.csv(final_matrix,file="combined_ordering_matrix_final.csv",quote=F)
