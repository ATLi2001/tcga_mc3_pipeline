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

#blue to yellow to red
col_heatmap <- c("#295A9E","#FFF1CF","#A80505")

#custom heatmap function
.heatmap <- function(expres){
  library("RColorBrewer")
  library("gplots")
  expres = as.matrix(expres)
  temp = c(expres)
  n <- length(temp) - sqrt(length(temp))
  col_byr = colorRampPalette(colors = col_heatmap)(n=n)
  vv = seq(min(temp, na.rm = T), max(temp, na.rm = T), length.out = n+1)
  heatmap.2(expres, col = col_byr, scale="none",
            key = TRUE, symm=F,symkey=F,symbreaks=T,density.info="none",
            trace="none",cexRow=1,dendrogram=c("none"), cexCol=1,
            breaks=vv,Colv=F,Rowv=F,margins=c(6,8),lhei=c(1,4))
  
  dev.off()
}

out_path <- "ordering_matrix_all_results_final"

#set overall output directory
if(!dir.exists(out_path)){
  dir.create(out_path,showWarnings = F)
}
setwd(out_path)

for(i in 1:length(datalist)){
  rownames(datalist[[i]]) <- colnames(datalist[[i]])
  
  #write out the ordered genes in the heatmap
  write.table(colnames(datalist[[i]]),file = paste0(cancer_type[i],"_ordered_genes",".txt"),
              sep = "\n",quote = F,row.names = F,col.names = F)
  
  #save to jpeg file
  jpeg(filename = paste0(cancer_type[i],"_heatmap",".jpeg"),
       width = 1000,height = 1000,res = 100)
  .heatmap(datalist[[i]])
}


