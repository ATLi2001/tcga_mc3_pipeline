#manifest file names
non_pan12_manifest <- list.files("non_pan12_manifest",full.names = TRUE)
#manifest files
datalist <- lapply(non_pan12_manifest, FUN = read.delim,header=TRUE,stringsAsFactors=FALSE)

#identify rows of clinical and recurrence data files
out <- as.data.frame(matrix(nrow=0,ncol = ncol(datalist[[1]])))
for(i in 1:length(datalist)){
  #find for the file name strings
  clinical <- grep("clinical_patient",datalist[[i]]$filename)
  #many versions, keep version without nte
  recur <- grepl("follow_up",datalist[[i]]$filename) & !grepl("nte",datalist[[i]]$filename)
  recur <- which(recur == TRUE)
  #more than 1 version
  if(length(recur) > 1){
    #keep the most updated
    versions <- c()
    for(r in 1:length(recur)){
      #filename
      fname <- datalist[[i]]$filename[recur[r]]
      #extract version number
      versions <- c(versions,as.numeric(gsub("\\D", "", fname)))
    }
    #keep only largest, or most recent version
    recur <- recur[which.max(versions)]
  }
  
  out <- rbind(out,datalist[[i]][c(clinical,recur),])
}

write.table(out,file = "non_pan12_clinical_manifest.txt",col.names = T,row.names = F,sep = '\t',quote = F)
