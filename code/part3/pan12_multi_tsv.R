#list of patients with successful parsing
patient_list <- read.table("pan12_phylowgs.txt",stringsAsFactors = FALSE,header = FALSE)
patient_list <- patient_list$V1

#leave only TCGA Barcode
patient_list <- substring(patient_list,16,31)

#create a blank dataframe for the tsv
tsv <- as.data.frame(matrix("", nrow = 0, ncol = 9), stringsAsFactors = FALSE)
#column names for tsv
colnames(tsv) <- c("--input SSM", "--input CNV", "--output-recursive OUTPUT_DIR",
                   "n","r","I","B","s","i")

#fill out the tsv
for(i in 1:length(patient_list)){
  #--input SSM
  tsv_input <- paste("gs://austin_pipeline_mc3/phylowgs/pan12/INPUT/phylowgs_input",patient_list[i],
                     "ssm_data.txt", sep = "/")
  #--input CNA
  tsv_cna <- paste("gs://austin_pipeline_mc3/phylowgs/pan12/INPUT/phylowgs_input", patient_list[i],
                   "cnv_data.txt", sep = "/")
  #--output OUTPUT_DIR
  tsv_output <- paste("gs://austin_pipeline_mc3/phylowgs/pan12/OUTPUT",patient_list[i],
                      sep = "/")
  #note that default n=4,B=1000,s=2500,i=5000
  #r=-1 means let phylowgs choose a random seed
  #I=inf means keep all trees
  tsv_add <- c(tsv_input, tsv_cna, tsv_output,4,-1,"inf",1000,2500,5000)
  #make into a dataframe to then attach to the main output
  tsv_add <- as.data.frame(tsv_add)
  tsv_add <- t(tsv_add)
  colnames(tsv_add) <- c("--input SSM", "--input CNV", "--output-recursive OUTPUT_DIR",
                         "n","r","I","B","s","i")
  
  tsv <- rbind.data.frame(tsv, tsv_add, stringsAsFactors = FALSE)
}

write.table(tsv, file="pan12_multi.tsv",sep="\t",row.names=F,col.names = T,quote=F)
