#!/usr/bin/env Rscript

#manifest with clinical and recur data
pan12_manifest <- read.delim("pan12_clinical_manifest.txt",stringsAsFactors = FALSE)

#need to identify the cancer types
cancer_type <- rep("",nrow(pan12_manifest))
for(i in 1:nrow(pan12_manifest)){
  #start after last _
  st <- gregexpr("[:_:]",pan12_manifest$filename[i])[[1]]
  st <- st[length(st)]
  #stop before the .txt extension
  stp <- nchar(pan12_manifest$filename[i])-4
  cancer_type[i] <- substr(pan12_manifest$filename[i],st+1,stp)
}

#start to loop through and run the ordering matrix for each cancer
cancer_type <- unique(cancer_type)

for(i in 1:length(cancer_type)){
  print(cancer_type[i])
  
  #run ordering_matrix
  system(paste("Rscript ordering_matrix_num.R",
               cancer_type[i],
               "summ_stats_pan12",
               "/home/austintianli/ordering_matrix_num/"))
  
}

#manifest with clinical and recur data
non_pan12_manifest <- read.delim("non_pan12_clinical/non_pan12_clinical_manifest.txt",stringsAsFactors = FALSE)

#need to identify the cancer types
cancer_type <- rep("",nrow(non_pan12_manifest))
for(i in 1:nrow(non_pan12_manifest)){
  #start after last _
  st <- gregexpr("[:_:]",non_pan12_manifest$filename[i])[[1]]
  st <- st[length(st)]
  #stop before the .txt extension
  stp <- nchar(non_pan12_manifest$filename[i])-4
  cancer_type[i] <- substr(non_pan12_manifest$filename[i],st+1,stp)
}

#start to loop through and run the ordering matrix for each cancer
cancer_type <- unique(cancer_type)

for(i in 1:length(cancer_type)){
  print(cancer_type[i])

  #run ordering_matrix
  system(paste("Rscript ordering_matrix_num.R",
               cancer_type[i],
               "summ_stats_non_pan12",
               "/home/austintianli/ordering_matrix_num/"))

}
