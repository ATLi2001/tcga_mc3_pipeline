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

#start to loop through and run the summ stat analysis for each cancer
cancer_type <- unique(cancer_type)

for(i in 1:length(cancer_type)){
  #identify clinical data file 
  clinical <- grep(paste0("patient_",cancer_type[i]),pan12_manifest$filename)
  clinical <- paste0("pan12_clinical/",pan12_manifest$filename[clinical])
  
  print(cancer_type[i])
  
  #run summ_stats_multi_analysis
  system(paste("Rscript summ_stats_analysis_final.R",
               cancer_type[i],
               clinical,
               "summ_stats_pan12",
               paste0("/home/austintianli/summ_stats_analysis_final/",cancer_type[i])))
  
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

#start to loop through and run the summ stat analysis for each cancer
cancer_type <- unique(cancer_type)

for(i in 1:length(cancer_type)){
  #identify clinical data file 
  clinical <- grep(paste0("patient_",cancer_type[i]),non_pan12_manifest$filename)
  clinical <- paste0("non_pan12_clinical/",non_pan12_manifest$filename[clinical])
  
  print(cancer_type[i])
  
  #run summ_stats_multi_analysis
  system(paste("Rscript summ_stats_analysis_final.R",
               cancer_type[i],
               clinical,
               "summ_stats_non_pan12",
               paste0("/home/austintianli/summ_stats_analysis_final/",cancer_type[i])))
  
}

