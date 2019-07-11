#!/usr/bin/env Rscript

#script to create most recent survival clinical file

path <- "non_pan12_clinical"

#follow up files
follow_up_files <- list.files(path = path,pattern = "follow_up",full.names = TRUE,recursive = TRUE)
#don't want parcel files
follow_up_files <- follow_up_files[!grepl("parcel",follow_up_files)]
follow_up <- lapply(follow_up_files, FUN = read.delim,header=TRUE,stringsAsFactors=FALSE)

#overall clinical files
clinical_files <- list.files(path = path,pattern = "patient",full.names = TRUE,recursive = TRUE)
#don't want parcel files
clinical_files <- clinical_files[!grepl("parcel",clinical_files)]
clinical <- lapply(clinical_files, FUN = read.delim,header=TRUE,stringsAsFactors=FALSE)

#for each follow_up file, look to update survival times/status
for(i in 1:length(follow_up)){
  #current follow up file
  curr_follow_up <- follow_up[[i]]
  
  #identify cancer type
  #start after last _
  st <- gregexpr("[:_:]",follow_up_files[i])[[1]]
  st <- st[length(st)]
  #stop before the .txt extension
  stp <- nchar(follow_up_files[i])-4
  cancer_type <- substr(follow_up_files[i],st+1,stp)
  
  #current clinical file
  curr_clinical <- grep(pattern=cancer_type,x=clinical_files)
  curr_clinical <- clinical[[curr_clinical]]
  
  #identify all patients with new follow up data
  follow_up_patients <- curr_follow_up$bcr_patient_barcode
  follow_up_patients <- follow_up_patients[grep("TCGA",follow_up_patients)]
  follow_up_patients <- unique(follow_up_patients)
  
  #for each patient with new follow up data, change original file
  for(j in 1:length(follow_up_patients)){
    #current patient
    curr_patient <- follow_up_patients[j]
    curr_patient_rows <- which(curr_follow_up$bcr_patient_barcode %in% curr_patient)
    #get most recent data
    curr_patient_rows <- curr_patient_rows[length(curr_patient_rows)]
    
    #current vital status
    curr_vital <- curr_follow_up$vital_status[curr_patient_rows]
    
    #current time
    curr_time <- 0
    
    #identify which row in the clinical data file
    curr_clinical_row <- which(curr_clinical$bcr_patient_barcode %in% curr_patient)
    
    #some clinical files have different names for columns
    if(!"last_contact_days_to" %in% colnames(curr_clinical)){
      col_name <- which(colnames(curr_clinical) == "days_to_last_followup")
      colnames(curr_clinical)[col_name] <- "last_contact_days_to"
    }
    if(!"last_contact_days_to" %in% colnames(curr_follow_up)){
      col_name <- which(colnames(curr_follow_up) == "days_to_last_followup")
      colnames(curr_follow_up)[col_name] <- "last_contact_days_to"
    }
    if(!"death_days_to" %in% colnames(curr_clinical)){
      col_name <- which(colnames(curr_clinical) == "days_to_death")
      colnames(curr_clinical)[col_name] <- "death_days_to"
    }
    if(!"death_days_to" %in% colnames(curr_follow_up)){
      col_name <- which(colnames(curr_follow_up) == "days_to_death")
      colnames(curr_follow_up)[col_name] <- "death_days_to"
    }
    #replace the data in the original file
    if(curr_vital == "Alive"){
      #get last follow up time
      curr_time <- curr_follow_up$last_contact_days_to[curr_patient_rows]
      #replace data from original clinical data file
      curr_clinical$vital_status[curr_clinical_row] <- curr_vital
      curr_clinical$last_contact_days_to[curr_clinical_row] <- curr_time
    }else{
      #get last follow up time
      curr_time <- curr_follow_up$death_days_to[curr_patient_rows]
      #replace data from original clinical data file
      curr_clinical$vital_status[curr_clinical_row] <- curr_vital
      curr_clinical$death_days_to[curr_clinical_row] <- curr_time
    }
    
  }
  
  #want to write new clinical files into master folder
  filename <- clinical_files[grep(pattern=cancer_type,x=clinical_files)]
  filename <- strsplit(filename,"/")[[1]]
  #write the clinical file again
  write.table(curr_clinical,file = paste(filename[1],filename[3],sep = "/"),
              row.names = F,col.names = T,quote = F,sep = "\t")
  
  filename <- follow_up_files[grep(pattern=cancer_type,x=follow_up_files)]
  filename <- strsplit(filename,"/")[[1]]
  #write the follow up file again
  write.table(curr_follow_up,file = paste(filename[1],filename[3],sep = "/"),
              row.names = F,col.names = T,quote = F,sep = "\t")
  
}


