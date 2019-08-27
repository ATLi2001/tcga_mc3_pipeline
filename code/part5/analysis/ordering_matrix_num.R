library(GenomicRanges)
library(dplyr)
# The annotation data
load("grch38_cosmic_census.rda")

#!/usr/bin/env Rscript

#input arguments should be
#1. cancer type
#2. input folder
#3. output folder
input_args <- commandArgs(trailingOnly = TRUE)

#test if there are right number of argument; if not, return an error
if(length(input_args) != 3){
  stop("Improper number of input arguments given", call.=FALSE)
}

print("Reading Input Arguments...")

#get input args
cancer_type <- input_args[1]
path <- input_args[2]
out_path <- input_args[3]


print("Loading ordering matrices...")

#depth file name strings
filelist <- list.files(path = path, pattern = ".ordering_matrix.csv")
#add on full path name
filelist <- paste(path,filelist, sep = "/")
#keep only those in the current cancer type
filelist <- filelist[grep(paste0("ssm_data_",cancer_type),filelist)]

#get list of the matrices
datalist <- lapply(filelist, FUN = read.csv, header=FALSE, 
                   stringsAsFactors = FALSE)

#get the patient ids from the input file names
patient_ids <- rep("",length(filelist))
for(i in 1:length(filelist)){
  #[[1]] as gregexpr returns a list
  st <- gregexpr("TCGA",filelist[i])[[1]] 
  #length of id is 12
  patient_ids[i] <- substr(filelist[i],st,st+11)
}


print("Loading SSM mapping files...")

#path to mapping files
if(grepl("non",path)){
  path_mapping <- "non_pan12_annot"
}else{
  path_mapping <- "pan12_annot"
}
#all mapping file names are patient_id.filtered.mapping.txt
filelist_mapping <- paste(patient_ids,"filtered.mapping.txt",sep = ".")
filelist_mapping <- paste(path_mapping,filelist_mapping,sep = "/")

#get mapping data
datalist_mapping <- lapply(filelist_mapping,FUN=read.table,header=TRUE,
                           stringsAsFactors=FALSE)

all_genes <- c()

print("Loading maf...")

#read in original maf data
maf <- read.delim("mc3.v0.2.8.PUBLIC.filtered.maf",
                  header = TRUE,
                  stringsAsFactors = FALSE,
                  comment.char = "#")

print("Done loading maf")

#combine the gene name,variant classification and HGVSp short 
#matches the mapping identification
maf$Mapping <- paste(maf$Hugo_Symbol, maf$HGVSp_Short, sep = "_")

#list of dataframes based by patient
patient_list <- list()
patient_list <- split(maf, f = maf$Tumor_Sample_Barcode)
#keep only the patients that are in the study
keep <- which(substr(names(patient_list),1,12) %in% patient_ids)
patient_list <- patient_list[keep]

#keep only relevant mutations
relevant_mutations <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")

#get original cnv files
path_cnv <- "phylowgs_input"

#need to trim down to only patients that ended up working
path_cnv_patients <- list.dirs(path = path_cnv,recursive = FALSE,full.names = FALSE)
path_cnv_patients_sub <- substr(path_cnv_patients,1,12)
matches <- which(path_cnv_patients_sub %in% patient_ids)
path_cnv_patients <- path_cnv_patients[matches]
#for some reason prad cancer has some paths with \r in them
#get rid of those
path_cnv_patients <- path_cnv_patients[!grepl("\r",path_cnv_patients)]

path_cnv <- paste(path_cnv,path_cnv_patients,"cnv_data.txt",sep = "/")


print("Creating CNV to SSM mapping files")

#clear any existing mapping files
path_cnv_mapping <- "cnv_data_mapping"
if(!dir.exists(path_cnv_mapping)){
  dir.create(path_cnv_mapping,showWarnings = FALSE)
}
unlink(paste0(path_cnv_mapping,"/*"))

#patients that will need to be deleted due to issues with cnv mapping
pat_del <- c()

#create mapping files
for(i in 1:length(path_cnv)){
  ensembl_gr = makeGRangesFromDataFrame(ensembl_df,
                                        seqnames.field = 'chromosome_name',
                                        start.field = 'start_position',
                                        end.field = 'end_position',
                                        keep.extra.columns=T)
  # Map physical_cnvs to genes
  #make sure path exists
  if(file.exists(path_cnv[i])){
    dt = read.delim(path_cnv[i],stringsAsFactors = F)
  }
  #if not, skip
  else{
    #delete later
    pat_del <- c(pat_del,i)  
    #write a blank file to not throw off indices later
    setwd(path_cnv_mapping)
    write.table("",paste0('cnv_data.',patient_ids[i],'.mapping.txt'),row.names=F,quote=F,sep='\t')
    setwd("..")
    next()
  }
  
  #if there are no cnvs, skip
  if(length(dt$physical_cnvs)==0){
    #delete later
    pat_del <- c(pat_del,i)  
    #write a nonsense file to not throw off indices later
    setwd(path_cnv_mapping)
    write.table(dt,paste0('cnv_data.',patient_ids[i],'.mapping.txt'),row.names=F,quote=F,sep='\t')
    setwd("..")
    next()
  }
  cnv = strsplit(dt$physical_cnvs,split = ';') 
  allgenes = NULL
  genes =sapply(cnv,function(x){
    
    item = strsplit(x,split = ',')
    
    y = lapply(item,function(x){
      gsub('[^=]+=','',x)
    })
    
    y = do.call(rbind.data.frame, y)
    colnames(y) = c("chr","start",'end','major','minor','prev')
    gr = makeGRangesFromDataFrame(y,seqnames.field = 'chr',start.field = 'start',end.field = 'end')
    mapped = findOverlaps(gr,ensembl_gr)
    hit = as.data.frame(ensembl_gr[subjectHits(mapped),]) %>% filter(hgnc_symbol!='')
    hit = unique(hit$hgnc_symbol)
    allgenes <<- c(allgenes,hit)
    return(paste(hit,collapse =';'))
  })
  
  message(sprintf("Number of Unique Genes mapped = %d",length(unique(allgenes))))
  # add genes column to the dt, the new column will have all the genes
  # in the physical_cnv regions
  dt$genes = genes
  #change directory so table writes into folder
  setwd(path_cnv_mapping)
  write.table(dt,paste0('cnv_data.',patient_ids[i],'.mapping.txt'),row.names=F,quote=F,sep='\t')
  setwd("..")
}

#path to mapped cnv files
path_cnv_mapping <- list.files(path_cnv_mapping,full.names = TRUE)

print("Starting combining process...")

#for each patient, edit their matrix
for(i in 1:length(datalist)){
  
  if(i %in% pat_del)
    next()
  
  print(paste("Starting patient",i,"of",length(datalist)))
  
  #edit the row and col names in datalist to the mutation names
  
  tmp <- datalist[[i]]
  curr_patient <- patient_list[[i]]
  
  #get the mutation names
  nm <- tmp[,1]
  nm <- nm[2:length(nm)]
  
  #get rid of mutation name columns
  tmp <- tmp[-1,]
  tmp <- tmp[,-1]
  
  #use mapping to convert to actual gene names
  for(mpf in 1:nrow(datalist_mapping[[i]])){
    if(substr(nm[mpf],1,1)!='c')
      nm[mpf] <- datalist_mapping[[i]]$Mutation[mpf]
  }
  
  
  #get rid of mutations that aren't relevant
  for(j in length(nm):1){
    if(substr(nm[j],1,1)!='c'){
      mut_ind <- which(curr_patient$Mapping==nm[j])
      if(!isTRUE(curr_patient$Variant_Classification[mut_ind] %in% relevant_mutations)){
        nm <- nm[-j]
        tmp <- tmp[-j,]
        tmp <- tmp[,-j]
      }
    }
  }
  
  #edit mutation names to only gene name 
  for(j in 1:length(nm)){
    if(substr(nm[j],1,1)!='c')
      #gene name is before the _
      nm[j] <- substr(nm[j],1,unlist(gregexpr(pattern = "[:_:]", nm[j]))[1]-1)
  }
  
  #read in cnv mapped file
  cnv_mapped <- read.delim(path_cnv_mapping[i],stringsAsFactors = FALSE)
  cnv_genes <- cnv_mapped$genes
  if(length(cnv_genes)==0 | is.na(cnv_genes)) #if cnv_genes is empty, go to next patient
    next()
  cnv_genes <- strsplit(cnv_genes, split = ";")
  
  empty_cnv <- c()
  
  if(length(dim(tmp))!=2) #if the ordering matrix is not a matrix??
    next()
  
  print("Expanding CNV's into genes...")
  
  #expand all cnvs to the genes they cover
  for(v in 1:length(cnv_genes)){
    #cnv_genes[[v]] has > 1 is normal
    #cnv_genes[[v]] has 1 is strange; is a skip anyways
    if(length(cnv_genes[[v]]) > 1){
      #convert the index to what it would be in the matrix
      #c0,c1, etc. start after all the ssms
      cnv_tmp_ind <- length(tmp)-(length(cnv_genes)-v)
      
      #dataframe before the current cnv
      pre_df <- tmp[1:(cnv_tmp_ind-1),1:(cnv_tmp_ind-1),drop=F]
      
      #column of data for the cnv
      cnv_tmp_col <- tmp[,cnv_tmp_ind]
      #repeat for each of the genes covered by the cnv
      pre_df_extend <- data.frame(matrix(rep(cnv_tmp_col,length(cnv_genes[[v]])),ncol = length(cnv_genes[[v]])),stringsAsFactors = FALSE)
      #add to the right of pre_df so get rid of rows below
      pre_df_extend <- pre_df_extend[1:(cnv_tmp_ind-1),]
      
      #if there would be stuff after the cnv
      if(cnv_tmp_ind<nrow(tmp)){
        #new matrix with expanded data; also data after the cnv
        tmp_expand <- cbind.data.frame(pre_df,pre_df_extend,tmp[1:(cnv_tmp_ind-1),(cnv_tmp_ind+1):ncol(tmp)])
        
        #row of data for the cnv
        cnv_tmp_row <- tmp[cnv_tmp_ind,]
        
        #repeated cnv data for each gene it covers
        row_extend <- cnv_tmp_row[rep(seq.int(1,nrow(cnv_tmp_row)),length(cnv_genes[[v]])),]
        #row_extend but before the cnv-cnv comparison
        row_extend1 <- row_extend[,1:(cnv_tmp_ind-1),drop=F]
        
        #have to match colnames for rbind
        names(row_extend1) <- names(tmp[1:(cnv_tmp_ind-1)])
        #add on the data below the cnv
        row_extend1 <- rbind(row_extend1,tmp[(cnv_tmp_ind+1):nrow(tmp),1:(cnv_tmp_ind-1),drop=F])
        
        #cnv-cnv comparison section of the rows
        #all 0's for genes covered in same cnv
        row_extend2 <- data.frame(matrix(0,nrow = length(cnv_genes[[v]]),ncol = length(cnv_genes[[v]])))
        #match the colnames to what would be above it
        colnames(row_extend2) <- colnames(tmp_expand)[cnv_tmp_ind:(cnv_tmp_ind+length(cnv_genes[[v]])-1)]
        #data after the current cnv compared to cnv
        after_cnv_col <- tmp[(cnv_tmp_ind+1):nrow(tmp),cnv_tmp_ind]
        row_extend2_below <- data.frame(matrix(rep(after_cnv_col,length(cnv_genes[[v]])),ncol = length(cnv_genes[[v]])))
        #match colnames for rbind
        colnames(row_extend2_below) <- colnames(row_extend2)
        row_extend2 <- rbind(row_extend2,row_extend2_below)
        
        #after cnv-after cnv comparison section
        row_extend3 <- as.data.frame(tmp[(cnv_tmp_ind+1):nrow(tmp),(cnv_tmp_ind+1):ncol(tmp)],stringsAsFactors=FALSE)
        #data of cnv compared to after cnv
        after_cnv_row <- as.data.frame(tmp[cnv_tmp_ind,(cnv_tmp_ind+1):ncol(tmp)],stringsAsFactors=FALSE)
        #cnv-after cnv comparison
        row_extend3_above <- as.data.frame(after_cnv_row[rep(seq.int(1,nrow(after_cnv_row)),length(cnv_genes[[v]])),],stringsAsFactors=FALSE)
        
        #match colnames for rbind
        colnames(row_extend3_above) <- colnames(row_extend3)
        row_extend3 <- rbind(row_extend3_above,row_extend3)
        
        #overall rows to extend the matrix by
        row_extend <- cbind.data.frame(row_extend1,row_extend2,row_extend3)
      }
      #no stuff after cnv
      else{
        #new matrix with expanded data; also data after the cnv
        tmp_expand <- cbind.data.frame(pre_df,pre_df_extend)
        
        #row of data for the cnv
        cnv_tmp_row <- tmp[cnv_tmp_ind,]
        
        #repeated cnv data for each gene it covers
        row_extend <- cnv_tmp_row[rep(seq.int(1,nrow(cnv_tmp_row)),length(cnv_genes[[v]])),]
        #row_extend but before the cnv-cnv comparison
        row_extend1 <- row_extend[,1:(cnv_tmp_ind-1),drop=F]
        
        #have to match colnames for rbind
        colnames(row_extend1) <- colnames(tmp)[1:(cnv_tmp_ind-1)]
        
        #cnv-cnv comparison section of the rows
        #all 0's for genes covered in same cnv
        row_extend2 <- data.frame(matrix(0,nrow = length(cnv_genes[[v]]),ncol = length(cnv_genes[[v]])))
        #match the colnames to what would be above it
        colnames(row_extend2) <- colnames(tmp_expand)[cnv_tmp_ind:(cnv_tmp_ind+length(cnv_genes[[v]])-1)]
        
        #overall rows to extend the matrix by
        row_extend <- cbind(row_extend1,row_extend2)
      }
      #updated matrix
      colnames(row_extend) <- colnames(tmp_expand)
      tmp_expand <- rbind.data.frame(tmp_expand,row_extend,stringsAsFactors = FALSE)
      tmp <- tmp_expand
    }
    #cnv_genes[[v]] has 0 is a skip too; deleted later
    if(length(cnv_genes[[v]])==0){
      cnv_tmp_ind <- length(tmp)-(length(cnv_genes)-v)
      empty_cnv <- c(empty_cnv,cnv_tmp_ind)
    }      
  }
  
  #convert all data to numeric
  tmp <- data.frame(sapply(tmp,function(x) as.numeric(as.character(x))))
  
  #delete columns of empty cnv data
  if(length(empty_cnv) > 0){
    tmp <- tmp[-empty_cnv,]
    tmp <- tmp[,-empty_cnv]
  }
  
  #look for repeats
  
  #combine ssm and cnv genes
  nm <- c(nm,unlist(cnv_genes))
  for(n in length(nm):1){
    if(substr(nm[n],1,1)=='c'){
      nm <- nm[-n]
    }
  }
  
  #get the duplicated gene names
  dup <- duplicated(nm)
  dup_values <- nm[dup]
  dup_values <- unique(dup_values)
  
  
  #if there are duplicated values
  if(length(dup_values) > 0){
    print("Combining duplicated genes...")
    #for each duplicated value
    for(k in 1:length(dup_values)){
      #get the rows of all genes matching the duplicated value
      tmp_ind <- which(nm == dup_values[k])
      
      #leave only the row that has the min average 
      #use unname to prevent extra information
      means <- unname(rowMeans(tmp, na.rm = TRUE)[tmp_ind])
      #use min to get earliest occurence
      tmp_ind_min <- which.min(means)
      
      #indices to delete
      tmp_ind_del <- tmp_ind[-tmp_ind_min]
      #order in decreasing to prevent skips
      tmp_ind_del <- tmp_ind_del[order(tmp_ind_del, decreasing = TRUE)]
      
      #delete the rows and columns
      for(l in 1:length(tmp_ind_del)){
        tmp <- tmp[-tmp_ind_del[l],]
        tmp <- tmp[,-tmp_ind_del[l]]
      }
      
      nm <- nm[-tmp_ind_del]
      
    }
  }
  
  #set the row and col names
  colnames(tmp) <- nm
  rownames(tmp) <- nm
  
  datalist[[i]] <- tmp
  
  #union the overall gene list
  all_genes <- union(nm, all_genes)
}

# #cut down to COSMIC gene list
# cosmic <- read.csv("~/Desktop/Austin/MSTC/Capstone Project/PhyloWGS/summ_stats_ordering_matrix/Census_allWed Aug  8 19_43_26 2018.csv",
#                    stringsAsFactors = FALSE,
#                    header = TRUE)
# cosmic_genes <- cosmic$Gene.Symbol
# final_genes <- intersect(all_genes, cosmic_genes)

#cut down to significant gene list
sig_genes <- read.table(paste0("sig_genes/",cancer_type,"_sig_genes.txt"),stringsAsFactors = FALSE)
sig_genes <- sig_genes$V1
final_sig_genes <- intersect(all_genes,sig_genes)

#create the overall ordering array
#make it 3 dimensional to hold extra info
#row x col is still genes x genes
#each element in first layer is the sum of all appearances
#each element in second layer is the number of appearances
#can collapse to 2 dimensional in the end by taking average
ordering_array <- array(NA, c(length(final_sig_genes),length(final_sig_genes),2))
ordering_array[,,2] <- 0
rownames(ordering_array) <- final_sig_genes
colnames(ordering_array) <- final_sig_genes

#ordering matrix will be array collapsed down
ordering_matrix <- ordering_array[,,1]

print("Creating ordering matrix...")

#loop through the entire array
for(r in 1:nrow(ordering_array)){
  for(c in 1:ncol(ordering_array)){
    #go through the patient to find any matches to rowname r and colname c
    for(l in 1:length(datalist)){
      if(is.element(final_sig_genes[r],rownames(datalist[[l]]))
         && is.element(final_sig_genes[c], colnames(datalist[[l]]))){
        
        #which row and col to add data from
        row_to_add <- which(rownames(datalist[[l]]) == final_sig_genes[r])
        col_to_add <- which(colnames(datalist[[l]]) == final_sig_genes[c])
        
        #put that value into the 3d array
        ordering_array[r,c,1] <- sum(ordering_array[r,c,1],datalist[[l]][row_to_add,col_to_add],
                                     na.rm = TRUE)
        #increase appearances
        ordering_array[r,c,2] <- ordering_array[r,c,2] + 1
      }#end if
    }#end for l
    
    #each element of the ordering matrix is sum/num_apps
    ordering_matrix[r,c] <- ordering_array[r,c,1]/ordering_array[r,c,2]
    #do weighting by cooccurence
    actual_prop <- ordering_array[r,c,2]/length(datalist)
    expected_prop <- ordering_array[r,r,2]/length(datalist) * ordering_array[c,c,2]/length(datalist)
    co_weight <- actual_prop - expected_prop
    if(co_weight<0){
      message(sprintf("Co weight < 0 at %d",c(r,c)))
      co_weight <- 0
    }
    #multiply by weight
    ordering_matrix[r,c] <- ordering_matrix[r,c]*co_weight
    
  }#end for c
}#end for r

#order based on unweighted rows
unweighted_rows <- rep(0, nrow(ordering_matrix))

for(m in 1:nrow(ordering_matrix)){
  unweighted_rows[m] <- mean(ordering_matrix[m,], na.rm = TRUE)
}

#overall formula
#each individual element is an unweighted average
#multiply each element by the cooccurence weight
#cooccurence weight = actual cooccurence prop - expected
#expected = product of proportion of appearances each gene has
#however, order the rows based on the unweighted averages
ordering_matrix <- ordering_matrix[order(unweighted_rows),]
ordering_matrix <- ordering_matrix[,order(unweighted_rows)]


#set overall output directory
if(!dir.exists(out_path)){
  dir.create(out_path,showWarnings = F)
}
setwd(out_path)

#write out the ordered genes in the heatmap
write.csv(ordering_matrix,file = paste0(cancer_type,"_ordering_matrix.csv"),
          quote = F,row.names = F,col.names = F)


print("Ordering matrix analysis complete!")
