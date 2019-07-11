# source("https://bioconductor.org/biocLite.R")
# biocLite("maftools")
library(maftools)
library(vcfR)
library(dplyr)


print("Reading maf")

maf = read.maf("mc3.v0.2.8.PUBLIC.filtered.maf",isTCGA = T,verbose = F)

print("Done reading maf")

# subset maf for patient
vcfs = read.delim("vcf_pan12.tsv",stringsAsFactors = F)

#because there are two vcfs per patient, go by patient
patients <- unique(vcfs$patient)

#directory for pan12
of = "phylowgs_input/"

#directory for ssm mapping files
if(!dir.exists('pan12_annot/')){
  dir.create('pan12_annot/',showWarnings = F)
}

######################### Part3
# The phylowgs parser will remove some snps from the vcf generated above, 
# so the mapping will not work. The code below remap annotate the mapping between 
# the ssm_id and the maf gene annotation
# Store the ssm_data file created by phylowgs somewhere
ssm_files = list.files(path=of,pattern="ssm_data.txt",full.names = T,recursive = T)

for(i in 1:length(patients)){
  pt = patients[i]
  
  #only run if file isn't already there
  if(!file.exists(paste0('pan12_annot/',pt,'.filtered.mapping.txt'))){
    #use trycatch to continue after errors
    tryCatch({
      print(i)
      print(pt)
      
      mafsub = subsetMaf(maf,tsb=pt,includeSyn = F)
      
      #on the cloud, mafsub is an maf class object rather than dataframe
      if(!is.data.frame(mafsub)){
        mafsub <- mafsub@data
      }
      
      f = grep(pt,ssm_files,value=T)
      if(length(f)==0){
        next()
      }
      ssm = read.delim(f,stringsAsFactors = F)
      
      idx = match(ssm$gene ,
                  paste(gsub("chr",'',mafsub$Chromosome),mafsub$Start_Position,sep = '_'))
      
      annot = mafsub[idx[!is.na(idx)],]
      if(nrow(annot)!=nrow(ssm)){
        stop()
      }
      tmp = cbind("ID"=paste0("s",(1:nrow(annot)-1)),"Mutation"=paste(annot$Hugo_Symbol,annot$HGVSp_Short,sep='_'))
      write.table(tmp,paste0('pan12_annot/',pt,'.filtered.mapping.txt'),row.names = F,quote=F,sep='\t' )
      
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}
