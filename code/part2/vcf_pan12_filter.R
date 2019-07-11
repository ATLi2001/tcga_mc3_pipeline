#!/usr/bin/env Rscript

# source("https://bioconductor.org/biocLite.R")
# biocLite("maftools")
library(maftools)
print("maftools loaded")
library(vcfR)
print("vcfR loaded")
library(dplyr)
print("dplylr loaded")
###################### Part 2

print("Reading maf")

maf = read.maf("mc3.v0.2.8.PUBLIC.filtered.maf",isTCGA = T,verbose = F)

print("Done reading maf")

# subset maf for patient
vcfs = read.delim("vcf_pan12.tsv",stringsAsFactors = F)

print("Creating directory")

#directory for output
of = "tcga_phylowgs/"
if(!dir.exists(of)){
  dir.create(of,showWarnings = F)
}


print("Done creating directory")

#because there are two vcfs per patient, go by patient
patients <- unique(vcfs$patient)

print("Generating filtered vcfs")

for(i in 1:length(patients)){
  pt = patients[i]
  #only run if file isn't already there
  if(!file.exists(paste0(of,pt,'.filtered.vcf.gz'))){
    #use trycatch to continue after errors
    tryCatch({
      print(i)
      print(pt)
      
      mafsub = subsetMaf(maf,tsb=pt,includeSyn = F) # includeSyn=F to filter silent mutations
      
      #on the cloud, mafsub is an maf class object rather than dataframe
      if(!is.data.frame(mafsub)){
        mafsub <- mafsub@data
      }
      
      #find the rows in vcfs that correspond to patient id
      vcfs_rows <- which(vcfs$patient == pt)
      
      vcfFile = paste0("vcfs/",vcfs$filename[vcfs_rows]) #vcf file path
      
      #every patient has 2 vcfs
      vcf1 = read.vcfR(file=vcfFile[1],verbose = F)
      vcf2 = read.vcfR(file=vcfFile[2],verbose = F)
      #rbind2 combines them
      vcf <- rbind2(vcf1,vcf2)
      # get pass vcf
      vcf.pass = vcf[(vcf@fix[,"FILTER"])=="PASS",]
      
      # get intersection by chrom position and ref allele
      idx = match(apply(getFIX(vcf.pass)[,c("CHROM","POS","REF")],1,paste,collapse='-') ,
                 paste(mafsub$Chromosome,mafsub$Start_Position, mafsub$Reference_Allele,sep = '-'))
      
      vcf.pass = vcf.pass[!is.na(idx),]
      vcf.pass.annot = mafsub[idx[!is.na(idx)],]
      tmp = cbind("ID"=paste0("s",(1:nrow(vcf.pass.annot)-1)),"Mutation"=paste(vcf.pass.annot$Hugo_Symbol,vcf.pass.annot$Variant_Classification,vcf.pass.annot$HGVSp_Short,sep='_'))
      write.vcf(vcf.pass,paste0(of,pt,'.filtered.vcf.gz'),APPEND = F)
      #write.table(tmp,paste0('phylowgs/',pt,'.filtered.mapping.txt'),row.names = F,quote=F,sep='\t' )
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  }
}
# The next step is to use the filtered vcf file together with the cnv file to create phylowgs input
