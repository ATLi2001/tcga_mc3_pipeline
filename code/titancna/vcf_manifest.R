#meta data
pan12_meta <- read.delim("phylowgs_mc3-pan12-meta.pan12.paired_bam_hg38_chr.tsv",
                         sep = "\t",
                         stringsAsFactors = FALSE,
                         header = TRUE)
non_pan12_meta <- read.delim("phylowgs_mc3-non_pan12-non.pan12.paired_bam_hg38_chr.tsv",
                             sep = "\t",
                             stringsAsFactors = FALSE,
                             header = TRUE)

#all vcf file names
vcf_files <- read.delim("vcf_files.txt",
                        sep = "\t",
                        stringsAsFactors = FALSE,
                        header = FALSE,
                        col.names = "filename")

#case uuid is in first 36 char of the filename
vcf_files$case_uuid <- substr(vcf_files$filename,1,36)

#subset to pan12
vcf_files_pan12 <- vcf_files[vcf_files$case_uuid %in% pan12_meta$gdc_case_uuid,]

#get patient id
vcf_files_pan12$patient <- vcf_files_pan12$case_uuid
for(i in 1:nrow(vcf_files_pan12)){
  curr_uuid <- vcf_files_pan12$case_uuid[i]
  patient <- pan12_meta$case_barcode[which(pan12_meta$gdc_case_uuid == curr_uuid)]
  vcf_files_pan12$patient[i] <- patient
}

write.table(vcf_files_pan12,"vcf_pan12.tsv",row.names = F,quote=F,sep='\t')

#subset to non_pan12
vcf_files_non_pan12 <- vcf_files[vcf_files$case_uuid %in% non_pan12_meta$gdc_case_uuid,]

#get patient id
vcf_files_non_pan12$patient <- vcf_files_non_pan12$case_uuid
for(i in 1:nrow(vcf_files_non_pan12)){
  curr_uuid <- vcf_files_non_pan12$case_uuid[i]
  patient <- non_pan12_meta$case_barcode[which(non_pan12_meta$gdc_case_uuid == curr_uuid)]
  vcf_files_non_pan12$patient[i] <- patient
}

write.table(vcf_files_non_pan12,"vcf_non_pan12.tsv",row.names = F,quote=F,sep='\t')
