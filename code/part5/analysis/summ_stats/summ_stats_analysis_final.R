#!/usr/bin/env Rscript

#input arguments should be
#1. cancer type
#2. clinical data file
#3. input folder
#4. output folder
#2 should be nationwidechildrens.org.*.txt files
input_args <- commandArgs(trailingOnly = TRUE)

#test if there are right number of argument; if not, return an error
if(length(input_args) != 4){
  stop("Improper number of input arguments given", call.=FALSE)
}

print("Reading Input Arguments...")

#get input args
cancer_type <- input_args[1]
clinical_path <- input_args[2]
path <- input_args[3]
out_path <- input_args[4]

print("Loading summary statistics...")

#depth file name strings
depth_filelist <- list.files(path = path,pattern = ".depth.txt")
#add on full path name
depth_filelist <- paste(path,depth_filelist, sep = "/")
#keep only those in the current cancer type
depth_filelist <- depth_filelist[grep(paste0("ssm_data_",cancer_type),depth_filelist)]

#get list of the depths
depth_datalist <- lapply(depth_filelist, FUN = read.table, header=FALSE)

#vector for the averages of the depths of the trees of each patient
depth_avg <- rep(0, length(depth_datalist))
for(i in 1:length(depth_datalist)){
  #because the list has dataframes use [,1] to get actual numbers
  depth_avg[i] <- mean(depth_datalist[[i]][,1])
}

#########################################################################

#max_muts file name strings
max_muts_filelist <- list.files(path = path,pattern = ".maxmuts.txt")
#add on full path name
max_muts_filelist <- paste(path,max_muts_filelist, sep = "/")
#keep only those in the current cancer type
max_muts_filelist <- max_muts_filelist[grep(paste0("ssm_data_",cancer_type),max_muts_filelist)]

#get list of the max_muts
max_muts_datalist <- lapply(max_muts_filelist, FUN = read.table, header=FALSE)

#vector for the averages of the max_muts of the trees of each patient
max_muts_avg <- rep(0, length(max_muts_datalist))
for(i in 1:length(max_muts_datalist)){
  #because the list has dataframes use [,1] to get actual numbers
  max_muts_avg[i] <- mean(max_muts_datalist[[i]][,1])
}

#########################################################################

#num_branches file name strings
num_branches_filelist <- list.files(path = path,pattern = ".numbranches.txt")
#add on full path name
num_branches_filelist <- paste(path,num_branches_filelist, sep = "/")
#keep only those in the current cancer type
num_branches_filelist <- num_branches_filelist[grep(paste0("ssm_data_",cancer_type),num_branches_filelist)]

#get list of the num_branches
num_branches_datalist <- lapply(num_branches_filelist, FUN = read.table, header=FALSE)

#vector for the averages of the num_branches of the trees of each patient
num_branches_avg <- rep(0, length(num_branches_datalist))
for(i in 1:length(num_branches_datalist)){
  #because the list has dataframes use [,1] to get actual numbers
  num_branches_avg[i] <- mean(num_branches_datalist[[i]][,1])
}

#########################################################################

#num_nodes file name strings
num_nodes_filelist <- list.files(path = path,pattern = ".numnodes.txt")
#add on full path name
num_nodes_filelist <- paste(path,num_nodes_filelist, sep = "/")
#keep only those in the current cancer type
num_nodes_filelist <- num_nodes_filelist[grep(paste0("ssm_data_",cancer_type),num_nodes_filelist)]

#get list of the num_nodes
num_nodes_datalist <- lapply(num_nodes_filelist, FUN = read.table, header=FALSE)

#vector for the averages of the num_nodes of the trees of each patient
num_nodes_avg <- rep(0, length(num_nodes_datalist))
for(i in 1:length(num_nodes_datalist)){
  #because the list has dataframes use [,1] to get actual numbers
  num_nodes_avg[i] <- mean(num_nodes_datalist[[i]][,1])
}

#########################################################################

#total_muts file name strings
total_muts_filelist <- list.files(path = path,pattern = ".totalmuts.txt")
#add on full path name
total_muts_filelist <- paste(path,total_muts_filelist, sep = "/")
#keep only those in the current cancer type
total_muts_filelist <- total_muts_filelist[grep(paste0("ssm_data_",cancer_type),total_muts_filelist)]

#get list of the total_muts
total_muts_datalist <- lapply(total_muts_filelist, FUN = read.table, header=FALSE)

#vector for the averages of the total_muts of the trees of each patient
total_muts_avg <- rep(0, length(total_muts_datalist))
for(i in 1:length(total_muts_datalist)){
  #because the list has dataframes use [,1] to get actual numbers
  total_muts_avg[i] <- mean(total_muts_datalist[[i]][,1])
}

#########################################################################

#trunk_prop file name strings
trunk_prop_filelist <- list.files(path = path,pattern = ".trunkprop.txt")
#add on full path name
trunk_prop_filelist <- paste(path,trunk_prop_filelist, sep = "/")
#keep only those in the current cancer type
trunk_prop_filelist <- trunk_prop_filelist[grep(paste0("ssm_data_",cancer_type),trunk_prop_filelist)]

#get list of the trunk_prop
trunk_prop_datalist <- lapply(trunk_prop_filelist, FUN = read.table, header=FALSE)

#vector for the averages of the trunk_prop of the trees of each patient
trunk_prop_avg <- rep(0, length(trunk_prop_datalist))
for(i in 1:length(trunk_prop_datalist)){
  #because the list has dataframes use [,1] to get actual numbers
  trunk_prop_avg[i] <- mean(trunk_prop_datalist[[i]][,1])
}

#########################################################################

#ce file name strings
ce_filelist <- list.files(path = path,pattern = ".CE.txt")
#add on full path name
ce_filelist <- paste(path,ce_filelist, sep = "/")
#keep only those in the current cancer type
ce_filelist <- ce_filelist[grep(paste0("ssm_data_",cancer_type),ce_filelist)]

#get list of the ce
ce_datalist <- lapply(ce_filelist, FUN = read.table, header=FALSE)

#vector for the averages of the ce of the trees of each patient
ce_avg <- rep(0, length(ce_datalist))
for(i in 1:length(ce_datalist)){
  #because the list has dataframes use [,1] to get actual numbers
  ce_avg[i] <- mean(ce_datalist[[i]][,1], na.rm = TRUE)
}

#########################################################################

#trunk_prop_cnv file name strings
trunk_prop_cnv_filelist <- list.files(path = path,pattern = ".trunkprop_cnv.txt")
#add on full path name
trunk_prop_cnv_filelist <- paste(path,trunk_prop_cnv_filelist, sep = "/")
#keep only those in the current cancer type
trunk_prop_cnv_filelist <- trunk_prop_cnv_filelist[grep(paste0("ssm_data_",cancer_type),trunk_prop_cnv_filelist)]

#get list of the trunk_prop_cnv
trunk_prop_cnv_datalist <- lapply(trunk_prop_cnv_filelist, FUN = read.table, header=FALSE)

#vector for the averages of the trunk_prop_cnv of the trees of each patient
trunk_prop_cnv_avg <- rep(0, length(trunk_prop_cnv_datalist))
for(i in 1:length(trunk_prop_cnv_datalist)){
  #because the list has dataframes use [,1] to get actual numbers
  trunk_prop_cnv_avg[i] <- mean(trunk_prop_cnv_datalist[[i]][,1])
}

print("Done loading summary statistics!")

#########################################################################

#all summary statistic data into one dataframe
summ_stats_data <- cbind.data.frame(depth_avg,max_muts_avg,num_branches_avg,num_nodes_avg,total_muts_avg,trunk_prop_avg,ce_avg,trunk_prop_cnv_avg)

#get the patient ids from the input file names
patient_ids <- rep("",length(depth_filelist))
for(i in 1:length(depth_filelist)){
  #[[1]] as gregexpr returns a list
  st <- gregexpr("TCGA",depth_filelist[i])[[1]] 
  #length of id is 12
  patient_ids[i] <- substr(depth_filelist[i],st,st+11)
}

#add patient_ids as a column
summ_stats_data <- cbind.data.frame(patient_ids, summ_stats_data)

print("Loading clinical data...")

#read in clinical data
clinical_data <- read.delim(clinical_path,
                            stringsAsFactors = FALSE,
                            header = TRUE)
#get rid of first two rows
clinical_data <- clinical_data[-c(1,2),]

#alphabetize
clinical_data <- clinical_data[order(clinical_data$bcr_patient_barcode),]
summ_stats_data <- summ_stats_data[order(summ_stats_data$patient_ids),]

#keep only patients that have both summary stats and clinical data
matches_index_summ_stats <- which(patient_ids %in% clinical_data$bcr_patient_barcode)
matches_index_clinical <- which(clinical_data$bcr_patient_barcode %in% patient_ids)
summ_stats_data <- summ_stats_data[matches_index_summ_stats,]
clinical_data <- clinical_data[matches_index_clinical,]

#########################################################################

#set up colors of yellow and blue
col_yb <- c("#FFC400", "#0033A0")

#corresponding summary statistic names
ss_names <- c("Depth", "Maximum Mutations", "Number of Branches", "Number of Nodes",
              "Total Mutations", "Trunk Proportion", "Clonal Expansion Index",
              "Trunk Proportion CNV")

print("Conducting survival analysis...")

#association with overall survival time

#set overall output directory
if(!dir.exists(out_path)){
  dir.create(out_path,showWarnings = F, recursive = T)
}
setwd(out_path)

#save overall summary statistic data
write.csv(summ_stats_data,file=paste0(cancer_type,"_summ_stats_data.csv"),quote = F,row.names=F,col.names=T)


#create subdir for survival analysis
surv_dir <- "survival_plots/"
if(!dir.exists(surv_dir)){
  dir.create(surv_dir,showWarnings = F)
}
#move to survival directory
setwd(surv_dir)

#load the survival package
library(survival)

tryCatch({
  #event (alive/dead = 1/2)
  data_event <- clinical_data$vital_status
  
  #get rid of na and not available stage data
  not <- which(data_event == "[Not Available]" | data_event == "[Discrepancy]")
  if(length(not) > 0){
    data_event <- data_event[-not]
  }
  data_event <- as.numeric(as.factor(data_event))
  
  #dataframe of clinical data with survival
  if(length(not) > 0)
    clinical_data_survival <- clinical_data[-not,]
  else
    clinical_data_survival <- clinical_data
  
  #followup time in days
  data_time <- rep(NA,length(data_event))
  for(i in 1:length(data_time)){
    #alive
    if(data_event[i] == 1)
      data_time[i] = clinical_data_survival$last_contact_days_to[i]
    #dead
    else
      data_time[i] = clinical_data_survival$death_days_to[i]
  }
  data_time <- as.numeric(data_time)
  
  #do an analysis for each summary stat
  #start at 2 to avoid the patient id column
  #quartiles
  q_dir <- "quartile"
  if(!dir.exists(q_dir)){
    dir.create(q_dir,showWarnings = F)
  }
  setwd(q_dir)
  
  #save logrank p, cox p, hr, hr lowerbound, hr upperbound into one dataframe
  .lp <- c()
  .cp <- c()
  .hr <- c()
  .hrl <- c()
  .hrh <- c()
  
  for(j in 2:ncol(summ_stats_data)){
    #current statistic
    curr_stat <- as.vector(summ_stats_data[,j])
    #get rid of patients without survival data
    if(length(not) > 0)
      curr_stat <- curr_stat[-not]
    
    #cutoff for low vs high
    cutoffl <- quantile(curr_stat, .25, na.rm = TRUE)
    cutoffh <- quantile(curr_stat, .75, na.rm = TRUE)
    # NA = not used, 1 = lower expression, 2 = higher expression
    teller <- rep(NA, length(curr_stat))
    teller[curr_stat <= cutoffl] <- 1
    teller[curr_stat >= cutoffh] <- 2
    
    #survival fitted data based on low vs high expression
    survplot <- survfit(Surv(data_time, data_event) ~ as.factor(teller))
    #calculate the logrank based upon the dichotomized variable
    logrank <- survdiff(Surv(data_time, data_event) ~ as.factor(teller))
    logrank_p = 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
    .lp <- c(.lp, logrank_p)
    #calculate the cox as a result of the continuous variable
    cox <- coxph(Surv(data_time, data_event) ~ curr_stat)
    cox_p = summary(cox)$coefficients[1,5]
    .cp <- c(.cp, cox_p)
    #hazard ratio interval
    hr = paste(round(summary(cox)$conf.int[1,1], digits=2), " (", round(summary(cox)$conf.int[1,3], digits=2), ", ", round(summary(cox)$conf.int[1,4], digits=2), ")", sep="")
    .hr <- c(.hr, summary(cox)$conf.int[1,1])
    .hrl <- c(.hrl, summary(cox)$conf.int[1,3])
    .hrh <- c(.hrh, summary(cox)$conf.int[1,4])
    
    #save to jpeg file
    jpeg(filename = paste0(substr(colnames(summ_stats_data)[j],1,nchar(colnames(summ_stats_data)[j])-4),".jpeg"),
         width = 1000,height = 1000,res = 100)
    #make margins for x, y axis larger to fit labels
    par(mar = c(8,8,6,4), mgp = c(3,1.5,0))
    
    main_str <- paste0(toupper(cancer_type),"\n",ss_names[j-1])
    plot(survplot, conf.int = FALSE, mark.time = TRUE, col = col_yb,
         main = main_str, cex.main = 2.5,xlab = "", ylab = "",
         #line width of survival curves thicker
         lwd = 4, 
         #text of axis labels larger
         cex.lab = 2.5, cex.axis = 2.5, 
         #text of all axis labels (not titles) horizontal
         las=1)
    #set titles of y and x axis
    #use line parameter to distance from axis labels
    title(ylab = "Survival Probability",line = 6, cex.lab = 2.5)
    title(xlab = "Days",line = 4, cex.lab = 2.5)
    #make the overall plot box thicker
    box(lwd = 3)
    #make tick marks thicker but don't redo labels
    axis(side = 1, lwd.ticks = 3, labels = FALSE)
    axis(side = 2, lwd.ticks = 3, labels = FALSE)
    #legend in the top right with no box around
    legend('topright', c("Low", "High"),
           col = col_yb, lty = 1, lwd = 4,cex = 2, bty = 'n')
    
    #create string for logrank p value to put in plot box
    logrank_p_str <- ""
    if(logrank_p < .001)
      logrank_p_str <- "P < .001"
    else
      logrank_p_str <- paste0("P = ", round(logrank_p, digits = 3))
    
    #add text inside plot at bottom left
    text(x=0, y=.08, pos = 4, cex = 2, offset = 2,
         labels = logrank_p_str)
    dev.off() #plot is now closed
  }
  
  #save p-value stuff 
  qnum_df <- data.frame(.lp, .cp, .hr, .hrl, .hrh, row.names = colnames(summ_stats_data)[2:ncol(summ_stats_data)])
  write.csv(qnum_df, file="survival_quartile_p.csv",quote = F,row.names = T)

  #reset directory to survival
  setwd("..")
  
  #median split
  m_dir <- "median"
  if(!dir.exists(m_dir)){
    dir.create(m_dir,showWarnings = F)
  }
  setwd(m_dir)
  
  #save logrank p, cox p, hr, hr lowerbound, hr upperbound into one dataframe
  .lp <- c()
  .cp <- c()
  .hr <- c()
  .hrl <- c()
  .hrh <- c()
  
  for(j in 2:ncol(summ_stats_data)){
    #current statistic
    curr_stat <- as.vector(summ_stats_data[,j])
    #get rid of patients without survival data
    if(length(not) > 0)
      curr_stat <- curr_stat[-not]
    
    #cutoff for low vs high
    cutoff <- quantile(curr_stat, .5, na.rm = TRUE)
    
    #1 = lower expression, 2 = higher expression
    teller <- rep(1, length(curr_stat))
    teller[curr_stat >= cutoff] <- 2
    
    #survival fitted data based on low vs high expression
    survplot <- survfit(Surv(data_time, data_event) ~ as.factor(teller))
    #calculate the logrank based upon the dichotomized variable
    logrank <- survdiff(Surv(data_time, data_event) ~ as.factor(teller))
    logrank_p = 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
    .lp <- c(.lp, logrank_p)
    #calculate the cox as a result of the continuous variable
    cox <- coxph(Surv(data_time, data_event) ~ curr_stat)
    cox_p = summary(cox)$coefficients[1,5]
    .cp <- c(.cp, cox_p)
    #hazard ratio interval
    hr = paste(round(summary(cox)$conf.int[1,1], digits=2), " (", round(summary(cox)$conf.int[1,3], digits=2), ", ", round(summary(cox)$conf.int[1,4], digits=2), ")", sep="")
    .hr <- c(.hr, summary(cox)$conf.int[1,1])
    .hrl <- c(.hrl, summary(cox)$conf.int[1,3])
    .hrh <- c(.hrh, summary(cox)$conf.int[1,4])
    
    #save to jpeg file
    jpeg(filename = paste0(substr(colnames(summ_stats_data)[j],1,nchar(colnames(summ_stats_data)[j])-4),".jpeg"),
         width = 1000,height = 1000,res = 100)
    #make margins for x, y axis larger to fit labels
    par(mar = c(8,8,6,4), mgp = c(3,1.5,0))
    
    main_str <- paste0(toupper(cancer_type),"\n",ss_names[j-1])
    plot(survplot, conf.int = FALSE, mark.time = TRUE, col = col_yb,
         main = main_str, cex.main = 2.5,xlab = "", ylab = "",
         #line width of survival curves thicker
         lwd = 4, 
         #text of axis labels larger
         cex.lab = 2.5, cex.axis = 2.5, 
         #text of all axis labels (not titles) horizontal
         las=1)
    #set titles of y and x axis
    #use line parameter to distance from axis labels
    title(ylab = "Survival Probability",line = 6, cex.lab = 2.5)
    title(xlab = "Days",line = 4, cex.lab = 2.5)
    #make the overall plot box thicker
    box(lwd = 3)
    #make tick marks thicker but don't redo labels
    axis(side = 1, lwd.ticks = 3, labels = FALSE)
    axis(side = 2, lwd.ticks = 3, labels = FALSE)
    #legend in the top right with no box around
    legend('topright', c("Low", "High"),
           col = col_yb, lty = 1, lwd = 4,cex = 2, bty = 'n')
    
    #create string for logrank p value to put in plot box
    logrank_p_str <- ""
    if(logrank_p < .001)
      logrank_p_str <- "P < .001"
    else
      logrank_p_str <- paste0("P = ", round(logrank_p, digits = 3))
    
    #add text inside plot at bottom left
    text(x=0, y=.08, pos = 4, cex = 2, offset = 2,
         labels = logrank_p_str)
    dev.off() #plot is now closed
  }
  
  mnum_df <- data.frame(.lp, .cp, .hr, .hrl, .hrh,row.names = colnames(summ_stats_data)[2:ncol(summ_stats_data)])
  write.csv(mnum_df, file="survival_median_p.csv",quote = F,row.names = T)
  
  
},error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

#reset to outpath
setwd(out_path)

#########################################################################

print("Conducting gender analysis...")

#gender

#create directory for gender plots
gender_dir <- "gender_plots/"
if(!dir.exists(gender_dir)){
  dir.create(gender_dir,showWarnings = F)
}
setwd(gender_dir)

#gender is the factor
gender <- clinical_data$gender
gender <- as.factor(gender)

tryCatch({
  
  #save t test p and wilcox p values
  .tp <- c()
  .wp <- c()
  
  #get a graph for each summary stat
  for(n in 2:ncol(summ_stats_data)){
    #current statistic
    curr_stat <- as.vector(summ_stats_data[,n])
    
    #need to use log scale for max and total mutations
    if(grepl("mutation",ss_names[n-1],ignore.case = TRUE))
      curr_stat <- log10(curr_stat)
    
    #significance tests
    t_test_p <- t.test(curr_stat ~ gender)$p.value
    wilcox_p <- wilcox.test(curr_stat ~ gender)$p.value
    
    .tp <- c(.tp,t_test_p)
    .wp <- c(.wp,wilcox_p)
    
    #save to jpeg file
    jpeg(filename = paste0(substr(colnames(summ_stats_data)[n],1,nchar(colnames(summ_stats_data)[n])-4),".jpeg"),
         width = 1000,height = 1000,res = 100)
    
    #create the graph
    par(mar = c(5,12,8,2))
    boxplot(curr_stat ~ gender,
            #no outlier points
            range = 0,
            ylab = "",
            xlab = toupper(cancer_type),
            main = "",
            las = 1,cex.lab = 2, cex.axis = 2,
            #make boxplot lines thicker
            boxlwd = 4,medlwd=6,whisklwd=3,staplelwd=3)
    #make overall plot box thicker
    box(lwd = 3)
    #make tick marks thicker but don't redo labels
    axis(side = 1, tick=FALSE, labels = FALSE)
    axis(side = 2, lwd.ticks = 3, labels = FALSE)
    
    #y label depends upon if there was log or not
    ylab_str <- ""
    if(grepl("mutation",ss_names[n-1],ignore.case = TRUE))
      ylab_str <- paste0(ss_names[n-1],"\n","Log10(mutation number)")
    else
      ylab_str <- ss_names[n-1]
    
    title(ylab = ylab_str,line = 6, cex.lab = 2.5)
    
    #create string for wilcox p value to put in plot box
    wilcox_p_str <- ""
    if(wilcox_p < .001)
      wilcox_p_str <- "P < .001"
    else
      wilcox_p_str <- paste0("P = ", round(wilcox_p, digits = 3))
    
    #add in Wilcoxon p value inside plot box
    text(x=.35, y=min(curr_stat, na.rm = T), pos = 4, cex = 1.5, offset = 2,
         labels = wilcox_p_str)
    
    #number of occurences for each gender
    occurences <- table(gender)
    
    #add in the specific data points
    points(x = jitter(rep(1, occurences[[1]]), amount = .05),
           curr_stat[gender == "FEMALE"],
           pch = 15, col = col_yb[1], cex = 1.5)
    points(x = jitter(rep(2, occurences[[2]]), amount = .05),
           curr_stat[gender == "MALE"],
           pch = 15, col = col_yb[2], cex = 1.5)
    
    dev.off() #close boxplot
  }
  
  gnum_df <- data.frame(.tp, .wp, row.names = colnames(summ_stats_data)[2:ncol(summ_stats_data)])
  write.csv(gnum_df,file="gender_p.csv",quote=F,row.names = T)
  
},error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

#reset to outpath
setwd(out_path)

#########################################################################

print("Conducting smoking analysis...")

#smoking anova

#smoking indicator
smoke <- clinical_data$tobacco_smoking_history_indicator

#if there is smoking data
if(!all(is.na(smoke))){
  
  #create directory for anova plots
  smoke_dir <- "smoking_plots/"
  if(!dir.exists(smoke_dir)){
    dir.create(smoke_dir,showWarnings = F)
  }
  setwd(smoke_dir)
  
  #get rid of na and not available stage data
  not <- which(smoke == "[Not Available]" | smoke == "[Unknown]" | smoke == "[Discrepancy]")
  smoke <- smoke[-not]
  
  #1 -> never smoked
  smoke[smoke != "1"] <- "Smoking"
  smoke[smoke == "1"] <- "Non-Smoking"
  
  
  smoke <- as.factor(smoke)
  
  #save t test p and wilcox p values
  .tp <- c()
  .wp <- c()
  
  #do an analysis for each summ stat
  #start at 2 to avoid the patient id column
  for(k in 2:ncol(summ_stats_data)){
    
    #current statistic
    curr_stat <- as.vector(summ_stats_data[,k])
    #get rid of patients without data
    curr_stat <- curr_stat[-not]
    
    #need to use log scale for max and total mutations
    if(grepl("mutation",ss_names[k-1],ignore.case = TRUE))
      curr_stat <- log10(curr_stat)
    
    #significance tests
    t_test_p <- t.test(curr_stat ~ smoke)$p.value
    wilcox_p <- wilcox.test(curr_stat ~ smoke)$p.value
    
    .tp <- c(.tp,t_test_p)
    .wp <- c(.wp,wilcox_p)
    
    #save to jpeg file
    jpeg(filename = paste0(substr(colnames(summ_stats_data)[k],1,nchar(colnames(summ_stats_data)[k])-4),".jpeg"),
         width = 1000,height = 1000,res = 100)
    
    #create the graph
    par(mar = c(5,12,8,2))
    boxplot(curr_stat ~ smoke,
            #no outlier points
            range = 0,
            names = c("Non-Smoking","Smoking"),
            ylab = "",
            xlab = toupper(cancer_type),
            main = "",
            las = 1,cex.lab = 2, cex.axis = 2,
            #make boxplot lines thicker
            boxlwd = 4,medlwd=6,whisklwd=3,staplelwd=3)
    #make overall plot box thicker
    box(lwd = 3)
    #make tick marks thicker but don't redo labels
    axis(side = 1, tick=FALSE, labels = FALSE)
    axis(side = 2, lwd.ticks = 3, labels = FALSE)
    
    #y label depends upon if there was log or not
    ylab_str <- ""
    if(grepl("mutation",ss_names[k-1],ignore.case = TRUE))
      ylab_str <- paste0(ss_names[k-1],"\n","Log10(mutation number)")
    else
      ylab_str <- ss_names[k-1]
    
    title(ylab = ylab_str,line = 6, cex.lab = 2.5)
    
    #create string for wilcox p value to put in plot box
    wilcox_p_str <- ""
    if(wilcox_p < .001)
      wilcox_p_str <- "P < .001"
    else
      wilcox_p_str <- paste0("P = ", round(wilcox_p, digits = 3))
    
    #add in Wilcoxon p value inside plot box
    text(x=.35, y=min(curr_stat, na.rm = T), pos = 4, cex = 1.5, offset = 2,
         labels = wilcox_p_str)
    
    #a table with the number of occurences for each
    occurences <- table(smoke)
    
    #add in the specific data points
    points(x = jitter(rep(1, occurences[[1]]), amount = .05),
           curr_stat[smoke == "Non-Smoking"],
           pch = 15, col = col_yb[1], cex = 1.5)
    points(x = jitter(rep(2, occurences[[2]]), amount = .05),
           curr_stat[smoke == "Smoking"],
           pch = 15, col = col_yb[2], cex = 1.5)
    
    dev.off() #close boxplot
  }
  
  snum_df <- data.frame(.tp, .wp, row.names = colnames(summ_stats_data)[2:ncol(summ_stats_data)])
  write.csv(snum_df,file="smoking_p.csv",quote=F,row.names = T)
  
} else{ print("No smoking data!") }

print("Completed all analyses!")
