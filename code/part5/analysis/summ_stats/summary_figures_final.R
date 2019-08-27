#path to summary statistic averaged data
path <- "summ_stats_analysis_final"

#file names
filelist <- list.files(path = path, pattern = ".csv", recursive = TRUE, full.names = TRUE)

#summ_stats_data csv files
ssd_filelist <- filelist[grepl("summ_stats_data", filelist)]

#get list of the summ_stats_data matrices
datalist <- lapply(ssd_filelist, FUN = read.csv, header=TRUE, 
                   stringsAsFactors = FALSE)

#combine all cancer summ stats data into one
all_summ_stats_data <- datalist[[1]]
for(i in 2:length(datalist)){
  curr_data <- datalist[[i]]
  all_summ_stats_data <- rbind(all_summ_stats_data,curr_data)
}

#smoking data
smoking_filelist <- filelist[grepl("smoking",filelist)]

#get list of the matrices
smoking_datalist <- lapply(smoking_filelist, FUN = read.csv, header=TRUE, 
                           stringsAsFactors = FALSE, row.names = 1)

#get cancer types with smoking data
smoking_cancer_types <- c()
for(i in 1:length(smoking_datalist)){
  smoking_cancer_types <- c(smoking_cancer_types,strsplit(smoking_filelist[i],"/")[[1]][2])
}

#survival data
survival_filelist <- filelist[grepl("survival",filelist)]

#get list of the matrices
survival_datalist <- lapply(survival_filelist, FUN = read.csv, header=TRUE, 
                            stringsAsFactors = FALSE,row.names = 1)

#get cancer types with survival data
survival_cancer_types <- c()
for(i in 1:length(survival_datalist)){
  survival_cancer_types <- c(survival_cancer_types,strsplit(survival_filelist[i],"/")[[1]][2])
}
#don't want duplicates
survival_cancer_types <- survival_cancer_types[!duplicated(survival_cancer_types)]


#overall output directory
out_dir <- "final_figures"
if(!dir.exists(out_dir))
  dir.create(out_dir)
setwd(out_dir)

#corresponding summary statistic names
ss_names <- c("Depth", "Maximum Mutations", "Number of Branches", "Number of Nodes",
              "Total Mutations", "Trunk Proportion", "Clonal Expansion Index",
              "Trunk Proportion CNV")

#set up colors of yellow and blue
col_yb <- c("#FFC400", "#0033A0")

################################################################################

#distribution graphs

#subdir for distribution graphs
d_dir <- "distribution"
if(!dir.exists(d_dir))
  dir.create(d_dir)
setwd(d_dir)

#need to match cancer types up with data
cancer_type <- c()
for(i in 1:length(datalist)){
  #cancer type starts after /
  curr_cancer <- strsplit(ssd_filelist[i],"/")[[1]][2]
  #need to repeat the cancer type once for every patient in the cancer type
  cancer_type <- c(cancer_type,rep(curr_cancer,nrow(datalist[[i]])))
}


#for each summary statistic
for(n in 2:ncol(all_summ_stats_data)){
  #set colors up
  col <-c("orangered", "firebrick4", "palegreen2", "coral3", "orange4", "red", 
          "blueviolet", "navajowhite1", "darkorchid", "red2", "palegreen4", 
          "lightgreen", "darkred", "orange2", "maroon4", "brown1", "gold", 
          "deepskyblue4", "darkolivegreen3", "lightskyblue3", "orange", "orchid3", 
          "deeppink", "navajowhite4", "navajowhite3", "maroon2", "lightsalmon1", 
          "red4", "darkorchid4", "lightsteelblue4", "cyan4", "chocolate1", "chocolate4")
  
  #need to do log10 graphs for mutation number
  if(colnames(all_summ_stats_data)[n] %in% c("max_muts_avg","total_muts_avg")){
    #get data relating to current summary statistic
    tcga.cohort <- all_summ_stats_data[,c(1,n)]
    tcga.cohort$cohort <- cancer_type
    colnames(tcga.cohort)  <- c("Tumor_Sample_Barcode", "total", "cohort")
    tcga.cohort$total = as.numeric(as.character(tcga.cohort$total))
    
    #remove na values
    tcga.cohort <- tcga.cohort[complete.cases(tcga.cohort),]
    
    # calculate median and then order projects based on median
    N <- c()
    med <- c()
    for(c in 1:length(unique(cancer_type))){
      curr_cancer <- unique(cancer_type)[c]
      N <- c(N, sum(cancer_type == curr_cancer))
      curr_stat <- tcga.cohort$total[which(tcga.cohort$cohort == curr_cancer)]
      med <- c(med, median(curr_stat))
    }
    tcga.cohort.med <- cbind.data.frame(unique(cancer_type), N, med)
    tcga.cohort.med = tcga.cohort.med[order(med, decreasing = F),]
    tcga.cohort$cohort = factor(x = tcga.cohort$cohort,levels = tcga.cohort.med[,1])
    colnames(tcga.cohort.med) = c('Cohort', 'Cohort_Size', 'Median')
    tcga.cohort$TCGA = 'TCGA'
    # split into list
    tcga.cohort = split(tcga.cohort, as.factor(tcga.cohort$cohort))
    
    #order colors based on median
    #this way the same color is always with the same cancer
    col <- col[order(med, decreasing = F)]
    
    # within each project order samples based on total
    plot.dat = lapply(seq_len(length(tcga.cohort)), function(i){
      x = tcga.cohort[[i]]
      x = data.table::data.table(rev(seq(i-1, i, length.out = nrow(x))),
                                 x[order(x$total, decreasing = T), "total"],
                                 x[,"TCGA"])
      x
    })
    
    names(plot.dat) = names(tcga.cohort)
    y_lims = range(log10(data.table::rbindlist(l = plot.dat)[,V2]))
    y_max = ceiling(max(y_lims))
    y_min = floor(min(y_lims))
    y_lims = c(y_min, y_max)
    y_at = pretty(y_lims)
    
    pdf(paste0(colnames(all_summ_stats_data)[n], ".pdf"))
    par(mar = c(4, 3, 1, 1))
    plot(NA, NA, xlim = c(0, length(plot.dat)), ylim = y_lims, axes = FALSE, xlab = NA, ylab = NA)
    rect(xleft = seq(0, length(plot.dat)-1, 1), ybottom = min(y_lims), xright = seq(1, length(plot.dat), 1),
         ytop = y_max, col = grDevices::adjustcolor(col = 'white', alpha.f = 0.2),
         border = NA)
    #abline(h = pretty(y_lims), lty = 2, col = "gray70")
    # col = c('gray70', 'black')
    # color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    # set.seed(1)
    # col = sample(color,33)
    lapply(seq_len(length(plot.dat)), function(i){
      x = plot.dat[[i]]
      points(x$V1, log10(x$V2), pch = 16, cex = 0.4, col = col[i])
    })
    axis(side = 2, at = y_at, las = 2, line = -1, tick = FALSE)
    axis(side = 1, at = seq(0.5, length(plot.dat)-0.5, 1), labels = toupper(names(plot.dat)),
         las = 2, tick = FALSE, line = -1, cex = .5)
    mtext(text = paste("log10", ss_names[n-1]), side = 2, line = 1.5)
    tcga.cohort.med$Median_log10 = log10(tcga.cohort.med$Median)
    medianCol='red'
    lapply(seq_len(nrow(tcga.cohort.med)), function(i){
      segments(x0 = i-1, x1 = i, y0 = tcga.cohort.med[i, "Median_log10"],
               y1 = tcga.cohort.med[i, "Median_log10"], col = col[i])
    })
    dev.off()
    
  }
  #no log10 necessary
  else{
    #get data relating to current summary statistic
    tcga.cohort <- all_summ_stats_data[,c(1,n)]
    tcga.cohort$cohort <- cancer_type
    colnames(tcga.cohort)  <- c("Tumor_Sample_Barcode", "total", "cohort")
    tcga.cohort$total = as.numeric(as.character(tcga.cohort$total))
    
    #remove na values
    tcga.cohort <- tcga.cohort[complete.cases(tcga.cohort),]
    
    # calculate median and then order projects based on median
    N <- c()
    med <- c()
    for(c in 1:length(unique(cancer_type))){
      curr_cancer <- unique(cancer_type)[c]
      N <- c(N, sum(cancer_type == curr_cancer))
      curr_stat <- tcga.cohort$total[which(tcga.cohort$cohort == curr_cancer)]
      med <- c(med, median(curr_stat))
    }
    tcga.cohort.med <- cbind.data.frame(unique(cancer_type), N, med)
    tcga.cohort.med = tcga.cohort.med[order(med, decreasing = F),]
    tcga.cohort$cohort = factor(x = tcga.cohort$cohort,levels = tcga.cohort.med[,1])
    colnames(tcga.cohort.med) = c('Cohort', 'Cohort_Size', 'Median')
    tcga.cohort$TCGA = 'TCGA'
    # split into list
    tcga.cohort = split(tcga.cohort, as.factor(tcga.cohort$cohort))
    
    #order colors based on median
    #this way the same color is always with the same cancer
    col <- col[order(med, decreasing = F)]
    
    # within each project order samples based on total
    plot.dat = lapply(seq_len(length(tcga.cohort)), function(i){
      x = tcga.cohort[[i]]
      x = data.table::data.table(rev(seq(i-1, i, length.out = nrow(x))),
                                 x[order(x$total, decreasing = T), "total"],
                                 x[,"TCGA"])
      x
    })
    
    names(plot.dat) = names(tcga.cohort)
    y_lims = range(data.table::rbindlist(l = plot.dat)[,V2])
    y_max = ceiling(max(y_lims))
    y_min = floor(min(y_lims))
    #some individual adjustments to y_lims for some summ_stats
    if(colnames(all_summ_stats_data)[n] == "depth_avg"){
      y_lims <- c(y_min, 3.9) #depth
    }else if(colnames(all_summ_stats_data)[n] == "num_branches_avg"){
      y_lims <- c(y_min, 6.9) #num_branches
    }else if(colnames(all_summ_stats_data)[n] == "num_nodes_avg"){
      y_lims <- c(y_min, 7.9) #num_nodes
    }else if(colnames(all_summ_stats_data)[n] == "ce_avg"){
      y_lims <- c(.5, 1.5) #ce
    }else{
      y_lims = c(y_min, y_max)
    }
    y_at = pretty(y_lims)
    
    pdf(paste0(colnames(all_summ_stats_data)[n], ".pdf"))
    par(mar = c(4, 4, 1, 1))
    plot(NA, NA, xlim = c(0, length(plot.dat)), ylim = y_lims, axes = FALSE, xlab = NA, ylab = NA)
    rect(xleft = seq(0, length(plot.dat)-1, 1), ybottom = min(y_lims), xright = seq(1, length(plot.dat), 1),
         ytop = y_max, col = grDevices::adjustcolor(col = 'white', alpha.f = 0.2),
         border = NA)
    #abline(h = pretty(y_lims), lty = 2, col = "gray70")
    # col = c('gray70', 'black')
    # color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    # set.seed(1)
    # col = sample(color,33)
    lapply(seq_len(length(plot.dat)), function(i){
      x = plot.dat[[i]]
      points(x$V1, x$V2, pch = 16, cex = 0.4, col = col[i])
    })
    axis(side = 2, at = y_at, las = 2, line = -1, tick = FALSE)
    axis(side = 1, at = seq(0.5, length(plot.dat)-0.5, 1), labels = toupper(names(plot.dat)),
         las = 2, tick = FALSE, line = -1, cex = .5)
    mtext(text = ss_names[n-1], side = 2, line = 1.5)
    #tcga.cohort.med[, Median_Mutations_log10 := log10(Median_Mutations)]
    medianCol='red'
    lapply(seq_len(nrow(tcga.cohort.med)), function(i){
      segments(x0 = i-1, x1 = i, y0 = tcga.cohort.med[i, "Median"],
               y1 = tcga.cohort.med[i, "Median"], col = col[i])
    })
    dev.off()
    
  }
  
}

#reset to overall out dir
setwd("..")


############################################################################

#pairwise correlation

pc_dir <- "pairwise_cor_results"
if(!dir.exists(pc_dir))
  dir.create(pc_dir)
setwd(pc_dir)

#loop through the cols and graph them pairwise
#first col is patient ids
for(i in 2:ncol(all_summ_stats_data)){
  for(j in 2:ncol(all_summ_stats_data)){
    #only want to do each comparison once
    #e.g. not depth vs nodes and nodes vs depth
    if(i > j){
      #two stats to compare
      x_stat <- all_summ_stats_data[,i]
      y_stat <- all_summ_stats_data[,j]
      
      #need to use log scale for max and total mutations
      if(grepl("mutation",ss_names[i-1],ignore.case = TRUE))
        x_stat <- log10(x_stat)
      if(grepl("mutation",ss_names[j-1],ignore.case = TRUE))
        y_stat <- log10(y_stat)
      
      #correlation value
      cor_value <- cor(x_stat,y_stat,use = "complete.obs", method = "spearman")
      
      #filename is summ_stat1_vs_summ_stat2.jpeg
      x_name <- substr(colnames(all_summ_stats_data)[i],1,nchar(colnames(all_summ_stats_data)[i])-4)
      y_name <- substr(colnames(all_summ_stats_data)[j],1,nchar(colnames(all_summ_stats_data)[j])-4)
      filename <- paste0(y_name,"_vs_",x_name,".jpeg")
      
      #save to jpeg file
      jpeg(filename = filename,width = 1000,height = 1000,res = 100)
      
      #scatter plot
      par(mar = c(12,12,8,2))
      plot(x_stat,y_stat,
           main = "",
           xlab = "",
           ylab = "",
           #extend the top by .2 * range of y
           ylim = c(min(y_stat,na.rm=T),max(y_stat,na.rm=T)+.1*(max(y_stat,na.rm=T)-min(y_stat,na.rm=T))),
           #text of axis labels larger
           cex.lab = 2, cex.axis = 2, cex.main = 2.5, 
           #points are dots and smaller
           pch = 20, cex = .5,
           #text of all axis labels (not titles) horizontal
           las=1,
           col = col_yb[2])
      
      #add correlation value inside plot box
      text(x = min(x_stat,na.rm=T), y = max(y_stat,na.rm=T)+.1*(max(y_stat,na.rm=T)-min(y_stat,na.rm=T)),
           labels = paste0("rho = ", signif(cor_value, digits = 3)),
           adj = c(0,0.5),cex = 2)
      
      #set titles of y and x axis
      #use line parameter to distance from axis labels
      #label depends upon if there was log or not
      xlab_str <- ss_names[i-1]
      ylab_str <- ss_names[j-1]
      if(grepl("mutation",ss_names[i-1],ignore.case = TRUE))
        xlab_str <- paste0(ss_names[i-1],"\n","Log10(mutation number)")
      if(grepl("mutation",ss_names[j-1],ignore.case = TRUE))
        ylab_str <- paste0(ss_names[j-1],"\n","Log10(mutation number)")
      title(ylab = ylab_str,line = 6, cex.lab = 2)
      title(xlab = xlab_str,line = 6, cex.lab = 2)
      #make the overall plot box thicker
      box(lwd = 3)
      #make tick marks thicker but don't redo labels
      axis(side = 1, lwd.ticks = 3, labels = FALSE)
      axis(side = 2, lwd.ticks = 3, labels = FALSE)
      #close plot
      dev.off()
    }
  }
}

#reset to overall out dir
setwd("..")


############################################################################

#smoking

#output directory
smoke_dir <- "smoking_p_values"
if(!dir.exists(smoke_dir)){
  dir.create(smoke_dir,showWarnings = F)
}
setwd(smoke_dir)

#get wilcox p values
wilcox_p <- smoking_datalist[[1]][,2,drop=F]
for(i in 2:length(smoking_datalist)){
  wilcox_p <- cbind(wilcox_p,smoking_datalist[[i]][,2])
}
#set colnames to be cancer types
colnames(wilcox_p) <- smoking_cancer_types

#need to use p-value adjustment
wilcox_p[] <- lapply(wilcox_p,p.adjust, method="hochberg")

#create waterfall plots for smoking p values
#edit each p-value to be the -log()
neg_log <- function(x){return (-log10(x))}
wilcox_p[] <- lapply(wilcox_p,neg_log)

#p < .05 is cutoff so > -log(.05) is new cutoff
cutoff <- neg_log(.05)

#loop through the rows and graph by summary statistic
for(r in 1:nrow(wilcox_p)){
  #current stat to plot
  curr_stat <- as.numeric(wilcox_p[r,])
  #set cancer type names as previous line strips those
  names(curr_stat) <- toupper(smoking_cancer_types)
  #order the curr stat from largest to smallest
  curr_stat <- curr_stat[order(curr_stat,decreasing = TRUE)]
  
  #filename is summ_stat.jpeg
  filename <- paste0(rownames(wilcox_p)[r],".jpeg")
  
  #save to jpeg file
  jpeg(filename = filename,width = 1000,height = 1000,res = 100)
  
  par(mar = c(8,8,8,2))
  #bar plot
  barplot(curr_stat,
          main = ss_names[r],
          xlab = "",
          ylab = "",
          ylim = c(0,max(2,curr_stat) + .2),
          cex.names = 2, cex.main = 2.5, cex.lab = 2, cex.axis = 2,
          las = 2,
          col = col_yb[2])
  #cutoff for significance
  abline(h=cutoff, lwd = 3, lty = 2)
  text(x=10,y=cutoff, labels = "P < .05", pos = 3, cex = 2)
  
  title(ylab = "-log10 (adjusted p-value)",line = 6, cex.lab = 2)
  #make the overall plot box thicker
  box(lwd = 3)
  #make tick marks thicker but don't redo labels
  axis(side = 1, tick = FALSE, labels = FALSE)
  axis(side = 2, lwd.ticks = 3, labels = FALSE)
  

  #close plot
  dev.off()
}

#reset to overall out dir
setwd("..")

############################################################################

#survival


#get logrank p values for median split
lr_p_med <- survival_datalist[[1]][,1,drop=F]
#quartile split
lr_p_q <- survival_datalist[[2]][,1,drop=F]
#get hazard ratio values as well
#note that hazard ratio is same regardless of median or quartile
hr <- survival_datalist[[2]][,3,drop=F]

for(i in 3:length(survival_datalist)){
  if(i %% 2 != 0){
    lr_p_med <- cbind(lr_p_med,survival_datalist[[i]][,1])
    hr <- cbind(hr,survival_datalist[[i]][,3])
  }else
    lr_p_q <- cbind(lr_p_q,survival_datalist[[i]][,1])
}
#set colnames to be cancer types
colnames(lr_p_med) <- survival_cancer_types
colnames(lr_p_q) <- survival_cancer_types
colnames(hr) <- survival_cancer_types

#need to use p-value adjustment
lr_p_med[] <- lapply(lr_p_med,p.adjust, method="hochberg")
lr_p_q[] <- lapply(lr_p_q,p.adjust, method="hochberg")

#edit each p-value to be the -log()
lr_p_med[] <- lapply(lr_p_med,neg_log)
lr_p_q[] <- lapply(lr_p_q,neg_log)

#output directory
surv_dir <- "survival_p_values"
if(!dir.exists(surv_dir)){
  dir.create(surv_dir,showWarnings = F)
}
setwd(surv_dir)

#loop through the rows and graph by summary statistic
#median split
m_dir <- "median"
if(!dir.exists(m_dir)){
  dir.create(m_dir,showWarnings = F)
}
setwd(m_dir)
for(r in 1:nrow(lr_p_med)){
  #need to identify colors to use
  #yellow = low > high for the survival (hr > 1)
  #blue = high > low for the survival (hr < 1)
  #FALSE + 1 -> 1 -> yellow;    TRUE + 1 -> 2 -> blue
  bar_colors <- col_yb[(as.numeric(hr[r,]) < 1) + 1]
  
  #current stat to plot
  curr_stat <- as.numeric(lr_p_med[r,])
  #set cancer type names as previous line strips those
  names(curr_stat) <- toupper(survival_cancer_types)
  #order the curr stat from largest to smallest
  #make sure to order bar_colors first so they match up
  bar_colors <- bar_colors[order(curr_stat,decreasing = TRUE)]
  curr_stat <- curr_stat[order(curr_stat,decreasing = TRUE)]
  
  #filename is summ_stat.jpeg
  filename <- paste0(rownames(lr_p_med)[r],".jpeg")
  
  #save to jpeg file
  jpeg(filename = filename,width = 1000,height = 1000,res = 100)
  
  par(mar = c(8,8,8,2))
  #bar plot
  #use ylim to ensure cutoff is always in the y scale
  barplot(curr_stat,
          main = ss_names[r],
          xlab = "",
          ylab = "",
          ylim = c(0,max(2,curr_stat) + .2),
          cex.names = 1.2, cex.main = 2.5, cex.lab = 2, cex.axis = 2,
          las = 2,
          col = bar_colors)
  
  #cutoff for significance
  abline(h=cutoff, lwd = 3, lty = 2)
  text(x=33,y=cutoff, labels = "P < .05", pos = 3, cex = 2)
  
  title(ylab = "-log10 (adjusted p-value)",line = 6, cex.lab = 2)
  #make the overall plot box thicker
  box(lwd = 3)
  #make tick marks thicker but don't redo labels
  axis(side = 1, tick = FALSE, labels = FALSE)
  axis(side = 2, lwd.ticks = 3, labels = FALSE)
  
  #legend in the top right
  legend('topright', c("Low survived better than high", "High survived better than low"),
         col = col_yb, pch = 15, cex = 2, bty = 'n')
  
  
  #close plot
  dev.off()
}

#reset directory
setwd("..")

#quartiles
q_dir <- "quartile"
if(!dir.exists(q_dir)){
  dir.create(q_dir,showWarnings = F)
}
setwd(q_dir)
for(r in 1:nrow(lr_p_q)){
  #need to identify colors to use
  #black = low > high for the survival (hr > 1)
  #red = high > low for the survival (hr < 1)
  #FALSE + 1 -> 1 -> "red";    TRUE + 1 -> 2 -> "black"
  bar_colors <- col_yb[(as.numeric(hr[r,]) < 1) + 1]
  
  #current stat to plot
  curr_stat <- as.numeric(lr_p_q[r,])
  #set cancer type names as previous line strips those
  names(curr_stat) <- toupper(survival_cancer_types)
  #order the curr stat from largest to smallest
  bar_colors <- bar_colors[order(curr_stat,decreasing = TRUE)]
  curr_stat <- curr_stat[order(curr_stat,decreasing = TRUE)]
  
  #filename is summ_stat.jpeg
  filename <- paste0(rownames(lr_p_q)[r],".jpeg")
  
  #save to jpeg file
  jpeg(filename = filename,width = 1000,height = 1000,res = 100)
  
  par(mar = c(8,8,8,2))
  #bar plot
  #use ylim to ensure cutoff is always in the y scale
  barplot(curr_stat,
          main = ss_names[r],
          xlab = "",
          ylab = "",
          ylim = c(0,max(2,curr_stat) + .2),
          cex.names = 1.2, cex.main = 2.5, cex.lab = 2, cex.axis = 2,
          las = 2,
          col = bar_colors)
  
  #cutoff for significance
  abline(h=cutoff, lwd = 3, lty = 2)
  text(x=33,y=cutoff, labels = "P < .05", pos = 3, cex = 2)
  
  title(ylab = "-log10 (adjusted p-value)",line = 6, cex.lab = 2)
  #make the overall plot box thicker
  box(lwd = 3)
  #make tick marks thicker but don't redo labels
  axis(side = 1, tick = FALSE, labels = FALSE)
  axis(side = 2, lwd.ticks = 3, labels = FALSE)
  
  #legend in the top right
  legend('topright', c("Low survived better than high", "High survived better than low"),
         col = col_yb, pch = 15, cex = 2, bty = 'n')
  
  
  #close plot
  dev.off()
}

#reset directory
setwd("..")
setwd("..")
