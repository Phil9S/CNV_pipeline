## CNV Pipeline - Olive Gold Volt
#### Input and environment setup ####
rm(list = ls())
## load data from cnvANNO.R script
load(file="cnvANNO.RData")

## load required libraries
require(methods)

if(!require(dplyr)){
  install.packages("dplyr",repos = "https://mirrors.ebi.ac.uk/CRAN/")
  library(dplyr)
}

if(!require(tidyr)){
  install.packages("tidyr",repos = "https://mirrors.ebi.ac.uk/CRAN/")
  library(tidyr)
}

if(!require(ggplot2)){
  install.packages("ggplot2",repos = "https://mirrors.ebi.ac.uk/CRAN/")
  library(ggplot2)
}

if(!require(stringr)){
  install.packages("stringr",repos = "https://mirrors.ebi.ac.uk/CRAN/")
  library(stringr)
}

if(!require(reshape2)){
  install.packages("reshape2", repos= "https://mirrors.ebi.ac.uk/CRAN/")
  library(reshape2)
}
if(!require(data.table)){
  install.packages("data.table", repos= "https://mirrors.ebi.ac.uk/CRAN/")
  library(data.table)
}
if(!require(ggrepel)){
  install.packages("ggrepel", repos= "https://mirrors.ebi.ac.uk/CRAN/")
  library(ggrepel)
}
if(!require(circlize)){
  install.packages("circlize", repos= "https://mirrors.ebi.ac.uk/CRAN/")
  library(circlize)
}
## 
## turn off use of scientific number notation
options(scipen=999)

#### p1 ####
## Bar plot of CNV types and counts per Sample called

## import a summation of the number of DELs or DUPs per sample using dplyr group_by and summarise
CNV_samp <- as.data.frame(x[14:ncol(x)] %>% group_by(x$CNV) %>% summarise_all(funs(sum(as.numeric(.)))))
## retain sample names for cbind downstream as names are lost on transpostion step
SAMPLES <- colnames(CNV_samp[-1])
CNV_samp <- as.data.frame(t(CNV_samp)[-1,]) ## transpose data frame from x,y to y,x
colnames(CNV_samp) <- c("DEL","DUP") ## restore col names
## remove factors and add col for sample names back to data.frame
CNV_samp %>% mutate_if(is.factor, as.numeric) -> CNV_samp
CNV_samp <- cbind(SAMPLES,CNV_samp)
## Melt alters the table from wide format to long format - easier for plotting
CNV_samp <- melt(CNV_samp,id.vars = "SAMPLES")
## Addition of an inx field consisting of the sample range for each line - used in facet_grid in the plot to 
## make labels for each child plot generated - basically just a series of nested seq.ints and reps
CNV_samp$inx_split <- rep(paste(rep("Samples", times= nrow(CNV_samp)/2),
                                paste(rep(seq.int(1,1001,10), times=1, each=10)[1:(nrow(CNV_samp)/2)],rep(" ",times=nrow(CNV_samp)/2),
                                rep("-", times=nrow(CNV_samp)/2),rep(" ",times=nrow(CNV_samp)/2),
                                rep(seq.int(10,1000,10), times=2, each=10)[1:(nrow(CNV_samp)/2)],sep = "")), times=2)
## form factor from inx_split col so that facets are ordered
CNV_samp$inx_split <- factor(CNV_samp$inx_split, levels = unique(CNV_samp$inx_split))
## sum of del/dups per sample - used for stat calcs 
total <- CNV_samp %>% group_by(SAMPLES) %>% summarise(sum(value))
total_sd <- sd(total$`sum(value)`) ## sdev calc
total_mean <- mean(total$`sum(value)`) ## mean calc
ylim <- max(ceiling(total$`sum(value)`/10)*10) ## ylim maximum calc

##ggplot line to generate the graph
p1 <- ggplot(CNV_samp, aes(SAMPLES,value, fill=variable)) + geom_col(colour = c("gray12")) + ## the actual data being plotted
## theme modifications to make things pretty
    theme(panel.border = element_blank(), axis.line = element_line(colour="black"), panel.grid.major.y = element_blank()) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 18,colour = "gray10",margin = margin(0,0,10,0))) +
    theme(legend.key = element_blank(), legend.position = "right", axis.title = element_text(size = 14)) +
    theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(), plot.caption = element_text(size = 11, colour = "gray10", face = "italic")) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "bottom") +
    theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14), legend.text.align=0.5) +
## the axis titles, graph title and caption
    labs(list(title = "Exons called as CNV per sample", y = "Count", fill = "CNV Type", caption = "Split into groups of 10")) +
## Bar colours for the data
    scale_fill_manual(values=c("steelblue3","goldenrod2")) +
## Modifying the x and y scales to customise break points and limits
    scale_x_discrete(labels = c(seq(1,10,1))) +
    scale_y_continuous(expand = c(0,0),limits = c(NA,ylim), breaks = seq(0,ylim,20)) +
## additional geom_lines used for mean and standard dev - includes formating and positioning
    geom_hline(aes(yintercept = total_mean + (total_sd *2), linetype="2sd (+)"), size = 0.8, color = "gray10", alpha=0.7) +
    geom_hline(aes(yintercept = total_mean, linetype="mean"), size = 0.8, color = "gray10", alpha=0.7) +
    geom_hline(aes(yintercept = total_mean - (total_sd *2), linetype="2sd (-)"), size = 0.8, color = "gray10", alpha=0.7) +
    scale_linetype_manual(name = "Metrics", breaks=c("2sd (+)","mean","2sd (-)"), values = c(4, 4, 1),
                          guide = guide_legend(override.aes = list(color = c("gray10", "gray10", "gray10"))))
## plot faceting on inx_split group to generate multiple plots split into groups of 10
p1 <- p1 + facet_wrap(~inx_split, scales = "free_x") + 
    ## Customising look of facet titles and backgrounds
    theme(strip.text = element_text(size = 14), strip.background = element_rect(fill = "gray95"))

## Plotting graph to .png file
png(file = "Exon-CNV_calls_per_sample.png", width = 32, height = 16, units = "in", res = 600)
print(p1)
dev.off()
rm(total,p1,SAMPLES,total_mean,total_sd,ylim) ## remove used and redundant variables

#### p2 ####
chr_vals <- paste("chr",c(seq.int(1,22,1),"X","Y"),sep = "")
sample_list <- colnames(x)[14:ncol(x)]

## add for loop for each sample 
for(s in sample_list){
  ## factor levels acquired for each autosomal and sex chromosome 
  y <- x[x[,s] == 1,]
  ## Selecting the required cols for analysis - positions, chr, and gene
  y <- y[c(5,3,2)]
  ##Remove any extra target values from instances of multiple targets per exon
  y[1] <- apply(y[1],1, function(x) unlist(str_split(x,","))[[1]][1])
  cnv_stat <- y[2]
  cnv_gene <- y[3]
  ## string spliting position into start and stop and reforming the dataframe
  y1 <- str_split(y$TARGET,"-",simplify = TRUE)
  y2 <- str_split(y1[,1],":",simplify = TRUE)
  y <- data.frame(chr = c(as.character(y2[,1])), start = c(as.integer(y2[,2])), stop = c(as.integer(y1[,2])))
  y <- cbind(y,cnv_stat,cnv_gene)
  ## Assigning data classes to dataframe - chr as factor for plotting is important for chr order on plot and sorting
  #y$chr <- factor(y$chr, levels = chr_vals)
  y$start <- as.numeric(y$start)
  y$stop <- as.numeric(y$stop)
  ## Ordering the dataframe sequentially by factorised chr value, start position, and then stop position
  y <- y[order(y$chr,y$start,y$stop),]
  y <- y[!duplicated(y),]
  
  ## addition of numerical values to represent del/dups
  y$CNVvalue[y$CNV == "DEL"] <- -1
  y$CNVvalue[y$CNV == "DUP"] <- 1
  
  y$chr <- paste("chr",y$chr,sep = "")
  
  png(file = paste("cnv_calls_",s,".png",sep=""), width = 8, height = 8, units = "in", res = 600)
  circos.clear()
  circos.par("track.height" = 0.15, start.degree = 90,cell.padding = c(0, 0, 0, 0),gap.degree = 5)
  
  circos.initializeWithIdeogram(species = "hg38",
                                chromosome.index = chr_vals,
                                plotType = c("ideogram","labels"),
                                ideogram.height =  .075)
  title(paste("Genomic Distribution of CNVs for sample",s,sep=" "))
  circos.genomicTrack(y[,c(1,2,3,6)],ylim = c(-1.5, 1.5),
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region, value, type = "segment", lwd = 5, 
                                            col = ifelse(value == 1,"goldenrod2","steelblue3"), ...)
                        circos.lines(CELL_META$cell.xlim,c(0,0),lty = 1, col = "#00000040")
                        circos.lines(CELL_META$cell.xlim,c(-1,-1),lty = 1, col = "#00000010")
                        circos.lines(CELL_META$cell.xlim,c(1,1),lty = 1, col = "#00000010")
  })
  for(i in chr_vals){
    circos.yaxis(side = "left",at = c(1,-1),labels = c("DUP","DEL"),tick = T,sector.index = i,labels.cex = 0.5,labels.niceFacing = T,tick.length = 0.1)
  }
  dev.off()
}
