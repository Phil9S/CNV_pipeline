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
png(file = "Exons called as CNV per sample.png", width = 32, height = 16, units = "in", res = 600)
print(p1)
dev.off()
rm(total,p1,SAMPLES,total_mean,total_sd,ylim) ## remove used and redundant variables



#### p2 ####
chr_vals <- c(seq.int(1,22,1),"X","Y")
sample_list <- colnames(x)[14:ncol(x)]

## add for loop for each sample 
for(s in sample_list){
  ## factor levels acquired for each autosomal and sex chromosome 
  y <- x[x[,s] == 1,]
  ## Selecting the required cols for analysis - positions, chr, and gene
  y <- y[c(6,3,1)]
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
  y$chr <- factor(y$chr, levels = chr_vals)
  y$start <- as.numeric(y$start)
  y$stop <- as.numeric(y$stop)
  ## Ordering the dataframe sequentially by factorised chr value, start position, and then stop position
  y <- y[order(y$chr,y$start,y$stop),]
  y <- y[!duplicated(y),]
  
  ##Form empty dataframe to hold rows describing the missing/diploid regions between called genes on each chromosome
  miss_x <- data.frame()
  
  ## for each row in y
  for(i in 1:nrow(y)){
    if(i  == nrow(y)){break} #break on last row as nothing to compare to
    if(y[i+1,1] == y[i,1] & y[i+1,5] != y[i,5] & y[i+1,2] - y[i,3] > 1000){ ## provided chr numbers match / gene names don't / size is > 1000 bases
      ## bind each new row to original data.frame - same format as the y ouput generated downstream
      miss_x <- rbind.data.frame(miss_x,c(y[i,1],        ## rbind new row to the exisiting frame iteratively
                                          y[i,3]+1,     ## position is calculated as a mean of the start and stop
                                          paste(y[i,4]),                     ## use the cnv type of the first row / not second
                                          NA,                                ## add column to be used later for "value" for plotting
                                          ## essentially a line for annotating the size of the region - if statement alters if it is scaled to kb or mb
                                          paste(if((y[i+1,2]-1) - (y[i,3]+1) > 1000000){
                                                  paste(round((y[i+1,2]-1 - y[i,3]+1)/1000000,digits = 2),"Mb")
                                                } else {
                                                  paste(round((y[i+1,2]-1 - y[i,3]+1)/1000,digits = 2),"Kb")
                                                })),
                                 ## exported without factors
                                 stringsAsFactors = FALSE)
    }
  }
  ## miss_x can be empty! Checking to see if there are any missing regions
  if(nrow(miss_x) > 0){
    ## formating, colnames, and class additions to miss_x
    colnames(miss_x) <- c("chr","pos","cnv","value","gene")
    ## coersion of X to 23 in miss_x causes NA chr values downstream
    miss_x$chr[miss_x$chr == 23] <- "X"
    miss_x$chr[miss_x$chr == 24] <- "Y"
    miss_x$value[is.na(miss_x$value) & miss_x$cnv == "DUP"] <- 0    # add in values cooresponding to cnv type
    miss_x$value[is.na(miss_x$value) & miss_x$cnv == "DEL"] <- 0   # add in values cooresponding to cnv type
    miss_x <- data.frame(chr=factor(miss_x$chr, levels = chr_vals),
                       pos=as.numeric(miss_x$pos),
                       cnv=as.character(miss_x$cnv),
                       value=as.numeric(miss_x$value),
                       gene=as.character(miss_x$gene), stringsAsFactors = FALSE)
  ## replicate lines so that the region is plotted in a viewable scale - not scaled based on region size though
    miss_x <- miss_x[rep(seq_len(nrow(miss_x)), each=500),]
  }
  ## removal of used or redunant variables
  rm(y1,y2,cnv_stat,cnv_gene,i)
  
  
  ## addition of numerical values to represent del/dups
  y$CNVvalue[y$CNV == "DEL"] <- -1
  y$CNVvalue[y$CNV == "DUP"] <- 1
  ## downstream seq generation of each base per exon retains only one extra column - all other information merged into chr col
  ## original cols that were condensed and pasted together are then removed
  y$chr <- paste(y$chr,y$CNV,y$CNVvalue,y$GENE,sep="#")
  y <- y[1:3]
  
  ## dpylr/data.table lines to take each exon start and stop and generate a line per bp/per exon
  ## expanding the size of each called exon on the plot
  ## specifies to perform the action rowwise and use each row's start and stop to generate a list of positions
  ## lists are then merged into a list of lists
  pos <- y %>% rowwise() %>% do(pos=seq(.$start, .$stop))
  pos$chr <- y$chr ## adding the merged chr col back to the list of lists
  pos <- data.table(pos) ## conversion to data.table format
  y <- pos[ , list(pos = unlist(pos)), by=chr] ## unlist the list of lists to generate a data.frame using chr information as a key
  ## split  merged information in the chr col to generate previous columns and format them
  y <- cbind(y$pos,as.data.frame(str_split(y$chr,"#",n = 4,simplify = TRUE), stringsAsFactors = FALSE))
  colnames(y) <- c("pos","chr","cnv","value","gene")
  y <- cbind(y[2],y[1],y[3:5])
  ## rm pos variable
  rm(pos)
  
  ## assigning class types for each data col in y
  y$pos <- as.numeric(y$pos)
  y$value <- as.numeric(y$value)
  y$chr <- factor(y$chr, levels = chr_vals)
  
  ## rows for missing data generated in miss_x above are rbind() to y dataframe and then sorted by chr (by factor) / pos (L to H)
  y <- rbind.data.frame(y,miss_x)
  y <- y[order(y$chr,y$pos),]
  ## addition of a int seqence along data to use as x-axis plot position
  y$int <- seq.int(1,nrow(y),1)
  
  ## Add spacing rows for each chromosome start and end - Min position
  chr_min <- y %>% group_by(chr) %>% summarise(min = min(pos)-1)
  chr_min <- chr_min[rep(seq_len(nrow(chr_min)), each=75),]
  chr_min[3:6] <- data.frame(NA,0,"#spacer#",NA)
  colnames(chr_min) <- c("chr","pos","cnv","value","gene","int")
  ## Add spacing rows for each chromosome start and end - Max position
  chr_max <- y %>% group_by(chr) %>% summarise(max = max(pos)+1)
  chr_max <- chr_max[rep(seq_len(nrow(chr_max)), each=75),]
  chr_max[3:6] <- data.frame(NA,0,"#spacer#",NA)
  colnames(chr_max) <- c("chr","pos","cnv","value","gene","int")
  ## rbind chr_min and chr_max positions to dataframe and order based on position
  y <- rbind.data.frame(y,chr_min,chr_max)
  y <- y[order(y$chr,y$pos),]
  ## reapply the int.seq to include the spacer points
  y$int <- seq.int(1,nrow(y),1)
  
  ## Generate the min and max values for each chr used in ploting geom_rect backgrounds
  min_max <- y %>% group_by(chr) %>% summarise(min = min(int), max = max(int))
  ## Generate the mean int value for each chr to align the x axis label with
  label_chr <- y %>% group_by(chr) %>% summarise(int=round(mean(int)))
  l <- y
  ## Seperate the dip values from actual calls to allow for seperate aes customisation
  ## Seperate the spacer values from actual calls to plot seperately
  miss_x <- y[grepl(pattern = " Kb",x = y$gene) | grepl(pattern = " Mb",x = y$gene),]
  ## removal of spacer and dip lines
  y <- y[!grepl(pattern = " Mb",x = y$gene),]
  y <- y[!grepl(pattern = " Kb",x = y$gene),]
  y <- y[!grepl(pattern = "#spacer#",x = y$gene),]
  ## Generate mid points of each gene region to anchor labels to
  label_y <- y %>% group_by(chr,gene) %>% summarise(int=round(mean(int)))
  ## drop labels for chr where the number of calls exceeds a value - default 20
  label_y_fil <- as.data.frame(label_y %>% group_by(chr) %>% summarise(sum=sum(length(gene))))
  label_y_fil <- as.vector(label_y_fil$chr[label_y_fil$sum > 30])
  if(length(label_y_fil) != 0){
    label_y <- label_y[!label_y$chr %in% label_y_fil,]
  }
  ## Generate mid points of each diploid region to anchor labels to
  label_miss_x <- miss_x %>% group_by(chr,gene) %>% summarise(int=round(mean(int)))
  ## drop labels for chr where the number of calls exceeds a value - default 20
  label_miss_x_fil <- as.data.frame(label_miss_x %>% group_by(chr) %>% summarise(sum=sum(length(gene))))
  label_miss_x_fil <- as.vector(label_miss_x_fil$chr[label_miss_x_fil$sum > 30])
  if(length(label_miss_x_fil) != 0){
  label_miss_x <- label_miss_x[!label_miss_x$chr %in% label_miss_x_fil,]
  }
  
  ## generate vector of colours for alternating chr geom_rects - same length as possible number of chr
  min_max_col <- rep_len(c("gray10","gray70"),length.out = 24)
  
  p2 <- ggplot() +
    ## addition of rect alternation on chr
    #geom_rect(data = min_max, aes(xmin= min, xmax= max, ymin= -Inf, ymax= Inf, fill=chr)) +
    ## plotting h-lines for del and dups
    geom_hline(yintercept = 1, colour = "darkgreen", linetype = "dotted", alpha = 0.3) +
    geom_hline(yintercept = 0, colour = "black",linetype = "dashed", alpha = 0.6) +
    geom_hline(yintercept = -1, colour = "darkred",linetype = "dotted", alpha = 0.3) +
    
    geom_line(data=l, aes(int,value)) +
    
    geom_point(data=miss_x, aes(int,value),colour = 'transparent') +
    geom_label_repel(data=subset(miss_x, miss_x$int %in% label_miss_x$int), 
                     aes(int,value,label = gene), size = 5, label.padding = unit(0.15,"lines"), nudge_y = -0.15) +
    ## adding point data for values in y & adding gene labels - nudge depending on type of CNV
    geom_point(data=y, aes(int,value,colour=gene)) +
    geom_label_repel(data=subset(y, y$int %in% label_y$int & y$cnv == "DEL"), label.padding = unit(0.15,"lines"),
                     aes(int,value,label = gene,colour=gene),nudge_y = -0.5, size = 5) +
    geom_label_repel(data=subset(y, y$int %in% label_y$int & y$cnv == "DUP"), label.padding = unit(0.15,"lines"),
                     aes(int,value,label = gene,colour=gene),nudge_y = 0.5, size = 5) +
    
    labs(list(title= paste("Genomic Distribution of CNVs for sample",s,sep=" "),
              caption = "Diploid Chromosomes not shown - not scaled accurately - spacer regions < 1000bp ignored")) +
    
    scale_y_continuous(breaks = c(-1,1),labels = c("LOSS","GAIN"),limits = c(-2,2),minor_breaks=NULL) +
    scale_x_continuous(expand = c(0,0),breaks = label_chr$int, labels = label_chr$chr) +
    scale_fill_manual(values = alpha(min_max_col,.1)) +
    
    theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 12),axis.title = element_blank()) +
    theme(plot.caption = element_text(size = 11, colour = "gray10", face = "italic")) +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0)) +
    theme(legend.position="none", panel.background = element_blank(), axis.ticks.x = element_blank()) +
    theme(axis.line = element_line(colour = "gray10"), panel.grid = element_blank())
  p2 <- p2 + facet_wrap(~chr, scales = "free_x",ncol = 1, labeller = label_both) +
             theme(strip.text = element_text(face = "bold", size = 16, hjust = 0.01))
            
  ## Plotting graph to .png file - scale png length with number of chromosomes so each are evenly sized
  png(file = paste("cnv_calls_",s,".png",sep=""), width = 15, height = length(unique(y$chr))*3, units = "in", res = 300)
  print(p2)
  dev.off()
  
  rm(chr_max,chr_min,label_chr,label_miss_x,label_y,min_max,y,min_max_col,p2,label_y_fil)
}
