rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
require(methods)

if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
if(!require(stringr)){
  install.packages("stringr")
  library(stringr)
}
if(!require(cowplot)){
  install.packages("cowplot")
  library(cowplot)
}


### set number of decimal places to use - mostly for formatting purposes
options(digits=3)

cnv <- read.table("xhmmCNV.xcnv", sep = "\t", header = TRUE)
intv <- read.table(args[1], sep = "\t")
colnames(intv) <- c("chr","start","stop","exon")
aux <- read.table("xhmmCNV.aux_xcnv", sep = "\t", header = TRUE)
freq.list <- read.table("BC1958_freqentCNVs_5pct.txt")
freq.list <- as.character(freq.list$V1)



f.aux <- aux[!aux$TARGET_IND == "U-2" & !aux$TARGET_IND == "U-1" & !aux$TARGET_IND == "D+1" & !aux$TARGET_IND == "D+2",]
intv$id <- paste(intv$chr,":", intv$start, "-", intv$stop, sep="")

x <- merge(f.aux, intv, by.x = "TARGET", by.y = "id", all.x = TRUE)
x$cnv_id <- as.numeric(as.factor(x$FULL_INTERVAL))
x <- cbind(x$cnv_id,x$SAMPLE,x$CHR,x$TARGET,x$FULL_INTERVAL,x$CNV,x$exon,x[,8:10])
colnames(x) <- c("CNV_ID","SAMPLE","CHR","TARGET","FULL_INTERVAL","CNV","EXON","POSTERIOR","RD","ORIG_RD")
###remove unmapped exons from interval files & drop unused factor levels left by merge
x <- x[!is.na(x$EXON),]
x <- droplevels.data.frame(x)

###QC check for sample/exon uniqueness - should  be a very low number
#for(i in unique(x$EXON)){
#  if(length(x$SAMPLE[x$EXON == i]) != length(unique(x$SAMPLE[x$EXON == i])))print(i)
#}

###add Q_SOME field from CNV file and number of targets per full interval originally called from CNV file
q.value <- data.frame(x=rep(0,nrow(x)))
tar.value <- data.frame(x=rep(0,nrow(x)))
for(i in 1:nrow(x)){
    q <- cnv[cnv$SAMPLE == as.character(x[i,2]) & cnv$INTERVAL == as.character(x[i,5]),]
    q.value[i,] <- q[10]
    tar.value[i,] <- q[8]
}
q.tar <- cbind(q.value, tar.value)
colnames(q.tar) <- c("Q_SOME", "TARGETS")
x <- cbind(x,q.tar)
rm(q.value,tar.value,q.tar,q,i)


###Generate exon call frequencys
af.table <- as.data.frame(table(x$EXON))
colnames(af.table) <- c("EXON", "FREQ_all")
af.table$AF_all <- af.table$FREQ_all / length(unique(x$SAMPLE))

af.dup <- as.data.frame(table(x$EXON[x$CNV == "DUP"]))
af.table$AF_dup <- af.dup$Freq / length(unique(x$SAMPLE))
af.del <- as.data.frame(table(x$EXON[x$CNV == "DEL"]))
af.table$AF_del <- af.del$Freq / length(unique(x$SAMPLE))
af.table$CNV_concensus <- as.integer((af.table$AF_dup / af.table$AF_all) * 100)

###Remove Posterior/RD/FULL_interval col
x <- merge(x, af.table, by = "EXON")
x <- x[,-(c(6,8:9))]

###FOR REF data sets - generate list of exons from CNV at more than percent freq
#freq.list <- as.vector(unique(x$EXON[x$AF_all > 0.05]))
#write.table(freq.list, file = "BC1958_freqentCNVs_5pct.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

###all exons with AF more than 0.05 in analysis cohort - removed for now - causes loss of all rows and errors if no output is < value
#x <- x[x$AF_all < 0.05,]
###remove commonly altered exons in BC1958 Cohort
if(exists("freq.list")){
  x <- x[!(x$EXON %in% freq.list),]
}


###string split Exon into gene and exon
gene_exon <- as.data.frame(str_split(as.character(x$EXON), "_", simplify = TRUE),stringsAsFactors = FALSE)
gene_exon <- gene_exon[-3]
colnames(gene_exon) <- c("GENE","EXON")
x <- x[-1]
x <- cbind(gene_exon, x)

###Filter out calls with both DUP and DEL calls - @ 10th and 90th percentile
x <- x[x$CNV_concensus > 80 | x$CNV_concensus < 20,]

write.table(x, file="cnv_xhmm_annotated.tsv", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "-9")

save.image(file="cnvANNO.RData")
