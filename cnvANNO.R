rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
require(methods)

if(!require(ggplot2)){
  install.packages("ggplot2",repos = "https://mirrors.ebi.ac.uk/CRAN/")
  library(ggplot2)
}
if(!require(stringr)){
  install.packages("stringr", repos = "https://mirrors.ebi.ac.uk/CRAN/")
  library(stringr)
}

if(!require(tidyr)){
  install.packages("tidyr", repos= "https://mirrors.ebi.ac.uk/CRAN/")
  library(tidyr)
}
if(!require(dplyr)){
  install.packages("dplyr", repos= "https://mirrors.ebi.ac.uk/CRAN/")
  library(dplyr)
}
#### VARIABLE OPTIONS! ####

int_af_value <- 0.05
ref_af_value <- 0.05

#### VARIABLE OPTIONS ####

# ---------------- DO NOT CHANGE BELOW THIS LINE UNLESS YOU KNOW WHAT YOU'RE DOING! ---------------- #

### set number of decimal places to use - mostly for formatting purposes
options(digits=3)
###Read in files for CNV annotation script
cnv <- read.table("xhmmCNV.xcnv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
intv <- read.table(args[1], sep = "\t", stringsAsFactors = FALSE)
colnames(intv) <- c("chr","start","stop","exon")
aux <- read.table("xhmmCNV.aux_xcnv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#NONREF#
ref.list <- read.table("BC1958_CNVs.txt", sep = "\t", stringsAsFactors = FALSE)
colnames(ref.list) <- c("EXON","CNV","AF_ref")



f.aux <- aux[!aux$TARGET_IND == "U-2" & !aux$TARGET_IND == "U-1" & !aux$TARGET_IND == "D+1" & !aux$TARGET_IND == "D+2",]
intv$id <- paste(intv$chr,":", intv$start, "-", intv$stop, sep="")

x <- merge(f.aux, intv, by.x = "TARGET", by.y = "id", all.x = TRUE)
x$cnv_id <- as.numeric(as.factor(x$FULL_INTERVAL))
x <- cbind(x$cnv_id,x$SAMPLE,x$CHR,x$TARGET,x$FULL_INTERVAL,x$CNV,x$exon,x[,8:10])
colnames(x) <- c("CNV_ID","SAMPLE","CHR","TARGET","FULL_INTERVAL","CNV","EXON","POSTERIOR","RD","ORIG_RD")
###remove unmapped exons from interval files & drop unused factor levels left by merge
x <- x[!is.na(x$EXON),]
x <- droplevels.data.frame(x)


###add Q_SOME field from CNV file and number of targets per full interval originally called from CNV file
q.value <- data.frame(x=rep(0,nrow(x)))
for(i in 1:nrow(x)){
    q <- cnv[cnv$SAMPLE == as.character(x[i,2]) & cnv$INTERVAL == as.character(x[i,5]),]
    q.value[i,] <- q[10]
}

q.tar <- q.value
colnames(q.tar) <- "Q_SOME"
x <- cbind(x,q.tar)
rm(q.value,q.tar,q,i)

##Remove redundant columns
x <- x[,-c(5,8:9)]

##Conversion to vcf sytyle genotype annotation
t <- unique(x[c("EXON","CNV")])
t <- cbind(t,seq.int(1,nrow(t),1))
colnames(t)[3] <- "EXON_CNV_ID"
x <- merge(x,t,by.y = c("CNV","EXON"), all.x = TRUE)
##adding genotype columns for each sample - spread function from dplyr - removal of factors to ease data handling

x$GT <- 1
x <- x[!duplicated(x),]
x <- spread(x,SAMPLE,GT,fill=0)
x %>% mutate_if(is.factor, as.character) -> x
#na replaced as 0 in ref genotype field
x[is.na(x)] <- 0

##selecting columns with constant values across rows & collapsing into single values
x_const <- x[,c(8,1,2,4,5)]
x_const <- x_const %>% group_by(EXON_CNV_ID) %>% summarise_all(funs(paste(unique(.), collapse=",")))
##Selecting columns with variable values and concatenating them  into single cells
x_var <- x[,c(8,3,6,7)]
x_var <- x_var %>% group_by(EXON_CNV_ID) %>% summarise_all(funs(paste(., collapse=",")))
##detaching genotype information and collapsing them into single row for each unique "site" - Idea by EF!
x_geno <- x[,c(8:ncol(x))]
x_geno <- as.data.frame(x_geno %>% group_by(EXON_CNV_ID) %>% summarise_all(funs(sum(as.numeric(.)))))
#confirming only 1 or 0 present
x_geno[-1][x_geno[-1] > 0] <- 1

##Reconstructing db into single dataframe
x <- cbind(x_const,x_var,x_geno)
x <- x[-c(6,10)]

rm(x_const,x_geno,x_var,t)

##adding mean Q_some for each row
x$Mean_Q_Some <- sapply(str_split(x$Q_SOME, ","),function(x) mean(as.numeric(x)))
x$Mean_Orig_RD <- sapply(str_split(x$ORIG_RD, ","),function(x) mean(as.numeric(x)))
x <- cbind(x[1:7],x[ncol(x)],x[9:ncol(x)-1])
x <- cbind(x[1:9],x[ncol(x)],x[11:ncol(x)-1])

##Addition of AF internal to file
AF_all <- apply(x[11:ncol(x)],1, function(y) (sum(y == 1)/sum(y == 0)))
x <- cbind(x[1:10],AF_all,x[11:ncol(x)])

#FOR REF# data sets - generate list of exons from CNV at more than percent freq
#freq.list <- x[,c("EXON","CNV","AF_all")]
#write.table(freq.list, file = "BC1958_freqentCNVs_5pct.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)



#NONREF#
#adding REF_AF
x <- merge(x, ref.list, by = c("EXON","CNV"), all.x = TRUE, fill = 0)
x[is.na(x)] <- 0
x <- cbind(x[1:11],x[ncol(x)],x[13:ncol(x)-1])
#remove commonly altered exons in BC1958 Cohort
x <- x[x$AF_ref < ref_af_value,]

#NONREF#
###all exons with AF more than 0.05 in analysis cohort - removed for now - causes loss of all rows and errors if no output is < value
x <- x[x$AF_all < int_af_value,]

###string split Exon into gene and exon
gene_exon <- as.data.frame(str_split(as.character(x$EXON), "_", simplify = TRUE),stringsAsFactors = FALSE)
gene_exon <- gene_exon[-3]
colnames(gene_exon) <- c("GENE","EXON")
x <- x[-1]
x <- cbind(gene_exon,x[1:ncol(x)])

## Make sure chr positions are unified as numeric - not containing "chr"
x$CHR <- gsub("chr","",x$CHR)
x$TARGET <- gsub("chr","",x$TARGET)
## write output
write.table(x, file="cnv_xhmm_annotated.tsv", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "-9")
#NONREF#
rm(aux,cnv,f.aux,gene_exon,intv,ref.list,AF_all,args,int_af_value,ref_af_value)
save.image(file="cnvANNO.RData")
