##CNV Pipeline - Olive Gold Volt
rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
###set working directory and import arguments and libraries
setwd(args[1])

PCA <- args[2]
PCA <- PCA * 100

if(!require(ggplot2)){
  install.packages("ggplot2",repos = "https://mirrors.ebi.ac.uk/CRAN/")
  library(ggplot2)
}

if(!require(data.table)){
  install.packages("data.table",repos = "https://mirrors.ebi.ac.uk/CRAN/")
  library(data.table)
}


###read in filtered and centered read depth data from xhmm
t <- fread("xhmmCNV.filtered_centered.RD.txt",sep = "\t",header = TRUE)
t1 <- t[,-1] ### remove header line for matrix PCA

p <- prcomp(t1) ### perform PCA on matrix of read depths

###generate SVD eigen values from SD data in PCA output
###Coerce into a dataframe with index values for each PC
scr <- as.data.frame(p$sdev^2/sum(p$sdev^2)*100)
scr$PC <- seq.int(nrow(scr))
names(scr)[1] <- "Eigen"
###find the PC at which 70% of all variance is accounted for
for(i in 1:nrow(scr)){
  c <- sum(scr$Eigen[1:i])
  if(c >= PCA){
    val <- i 
    break
    }
}
###plot the data of eigen value against PC including a marker of 70% contributed variance
png("PCA_Scree.png", width = 5, height = 5, units = 'in', res = 600)
ggplot(scr, aes(x=PC, y=Eigen)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept=val, linetype = "dashed", color="red") +
  scale_y_continuous(name="Eigen Value - Contributed Variance (%)", breaks = pretty(scr$Eigen, n = 10)) +
  scale_x_continuous(name="Princple Component", breaks = pretty(scr$PC, n = 10)) +
  ggtitle(label="PCA Scree Plot", subtitle="Cummulative Contributed Variance") +
  geom_text(data=NULL, x=val+3, y=max(scr$Eigen), label="Contributed Variance cut off", size=2.5) +
  theme(panel.border = element_blank(), axis.line = element_line(colour="black")) +
  theme(panel.background = element_blank(), panel.grid.major = element_line(size = 0.1,colour = "grey50"))
dev.off()

cat(capture.output(summary(p)), file = "PCA_summary.txt")
o <- paste("Variance threshold at Principle component:",val, sep=" ")
cat(o, file="PCA_summary.txt", append=TRUE)
