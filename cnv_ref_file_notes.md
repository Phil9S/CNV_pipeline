# Generation of BED4 file for CNV analysis & frequency lists from reference datasets

## Exon Interval file 
The interval file (specified by -v --interval) is a required file for the CNV Pipeline to run but is likely to be the cause of most common errors, run failures, and/or annotation failures. This file attempts to document how the interval file was generated and how a reference set of frequently altered genes/exons was created

#### Reference interval retrieval   
[Biomart](http://www.ensembl.org/biomart/martview/8ed472ded7ffb2e51ab8271b77fa4b96) is used as the source for the retrieval of genomic intervals for each recorded exon present in the ensembl.org databases

Genomic intervals are retrieved with wget using the biomart XML query below:

```sh
wget -c -O hg38_exons_Ensembl87.txt 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "chromosome_name" /><Attribute name = "exon_chrom_start" /><Attribute name = "exon_chrom_end" /><Attribute name = "external_gene_name" /><Attribute name = "rank" /></Dataset></Query>'
```
The subsequent file contains 5 columns, as such the 4th (gene name) & 5th (exon number/rank) are merged with underscore as a delimiter using [vim](http://www.vim.org/), and sorted in ascending order by chromosome and position using [sort](http://man7.org/linux/man-pages/man1/sort.1.html):
```sh
vim -c 's/\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)/\1\t\2\t\3\t\4_exon\5/' hg38_exons_Ensembl87.txt
sort -k1,1 -k2,2n hg38_exons_Ensembl87.txt > hg38_exons_Ensembl87.sorted.txt
```
The resulting tsv file in loaded into [R](https://www.r-project.org/) to isolate only canonical autosomes and allosomes, outputting the filtered tsv to a new file using the code below:
```R
t <- read.table("hg38_exons_Ensembl87.sorted.txt", sep="\t")
l <- unique(t$V1)
t1 <- t[t$V1 %in% l[c(1:22,349:350)],]
write.table(t1,"hg38_exons_Ensembl87.autoXY.sorted.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
```
#### Repeatmask files
An important component of the interval file used is that is appropriately filtered for regions that are **not** of interest; this included regions overlapping low-complexity, repeat masked, regions of the genome. These regions are not useful and add noise to CNV calling so can be justifibly removed from the targets. The site repeatmasker.org provides catagorical fasta files for each type and span of these regions.

Fasta files containing the type and genomic position of repeat-masked regions can be downloaded from [repeatmasker.org](http://www.repeatmasker.org/) using the following commands:
```sh
wget -c "http://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz"
gunzip hg38.fa.out.gz
```
The repeat-mask fasta is not immediately appropriate for use so some minor pre-processing steps are required to allow for interval comparisons (i.e. delimiter alterations, repeat type selection, and sorting by chromosome and position), shown below:
```sh
awk -v OFS="\t" '$1=$1' hg38.fa.out > hg38_tab.fa.out
grep "Simple_repeat" hg38_tab.fa.out >> lowcomplex_simpreps.hg38.bed
grep "Low_complexity" hg38_tab.fa.out >> lowcomplex_simpreps.hg38.bed
sort -k1,1 -k2,2n lowcomplex_simpreps.hg38.bed > lowcomplex_simpreps.sorted.hg38.bed
```
#### generating list of overlapping intervals 
Filtering out intervals with any amount of repeat-mask overlap would be overly stringent so only exome intervals harbouring an overlap of 25% or more are excluded - In line with the steps performed [here](http://www.nature.com/ng/journal/v48/n10/full/ng.3638.html)
```sh
bedtools intersect -wb -F 0.25 /
		   -a hg38_exons_Ensembl87.autoXY.sorted.bed /
		   -b lowcomplex_simpreps.sorted.hg38.bed > exon_repmask_overlap.txt
vim -c '%s/\S\+\t\S\+\t\S\+\t\(\S\+\)\t\S\+\t\S\+\t\S\+\t\S\+/\1/' exon_repmask_overlap.txt
sort -u exon_repmask_overlap.txt > exon_repmask_overlap.uniq.txt
```
```R
t <- read.table("hg38_exons_Ensembl87.autoXY.sorted.bed", sep="\t")
l <- read.table("exon_repmask_overlap.uniq.txt")
t1 <- t[!t$V4 %in% l$V1,]
write.table(t1, "hg38_exons_Ensembl87.masked.autoXY.sorted.bed", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
```
#### Removing duplicated exons and retaining exons of multiple transcripts with longest length 
###### **THIS STEP TAKES SOME TIME! SCALE OF 12+ HOURS NOT MINUTES!**
###### **This will retain duplicate exons of the same size - differing start and stop positions are present for a majority of these** 

```R
dat <- read.table("hg38_exons_Ensembl87.masked.autoXY.sorted.bed", sep = "\t")
dat$V5 <- dat[,3] - dat[,2]
exons <- unique(dat$V4)
table <- unique(dat)

tab <- data.frame()
for(i in exons){
	ex <- table[table$V4 == i,]
	if(nrow(ex) > 1){ tab <- rbind(tab, ex[ex$V5 == max(ex$V5),])
	} else {
	tab <- rbind(tab, ex)
	}
}
tab <- tab[,-5]
write.table(tab, "hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.bed", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
```

#### Adding padding to exon regions - 2bp upstream of start and 2bp downstream of end
```R
t <- read.table("hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.bed", sep="\t")
n <- t
n$V2 <- n$V2 - 2
n$V3 <- n$V3 + 2
write.table(n, "hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.pad.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
```
```sh
sort -k1,1 -k2,2n hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.pad.bed > hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.pad.bedsort
mv hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.pad.bedsort hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.pad.bed
```
Interval file for CNV analysis -> hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.pad.bed (8,875Kb)

##### Removal summary:
|File|Size|Change|Process|
|----|----|------|-------|
|hg38_exons_Ensembl87.sorted.txt| 1,285,529 targets | Original file |
|hg38_exons_Ensembl87.autoXY.sorted.bed| 1,181,871 targets| -103,658 | Removal of non-standard chromosomes|
hg38_exons_Ensembl87.masked.autoXY.sorted.bed| 969,614 targets| -212,257 | Removal of masked repeat regions|
hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.pad.bed| 259,691| -709,923 | Removal of multiple transcripts and duplicated rows|

## BC1958 Reference for frequently altered exons
file was generated after initial runnig of the CNV pipeline of the 1958BC cohort
Exons were selected from the output results from xhmm which occured at a frequency greater than 5%
across all 1000 samples
R (in cnvANNO - commented out):
```R
freq.list <- as.vector(unique(x$EXON[x$AF_all > 0.05]))
write.table(freq.list, file = "BC1958_freqentCNVs_5pct.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
```
