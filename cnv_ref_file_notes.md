# Generation of BED4 files & frequency lists from reference datasets

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

[Bedtools](http://bedtools.readthedocs.io/en/latest/) provides an excellent tool set for easily comparing sets of genomic intervals for overlap and/or intersections - the command below compares the exome interval set to the repeat-masked interval set generating a list of exons that are overlapped by repeat-masked regions by > 25%:
```sh
bedtools intersect -wb -F 0.25 /
		   -a hg38_exons_Ensembl87.autoXY.sorted.bed /
		   -b lowcomplex_simpreps.sorted.hg38.bed > exon_repmask_overlap.txt
vim -c '%s/\S\+\t\S\+\t\S\+\t\(\S\+\)\t\S\+\t\S\+\t\S\+\t\S\+/\1/' exon_repmask_overlap.txt
sort -u exon_repmask_overlap.txt > exon_repmask_overlap.uniq.txt
```
This list is then loaded into R with the exome interval BED file and non-matching exons are retained:
```R
t <- read.table("hg38_exons_Ensembl87.autoXY.sorted.bed", sep="\t")
l <- read.table("exon_repmask_overlap.uniq.txt")
t1 <- t[!t$V4 %in% l$V1,]
write.table(t1, "hg38_exons_Ensembl87.masked.autoXY.sorted.bed", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
```
#### Removal of duplicated exons 
###### NOTE: THIS STEP TAKES SOME TIME! SCALE OF 12+ HOURS NOT MINUTES!
###### _This will retain duplicate exons of the same size - differing start and stop positions are present for a majority of these_
###### _This step is optional - Downstream scripts have built-in error catching to deal with duplications_

Due to multiple transcripts, it is highly likely that many exons in the interval file will have been duplicated, all be with slightly differing start and stop positions. Instead of systematically identifying the canonical transcript for each gene, the most simple solution is to select the largest exon, to encompass the most genomic space

The following R script reads in the interval file and retains the gene_exon line with the largest size, rejecting all other duplicates (It will retain duplicates where both the start and stop positions differ entirely, as these are likely unique exons for a specific transcript):
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
_This is a long process but dramatically reduces the number of targets, resulting in a substantial downstream decrease in processing time at all stages_

#### Adding padding to exon regions - 2bp upstream of start and 2bp downstream of end
The final processing step is to add a number of padding bases to the intervals. This is intended to increase the fidelity of calls over target sets by incorporating more reads that align to the edges of target regions. Using R, code below adds a padding of 2bp to both the start and stop position of each exon present in the interval file:

```R
t <- read.table("hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.bed", sep="\t")
n <- t
n$V2 <- n$V2 - 2
n$V3 <- n$V3 + 2
write.table(n, "hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.pad.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
```
The subsequent file is then sorted and moved to have the desired name via the commandline:
```sh
sort -k1,1 -k2,2n hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.pad.bed > hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.pad.bedsort
mv hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.pad.bedsort hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.pad.bed
```

##### Removal summary:
Final interval file for CNV analysis: hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.pad.bed (8,875Kb)

|File|Size|Change|Process|
|----|----|------|-------|
|hg38_exons_Ensembl87.sorted.txt| 1,285,529 targets | Original file |
|hg38_exons_Ensembl87.autoXY.sorted.bed| 1,181,871 targets| -103,658 | Removal of non-standard chromosomes|
hg38_exons_Ensembl87.masked.autoXY.sorted.bed| 969,614 targets| -212,257 | Removal of masked repeat regions|
hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.pad.bed| 259,691| -709,923 | Removal of multiple transcripts and duplicated rows|

## BC1958 Reference for frequently altered exons
The analysis provides frequencies for the occurance of exon alterations both internally (on the analysis cohort) & externally (currently the ICR1958 birth control cohort) to allow for both filtering and comparison of 'allelic'frequencies in the analysis

With some minor alterations to the scripts, a new freq_list can be generated on any reference dataset provided - This section will attempt to detail some of those steps

#### Generating unfiltered CNV data
In order to create a list of frequently altered genes/exons from a  reference dataset, the CNV analysis Pipeline should be run in its entirety with minor ammendments to the script `cnvANNO.R`

**_Follow the instructions for data generation outlined in both the `README.md` and the previous section of this document until you have generated a complete data set_**

Once that is complete, move into the folder containing the final analysis results (this should include the `cnv_xhmm_annotated.tsv` file and the `cnvANNO.RData` for your reference dataset

By editing the cnvANNO.R script which was copied to the working directory, it is possible to generate a new reference gene/exon frequency list through commenting and uncommenting lines. It is also possible to re-run the annotation with differing parameters.

##### Comment out the following lines:
###### Lines 42-43
```R
ref.list <- read.table("BC1958_CNVs.txt", sep = "\t", stringsAsFactors = FALSE)
colnames(ref.list) <- c("EXON","CNV","AF_ref")
```
###### Lines 124-128
```R
x <- merge(x, ref.list, by = c("EXON","CNV"), all.x = TRUE, fill = 0)
x[is.na(x)] <- 0
x <- cbind(x[1:11],x[ncol(x)],x[13:ncol(x)-1])
#remove commonly altered exons in BC1958 Cohort
x <- x[x$AF_ref < ref_af_value,]
```
###### Line 132
```R
x <- x[x$AF_all < int_af_value,]
```
###### Lines 147-148
```R
rm(aux,cnv,f.aux,gene_exon,intv,ref.list,AF_all,args,int_af_value,ref_af_value)
save.image(file="cnvANNO.RData")
```
##### Remove comments from the following lines:
###### Lines 117-118
```R
#freq.list <- x[,c("EXON","CNV","AF_all")]
#write.table(freq.list, file = "ref_CNVs.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
```
Save this modified `cnvANNO.R` script and run it again using the command below with the file path for the interval file used in the analysis:
```sh
Rscript ./cnvANNO.R /data/my_exome_intervals.bed
```
Once this is complete a new file, in this case `ref_CNVs.txt`, is generated. Copying this file into the git folder will replace the existing reference file with one matching your reference dataset of choice
