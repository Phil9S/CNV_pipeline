# Generation of BED4 files & frequency lists from reference datasets

## Exon Interval file 
The interval file (specified by -v --interval) is a required file for the CNV Pipeline to run but is likely to be the cause of most common errors, run failures, and/or annotation failures. This file attempts to document how the interval file was generated and how a reference set of frequently altered genes/exons was created

#### Reference interval retrieval   

Whisperwind:
/data/Resources/Beds/All_exons_hg38.bed
/data/Resources/Beds/nexterarapidcapture_exome_targetedregions_v1.2_hg38.bed

[Bedtools](http://bedtools.readthedocs.io/en/latest/) provides an excellent tool set for easily comparing sets of genomic intervals for overlap and/or intersections
```sh
bedtools intersect -loj -wb -F 0.75 -b All_Exons_hg38.bed -a nextera_exome_targets_hg38.bed > intersect.txt
```

```R
intersect <- read.table("intersect.txt")
intersect.dedup <- intersect[!duplicated(intersect),]
intersect.blank.dedup <- intersect.dedup[!intersect.dedup$V4 == ".",]
intersect.blank.dedup <- intersect.blank.dedup[!duplicated(intersect.blank.dedup[c(1:3,7)]),]
intersect.blank.dedup.out <- intersect.blank.dedup[!duplicated(intersect.blank.dedup[c(4:6)]),]
write.table(intersect.blank.dedup.out[c(1:3,7)],"cnv_nextera_targets_hg38.bed",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
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

The command below compares the exome interval set to the repeat-masked interval set generating a list of exons that are overlapped by repeat-masked regions by > 50%:
```sh
bedtools intersect -wb -F 0.5 /
		   -a cnv_nextera_targets_hg38.bed /
		   -b lowcomplex_simpreps.sorted.hg38.bed > exon_repmask_overlap.txt
vim -c '%s/\S\+\t\S\+\t\S\+\t\(\S\+\)\t\S\+\t\S\+\t\S\+\t\S\+/\1/|wq' exon_repmask_overlap.txt
sort -u exon_repmask_overlap.txt > exon_repmask_overlap.uniq.txt
```
This list is then loaded into R with the exome interval BED file and non-matching exons are retained:
```R
t <- read.table("cnv_nextera_targets_hg38.bed", sep="\t")
l <- read.table("exon_repmask_overlap.uniq.txt")
t1 <- t[!t$V4 %in% l$V1,]
write.table(t1, "cnv_nextera_targets_masked_hg38.bed", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
```
#### Adding padding to exon regions - 2bp upstream of start and 2bp downstream of end
The final processing step is to add a number of padding bases to the intervals. This is intended to increase the fidelity of calls over target sets by incorporating more reads that align to the edges of target regions. Using R, code below adds a padding of 2bp to both the start and stop position of each exon present in the interval file:

```R
t <- read.table("cnv_nextera_targets_masked_hg38.bed", sep="\t")
n <- t
n$V2 <- n$V2 - 2
n$V3 <- n$V3 + 2
write.table(n, "cnv_nextera_targets_masked_pad_hg38.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
```
The subsequent file is then sorted and moved to have the desired name via the commandline:
```sh
sort -k1,1 -k2,2n cnv_nextera_targets_masked_pad_hg38.bed > cnv_nextera_targets_masked_pad_hg38.bedsort
mv cnv_nextera_targets_masked_pad_hg38.bedsort cnv_nextera_targets_masked_pad_hg38.bed
```

##### Removal summary:
Final interval file for CNV analysis: hg38_exons_Ensembl87.masked.autoXY.sorted.dedup.pad.bed (8,875Kb)

|File|Size|Change|Process|
|----|----|------|-------|
|All_Exons_hg38.bed| 1,271,759 targets | Original file - Including all transcripts |
|intersect.txt| 859,951 targets| -411,808 | Exon intersecting targets|
|cnv_nextera_targets_hg38.bed| 244,380 targets| -615,571 | Removal of duplicated and overlapping targets|
|cnv_nextera_targets_masked_pad_hg38.bed| 223,527 targets| -20853 | Repeat masked targets|

## Reference for frequently altered exons
The analysis provides frequencies for the occurance of exon alterations both internally (on the analysis cohort) & externally (currently the ICR1958 birth control cohort) to allow for both filtering and comparison of 'allelic' frequencies in the analysis

With some minor alterations to the scripts, a new freq_list can be generated on any reference dataset provided - This section will attempt to detail some of those steps

#### Generating unfiltered CNV data
In order to create a list of frequently altered genes/exons from a  reference dataset, the CNV analysis Pipeline should be run in its entirety with minor ammendments to the script `cnvANNO.R`

**_Follow the instructions for data generation outlined in both the `README.md` and the previous section of this document until you have generated a complete data set_**

Once that is complete, move into the folder containing the final analysis results (this should include the `cnv_xhmm_annotated.tsv` file and the `cnvANNO.RData` for your reference dataset

By editing the cnvANNO.R script which was copied to the working directory, it is possible to generate a new reference gene/exon frequency list through commenting and uncommenting lines. It is also possible to re-run the annotation with differing parameters.

##### Comment out the following lines:
###### Lines 42-43
```R
ref.list <- read.table("ref_CNVs.txt", sep = "\t", stringsAsFactors = FALSE)
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
