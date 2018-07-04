# Generation of BED4 files & frequency lists from reference datasets

## Exon Interval file 
The interval file (specified by -v --interval) is a required file for the CNV Pipeline to run but is likely to be the cause of most common errors, run failures, and/or annotation failures. This file attempts to document how the interval file was generated and how a reference set of frequently altered genes/exons was created.
### Required files:
#### Reference intervals

Exon bed file was downloaded from [Biomart](https://www.ensembl.org/biomart) in tsv format and reformatted to fit BED4 specfications and the specifications detailed in the main user guide.

The database attributes gathered were (unique entries only): 
- "Chromosome/scaffold name"
- "Exon region start (bp)"
- "Exon region end (bp)"
- "Gene name"
- "Exon rank in transcript"
		
Reformatting steps (Reference intervals):
```sh
sed -i '1d' exons_hg19.txt
sed -i 's/\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)/\1\t\2\t\3\t\4_\5/g' exons_hg19.txt
sort -k1,1 -k2,2n exons_hg19.txt > exons_hg19.sorted.bed
```
#### Probes
Nextera probe positions were downloaded from [illumina](http://emea.support.illumina.com/sequencing/sequencing_kits/nextera-rapid-capture-exome-kit/downloads.html) to match library preparation kit

#### Repeatmask files
An important component of the interval file used is that is appropriately filtered for regions that are **not** of interest; this included regions overlapping low-complexity, repeat masked, regions of the genome. These regions are not useful and add noise to CNV calling so can be justifibly removed from the targets. The site repeatmasker.org provides catagorical fasta files for each type and span of these regions.

Fasta files containing the type and genomic position of repeat-masked regions can be downloaded from [repeatmasker.org](http://www.repeatmasker.org/) using the following commands:
```sh
wget -c LINK_TO_REPEATMASKER_FA_FILE
gunzip hg19.fa.out.gz
```
The repeat-mask fasta is not immediately appropriate for use so some minor pre-processing steps are required to allow for interval comparisons (i.e. delimiter alterations, repeat type selection, and sorting by chromosome and position), shown below:
```sh
awk -v OFS="\t" '$1=$1' hg19.fa.out > hg19_tab.fa.out
grep "Simple_repeat" hg19_tab.fa.out >> lowcomplex_simpreps.hg19.bed
grep "Low_complexity" hg19_tab.fa.out >> lowcomplex_simpreps.hg19.bed
cut -f5-7 lowcomplex_simpreps.hg19.bed > lowcomplex_simpreps.cut.hg19.bed
sed -i 's/chr//g' lowcomplex_simpreps.cut.hg19.bed #Skip if using >=GRCh38 
sort -k1,1 -k2,2n lowcomplex_simpreps.cut.hg19.bed > lowcomplex_simpreps.sorted.hg19.bed

```
### Bed file generation:
#### Generating list of overlapping intervals 
Filtering out intervals with any amount of repeat-mask overlap would be overly stringent so only exome intervals harbouring an overlap of 25% or more are excluded - In line with the steps performed [here](http://www.nature.com/ng/journal/v48/n10/full/ng.3638.html)

The command below compares the exome interval set to the repeat-masked interval set generating a list of exons that are overlapped by repeat-masked regions by > 25%:
```sh
bedtools intersect -wb -v -F 0.25 \
				   -a exons_hg19.sorted.bed \
				   -b lowcomplex_simpreps.sorted.hg19.bed > exons_hg19.masked.bed

sort -k1,1 -k2,2n exons_hg19.masked.bed > exons_hg19.masked.sorted.bed
```

#### BED target filtering
Regions overlapping with 50% of matching targets from the library preparation probes were retained using [Bedtools](http://bedtools.readthedocs.io/en/latest/)  intersect

```sh
bedtools intersect -loj -wb -F 0.5 \
				   -b exons_hg19.masked.sorted.bed \
				   -a nextera_exome_targets_hg19.bed > intersect.txt
```
The preceding files are loaded into R for filtering and target merging using the following code:
```R
## Load libraries
library(dplyr)
library(stringr)
options(scipen = 999)
## Load data and remove empty fields
bed <- read.table("intersect.txt", sep="\t",comment.char = "",quote = "",fill = T,stringsAsFactors = F)
bed <- bed[bed$V1 != ".",]
bed <- bed[bed$V4 != ".",]
bed <- bed[,-c(4:6)]
## Keep unique and remove misencoded exons
bed <- unique(bed)
bed <- bed[which(!grepl(bed$V7,perl = T,pattern = "\\S+_\\S+_\\S+")),]
## Collapse on genomic position
bed_col <- bed %>% group_by(V1,V2,V3) %>% summarise_all(funs(paste(unique(.),collapse = ",")))
## Add padding
bed_col$V2 <- bed_col$V2 - 10
bed_col$V3 <- bed_col$V3 + 10
## Write output
write.table(bed_col, "cnv_targets_masked_pad_hg19.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
```
##### Adding padding to exon regions - 10bp upstream of start and 10bp downstream of end
The final processing step is to add a number of padding bases to the intervals. This is intended to increase the fidelity of calls over target sets by incorporating more reads that align to the edges of target regions. 

The subsequent file is then sorted:
```sh
sort -k1,1 -k2,2n cnv_targets_masked_pad_hg19.bed > cnv_targets_masked_pad_sort_hg19.bed
```
Lastly the file is merged to collapse overlapping intervals into a single interval
```sh
bedtools merge -i cnv_targets_masked_pad_sort_hg19.bed -c 4 -o collapse > COLLAPSE.cnv_targets_masked_pad_sort_hg19.bed
# Renamed
mv COLLAPSE.cnv_targets_masked_pad_sort_hg19.bed cnv_targets_masked_pad_sort_hg19.bed
```
##### Removal summary:
An example file for CNV analysis: cnv_targets_masked_pad_sort_hg19.bed (7,953 Kb)

|File|Targets|Change|
|----|----|------|
|Original file| 811,594 | 0 | 
|Repeat masking| 761,902 | -49,692 |
|Library probe overlap| 555,392 | -206,510 | 
|Unique/Complete targets| 191,617  | -363,775 |
|Collapsed| 178,861  | -12,756 |

## Reference for frequently altered exons
The analysis provides frequencies for the occurance of exon alterations both internally (on the analysis cohort) & externally (currently the ICR1958 birth control cohort) to allow for both filtering and comparison of 'allelic' frequencies in the analysis

With some minor alterations to the scripts, a new freq_list can be generated on any reference dataset provided - This section will attempt to detail some of those steps

#### Generating unfiltered CNV data
In order to create a list of frequently altered genes/exons from a  reference dataset, the CNV analysis Pipeline should be run in its entirety with minor ammendments to the script `cnvANNO.R`

**_Follow the instructions for data generation outlined in both the `README.md` and the previous section of this document until you have generated a complete data set_**

Once that is complete, move into the folder containing the final analysis results (this should include the `cnv_xhmm_annotated.tsv` file and the `cnvANNO.RData` for your reference dataset

By editing the cnvANNO.R script which was copied to the working directory, it is possible to generate a new reference gene/exon frequency list through commenting and uncommenting lines. It is also possible to re-run the annotation with differing parameters.

##### Comment out the following lines:
###### Lines 40-41
```R
ref.list <- read.table("ref_CNVs.txt", sep = "\t", stringsAsFactors = FALSE)
colnames(ref.list) <- c("EXON","CNV","AF_ref")
```
###### Lines 122-126
```R
x <- merge(x, ref.list, by = c("EXON","CNV"), all.x = TRUE, fill = 0)
x[is.na(x)] <- 0
x <- cbind(x[1:11],x[ncol(x)],x[13:ncol(x)-1])
#remove commonly altered exons in BC1958 Cohort
x <- x[x$AF_ref < ref_af_value,]
```
###### Line 130
```R
x <- x[x$AF_all < int_af_value,]
```
###### Lines 145-146
```R
rm(aux,cnv,f.aux,gene_exon,intv,ref.list,AF_all,args,int_af_value,ref_af_value)
save.image(file="cnvANNO.RData")
```
##### Remove comments from the following lines:
###### Lines 115-116
```R
#freq.list <- x[,c("EXON","CNV","AF_all")]
#write.table(freq.list, file = "ref_CNVs.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
```
Save this modified `cnvANNO.R` script and run it again using the command below with the file path for the interval file used in the analysis:
```sh
Rscript ./cnvANNO.R /data/my_exome_intervals.bed
```
Once this is complete a new file, in this case `ref_CNVs.txt`, is generated. Copying this file into the git folder will replace the existing reference file with one matching your reference dataset of choice
