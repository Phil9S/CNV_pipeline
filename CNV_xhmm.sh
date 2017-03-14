#!/bin/bash

##Programme and reference files
gatk="/data/Resources/Software/GenomeAnalysisTK.jar"
ref="/data/Resources/hg38_2/hg38_2MergeAll.fa"
##Command-line variables
inputfolder="NULL"
#refbed="NULL"
int="NULL"
params="NULL"
timestamp=$(date +"%m_%Y")
outputfolder="NULL"
temp="FALSE"

##Admin BLOCK - Help, temp files, default message
#help
for arg in "$@"; do
	if [[ "$arg" == "-h" ]] || [[ "$arg" == "--help" ]]; then
		echo -e "
## CNV analysis ## - HELP Documentation - v2.0 - Buridan (XHMM ONLY) ##

Script for running xhmm CNV algorithms on a set of bam files

GATK and Reference genome files are internally specifed for Whisperwind (Server)
Please set these to alternate values by editing the script

Required arguments:

ARUGMENT			TYPE				DESCRIPTION

-i --input			folder				A directory path containing the .bam files to
								be used in the CNV analysis. Must end in .bam
								& have index files present in the same folder	
		
-o --output			folder				A directory path for a preferrably empty folder
								in which all resulting files will be deposited

-p --params (xhmm only)		file				Path for the params.txt file specifically requ-
								-ired for the xhmm analysis. Default params.txt
								file is included in the xhmm download. Only re-
								-quired if xhmm is being run.

-v --interval			file				Interval file provided by illumina, typically a
								file named Nextera_Exome_hg38, in tsv bed format,
								BED4 without header.
								
-b --bed			file				Bed file consisting of the interval positions of
								each gene. Used as a reference file for called
								CNVs and annotation of their proximal positions.

-t --temp			argument			Providing the option -t or --temp, given without
								a proceeding file/folder/string, etc. results in
								the retention of all temporary files generated
								during the CNV analysis. Mostly for DEBUGGING or
								QC and error checks.
Dependencies:

Will only run on Server: Whisperwind (without modification)

Packages:
- xhmm		
- GATK	
- bedtools

Core utilities:
- vim		- cat
- find		- sort
- sed

Examples:

Running xhmm analysis on a folder containing bams:
./CNV_analysis.sh -i /data/BAMS/ -o /data/CNV_REULTS/ -p /xhmm/files/params.txt -v /data/ref/nextera_exome_targets.bed -b /data/ref/hg38_genes.bed
		"			               
		echo -e "\n"
		exit
	fi
done
#default
if [[ $# -eq 0 ]]; then
	echo -e "\n## CNV Pipeline ## - You need to provide at least SOME arguments! Try using -h / --help for documentation and examples!\n"
	exit
fi
#temp file handling
for arg in "$@"; do
  if [[ "$arg" == "-t" ]] || [[ "$arg" == "--temp" ]]; then
    echo -e "\n## CNV Pipeline ## - Temporary files are being retained"
    temp="TRUE"
  fi
done


##arugement parsing block
while [[ $# > 1 ]]
	do 
	key="$1"
	case $key in
		-i|--input)
		inputfolder=$2
		shift
		;;
		-o|--output)
		outputfolder=$2
		shift
		;;
		-p|--params)
		params=$2
		shift
		;;
		-v|--interval)
		int=$2
		shift
		;;
#		-b|--bed)
#  	refbed=$2
#		shift
#		;;    
	esac
	shift
done


##Argument Checking Block
##INPUT
if [[ ! -d ${inputfolder} ]]; then
	echo -e "\n## CNV Pipeline ##\nERROR: The input folder provided does not exist - Please confirm the input folder path\nExiting Now"
	exit
fi

if [[ `ls ${inputfolder} | wc -l` < '1' ]]; then
	echo -e "## CNV Pipeline ##\nERROR: The input folder provided contains no files - Please confirm the input folder path\nExiting Now"
	exit
fi
##OUTPUT
if [[ ! -d ${outputfolder} ]]; then
	echo -e "\n## CNV Pipeline ##\nERROR: Output folder provided does not exist, provided using -o / --output - Please use -h / --help for documentation\nExiting Now"
	exit
fi

if [[ ! -w ${outputfolder} ]]; then
	echo -e "\n ## CNV Pipeline ##\nERROR: Output folder provided is not writable - Please make sure you have the correct permissions and are not writing in a protected folder\nExiting Now"
	exit
fi


##REF BED
#if [[ "$refbed" == "NULL" ]]; then
#	echo -e "\n## CNV Pipeline ##\nERROR: No reference bed file (gene intverals) specified, provided by -b / --bed - Please use -h / --help for documentation\nExiting Now"
#	exit
#fi

#if [[ ! -f ${refbed} ]]; then
#	echo -e "\n## CNV Pipeline ##\nERROR: Specifed reference bed file does not exist, provided by -b / --bed - Please use -h / --help for documentation\nExiting Now"
#	exit
#fi

##INTERVAL
if [[ "$int" == "NULL" ]]; then
	echo -e "\n## CNV Pipeline ##\nERROR: No interval file (exome/targeted panel etc) specified, provided by -v / --interval - Please use -h / --help for documentation\nExiting Now"
	exit
fi

if [[ ! -f ${int} ]]; then
	echo -e "\n## CNV Pipeline ##\nERROR: Specified  interval file (exome/targeted panel etc) does not exist, provided by -v / --interval - Please use -h / --help for documentation\nExiting Now"
	exit
fi

##Variable reporting



##Generating work-environment folder
if [ -d "${outputfolder}cnv_analysis" ]; then
	echo -e "\n## CNV Pipeline ## - Root folder exists - folder not generated\n"
else
	mkdir ${outputfolder}cnv_analysis
fi
cd ${outputfolder}cnv_analysis



###XHMM output process block##
echo -e "## XHMM ANALYSIS ## - Started at: $(date)\n"

###Output directory
if [ -d "xhmm_analysis_${timestamp}" ]; then
	echo -e "## XHMM ANALYSIS ## - Analysis folder exists - folder not generated\n"
else
	mkdir xhmm_analysis_${timestamp}
fi

###Temp folder generation
if [ -d "xhmm_analysis_${timestamp}/temp" ]; then
	echo -e "## XHMM ANALYSIS ## - Temporary folder exists - folder not generated\n"
else 
	mkdir xhmm_analysis_${timestamp}/temp
fi
cd xhmm_analysis_${timestamp}/temp
cp ${int} xhmm.intervals
vim -c "%s/\(\S\+\)\t\(\S\+\)\t\(\S\+\)/\1:\2-\3/g|wq" xhmm.intervals 
interval="xhmm.intervals"		
find ${inputfolder} -name *.bam -type f > bam_list_xhmm

###XHMM Analysis

echo -e "## XHMM ANALYSIS ## - Bam files split into 4 sets...(Stage 1 of 7)\n"
split -a 1 --numeric-suffixes=1 --additional-suffix=.list -n l/4 bam_list_xhmm bam_chunk

echo -e "## XHMM ANALYSIS ## - Performing depth of coverage...(Stage 2 of 7)\n"

java -Xmx30g -jar ${gatk} -T DepthOfCoverage -I bam_chunk1.list -L ${interval} -R ${ref} -dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable --minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 --includeRefNSites --countType COUNT_FRAGMENTS -o bam_chunkOUT1 &

java -Xmx30g -jar ${gatk} -T DepthOfCoverage -I bam_chunk2.list -L ${interval} -R ${ref} -dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable --minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 --includeRefNSites --countType COUNT_FRAGMENTS -o bam_chunkOUT2 &

java -Xmx30g -jar ${gatk} -T DepthOfCoverage -I bam_chunk3.list -L ${interval} -R ${ref} -dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable --minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 --includeRefNSites --countType COUNT_FRAGMENTS -o bam_chunkOUT3 &

java -Xmx30g -jar ${gatk} -T DepthOfCoverage -I bam_chunk4.list -L ${interval} -R ${ref} -dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable --minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 --includeRefNSites --countType COUNT_FRAGMENTS -o bam_chunkOUT4 &

###Allow for all child processes in parallel to complete
wait
sleep 5

echo -e "## XHMM ANALYSIS ## - Merging depth of coverage files & Calculating GC content...(Stage 3 of 7)\n"
	
###Combines GATK Depth-of-Coverage outputs for multiple samples (at same loci):
xhmm --mergeGATKdepths -o xhmmCNV.mergeDepths.txt --GATKdepths bam_chunkOUT1.sample_interval_summary --GATKdepths bam_chunkOUT2.sample_interval_summary --GATKdepths bam_chunkOUT3.sample_interval_summary --GATKdepths bam_chunkOUT4.sample_interval_summary > /dev/null 2>&1

###calculates the GC Content of the exome intervals
java -Xmx3072m -jar ${gatk} -T GCContentByInterval -L ${interval} -R ${ref} -o DATA_GC_percent.txt > /dev/null 2>&1
echo -e "## XHMM ANALYSIS ## - Removing extreme GC regions and centering to mean read depth...(Stage 4 of 7)\n"
	
###Concatonates and asseses GC content (if less than 0.1 or more than 0.9 -> print to new file
cat DATA_GC_percent.txt | awk '{if ($2 < 0.1 || $2 > 0.9) print $1}' > extreme_gc_targets.txt
###Centers the data about the mean and filters high/low GC intervals out of analysis
xhmm --matrix -r xhmmCNV.mergeDepths.txt --centerData --centerType target -o xhmmCNV.filtered_centered.RD.txt --outputExcludedTargets xhmmCNV.filtered_centered.RD.txt.filtered_targets.txt --outputExcludedSamples xhmmCNV.filtered_centered.RD.txt.filtered_samples.txt --excludeTargets extreme_gc_targets.txt --minTargetSize 10 --maxTargetSize 10000 --minMeanTargetRD 10 --maxMeanTargetRD 500 --minMeanSampleRD 25 --maxMeanSampleRD 200 --maxSdSampleRD 150 > /dev/null 2>&1

echo -e "## XHMM ANALYSIS ## - Analysing PCA plot & Normalising data...(Stage 5 of 7)\n"
###Performs PCA to generate component variation - decreases data variability due to 1st-nth priciple components
xhmm --PCA -r xhmmCNV.filtered_centered.RD.txt --PCAfiles xhmmCNV.mergeDepths_PCA > /dev/null 2>&1

###Normalises the mean centered data using the PCA data
xhmm --normalize -r xhmmCNV.filtered_centered.RD.txt --PCAfiles xhmmCNV.mergeDepths_PCA --normalizeOutput xhmmCNV.PCA_normalized.txt --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7 > /dev/null 2>&1

###Generates and asseses z-score distribution of mean centered-normalised read depth data and filters inappropriate intervals
xhmm --matrix -r xhmmCNV.PCA_normalized.txt --centerData --centerType sample --zScoreData -o xhmmCNV.PCA_normalized.filtered.sample_zscores.RD.txt --outputExcludedTargets xhmmCNV.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt --outputExcludedSamples xhmmCNV.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt --maxSdTargetRD 30 > /dev/null 2>&1

###applies the normalisation and z-scoring to the standard non-normalised and centered data set
xhmm --matrix -r xhmmCNV.mergeDepths.txt --excludeTargets xhmmCNV.filtered_centered.RD.txt.filtered_targets.txt --excludeTargets xhmmCNV.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt --excludeSamples xhmmCNV.filtered_centered.RD.txt.filtered_samples.txt --excludeSamples xhmmCNV.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt -o xhmmCNV.same_filtered.RD.txt > /dev/null 2>&1
	
###performs assessment of the z-score to identify high levels of statistcal deviation in interval regions
xhmm --discover -p ${params} -r xhmmCNV.PCA_normalized.filtered.sample_zscores.RD.txt -R xhmmCNV.same_filtered.RD.txt -c xhmmCNV.xcnv -a xhmmCNV.aux_xcnv -s xhmmCNV > /dev/null 2>&1

echo -e "## XHMM ANALYSIS ## - Genotyping called CNVs...(Stage 6 of 7)\n"

###genotypes indentified CNV during prior discovery steps
xhmm --genotype -p ${params} -r xhmmCNV.PCA_normalized.filtered.sample_zscores.RD.txt -R xhmmCNV.same_filtered.RD.txt -g xhmmCNV.xcnv -F ${ref} -v xhmmCNV.vcf > /dev/null 2>&1

###Results annotation & formatting##
echo -e "## XHMM ANALYSIS ## - Annotating & formatting output...(Stage 7 of 7)\n"

###xcnv to bed format conversion
if (( $(cat xhmmCNV.xcnv | wc -l) < '2' )); then
	echo -e "## XHMM ANALYSIS ## - ERROR: No CNVs called - Likely too few samples\n"
	echo -e "## XHMM ANALYSIS ## - XHMM analysis exiting"
	exit
fi 
#sed '1d' xhmmCNV.xcnv > xhmmCNV.bed
#vim -c "%s/\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\S\+\t\(\S\+\)\t\S\+\t\S\+\t\S\+\t\S\+\t\S\+\t\S\+\t\S\+\t\S\+\t\S\+\t\S\+/\1\t\2\t\3/g|wq" xhmmCNV.bed
###subsitition of xhmm xcnv file for removing extra columns
#vim -c "%s/\(\S\+\)\t\(\S\+\)\t\(\S\+\):\(\S\+\)-\(\S\+\)/chr\3\t\4\t\5\t\2_\1/g|wq" xhmmCNV.bed
###rearranging trimmed columns into bed4 format (chr, start, end, info)
#sort -k 1,1 -k 2,2 xhmmCNV.bed > xhmmCNV.sort.bed
#mv xhmmCNV.sort.bed xhmmCNV.bed
###Intersect of found cnvs with known ref gene intervals - one result per interval
#intersectBed -a xhmmCNV.bed -b ${refbed} -wb > xhmmCNV.anno.tsv
###Output formating sed commands and header addition
#vim -c "%s/\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)/\4\t\1\t\2\t\3\t\8\t\5\t\6\t\7/g|wq" xhmmCNV.anno.tsv
#vim -c "%s/\(DEL\)_\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)/\2\t\1\t\3\t\4\t\5\t\6\t\7\t\8\t\9/g|wq" xhmmCNV.anno.tsv
#vim -c "%s/\(DUP\)_\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)/\2\t\1\t\3\t\4\t\5\t\6\t\7\t\8\t\9/g|wq" xhmmCNV.anno.tsv
#echo -e "sample\tcnv_type\tchr\tstart\tend\tgene\tgene_chr\tgene_start\tgene_end" > xhmmheader.tsv
#cat xhmmCNV.anno.tsv >> xhmmheader.tsv
#mv xhmmheader.tsv ../xhmmCNV.anno.tsv
mv xhmmCNV.vcf ../xhmmCNV.vcf
mv xhmmCNV.xcnv ../xhmmCNV.xcnv
mv bam_list_xhmm ../xhmm_samplelist.txt
#if [[ "$temp" == "FALSE" ]]; then
#	cd ../
#	rm -r temp
#fi
echo -e "## XHMM ANALYSIS ## - COMPLETE!"
