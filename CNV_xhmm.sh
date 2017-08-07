#!/bin/bash
date
## CNV Pipeline
##Programme and reference files
gatk="/data/Resources/Software/Javas/GenomeAnalysisTK.jar"
ref="/analysis/mtgroup_share/resources/gatk_bundle/b37/decompressed/human_g1k_v37.fasta"
##Command-line variables
version="Black Gold Radian"
inputfolder="NULL"
int="NULL"
params="NULL"
cohort="default"
outputfolder="NULL"
temp="TRUE"
minTargetSize=1
maxTargetSize=20000
minMeanTargetRD=10
maxMeanTargetRD=6000
minMeanSampleRD=30
maxMeanSampleRD=2000
maxSdSampleRD=500
PVE_mean_factor=0.7
maxSdTargetRD=100
call="FALSE"

#Admin BLOCK - Help, temp files, default message
#help
for arg in "$@"; do
	if [[ "$arg" == "-h" ]] || [[ "$arg" == "--help" ]]; then
		echo -e "
## CNV analysis ## - HELP Documentation - v.${version} ##

Script for running xhmm CNV algorithms on a set of bam files

Required arguments:

ARUGMENT			TYPE				DESCRIPTION

-i --input			folder				A directory path containing the .bam files to
								be used in the CNV analysis. Must end in .bam
								& have index files present in the same folder
	
-c --cohort			string				A string or identifying label to mark the ap-
								propraite output folder. e.g \"RCC\" for the
								kidney cancer dataset		

-o --output			folder				A directory path for a preferrably empty folder
								in which all resulting files will be deposited

-p --params (xhmm only)		file				Path for the params.txt file specifically requ-
								-ired for the xhmm analysis. Default params.txt
								file is included in the xhmm download. Only re-
								-quired if xhmm is being run.

-v --interval			file				Interval file provided by illumina, typically a
								file named Nextera_Exome_hg38, in tsv bed format,
								BED4 without header.

-tr --temp-remove		argument			Providing the option -tr or --temp-remove, given 
								without a proceeding file/folder/string results 
								in the DELETION of all temporary files generated
								during the CNV analysis, inlcuding Depth of cov-
								erage.

Optional Arguments:

ARGUMENT			TYPE				DESCRIPTION

-gatk --gatk			file				Location of the GenomeAnalysisTK.jar file used f-
								-or depth of coverage analysis - Default 

-r --ref			file				Reference fasta for the genome used in the analy-
								-sis. Reference must match contigs used in bed fi-
								-le provided - GATK will error if reference does
								not match the bed. Default reference provided is
															
-minTS --minTS			integer				Set minimum allowed target size for target filter-
								ing in xhmm --matrix step. Higher values are more
								stringent

-maxTS --maxTS			integer				Set maximum allowed target size for target filter-
                                                                ing in xhmm --matrix step. Lower values are more
                                                                stringent

-minMeanTRD --minMeanTRD	integer				Set minimum mean value read depth for each target
								in the interval file - targets with mean Read de-
								pth below this value are removed from the analysis

-maxMeanTRD --maxMeanTRD	integer				Set maximum mean value read depth for each target
                                                                in the interval file - targets with mean Read de-
                                                                pth above this value are removed from the analysis

-minMeanSRD --minMeanSRD        integer                         Set minimum mean value read depth for each sample
                                                                in the interval file - samples with mean Read de-
                                                                pth below this value are removed from the analysis

-maxMeanSRD --maxMeanSRD        integer                         Set maximum mean value read depth for each sample
                                                                in the interval file - samples with mean Read de-
                                                                pth above this value are removed from the analysis

-maxSampleSD --maxSampleSD	integer				Set the maximum allowed stardard deviation for e-
								ach sample - Samples are removed if their SD exc-
								eeds this value

-maxTargetSD --maxTargetSD	integer				Set the maximum allowed stardard deviation for e-
                                                                ach target - Targets are removed if their SD exc-
                                                                eeds this value
												
-pca --pca			float				PCA cut-off used for removing variance - Princip-
								le components exceeding this value as a percenta-
								ge are normalised and removed from the data to 
								reduce noise within the data

Dependencies:

Will only run on Server: Whisperwind (without modification)

Softwares:
- xhmm		
- GATK

- R:
Packages
- ggplot2	- dplyr
- tidyr		- data.table
- stringr	- ggrepel
- reshape2

Core utilities:
- vim		- cat
- find		- sort
- sed		- awk

Examples:

./CNV_analysis.sh [Required arguments] [optional arguments] [Mode tags]

Running xhmm analysis on a folder containing bams & remove temp files:
./CNV_analysis.sh -c My_project -i /data/BAMS/ -o /data/CNV_REULTS/ -p params.txt -v ref_exome.bed -tr
Running xhmm calling with new values:
./CNV_analysis.sh -c Existing_project -i /data/BAMS/ -o /data/CNV_REULTS/ -p params.txt -v ref_exome.bed -minTS 10 -maxTS 1000 -call
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
  if [[ "$arg" == "-tr" ]] || [[ "$arg" == "--temp-remove" ]]; then
    echo -e "\n## CNV Pipeline ## - Temporary files are being DELETED"
    temp="FALSE"
  fi
done

for arg in "$@"; do
  if [[ "$arg" == "-call" ]] || [[ "$arg" == "--call-only" ]]; then
    call="TRUE"
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
		-c|--cohort)
                cohort=$2
                shift
                ;;
		-gatk|--gatk)
                gatk=$2
                shift
                ;;
		-r|--ref)
                ref=$2
                shift
                ;;
		-minTS|--minTS)
                minTargetSize=$2
                shift
                ;;
		-maxTS|--maxTS)
                maxTargetSize=$2
                shift
                ;;
	 	-minMeanTRD|--minMeanTRD)
                minMeanTargetRD=$2
                shift
                ;;
		-maxMeanTRD|--maxMeanTRD)
                maxMeanTargetRD=$2
                shift
                ;;
		-minMeanSRD|--minMeanSRD)
                minMeanSampleRD=$2
                shift
                ;;
		-maxMeanSRD|--maxMeanSRD)
                maxMeanSampleRD=$2
                shift
                ;;
		-maxSampleSD|--maxSampleSD)
                maxSdSampleRD=$2
                shift
                ;;
		-maxTargetSD|--maxTargetSD)
               	maxSdTargetRD=$2 
                shift
                ;;
	        -pca|--pca)
                PVE_mean_factor=$2
                shift
                ;;
	esac
	shift
done


##Argument Checking Block

##Variable check
integer_check='^[1-9]$|^[1-9][0-9]+$'
float_check='^[0][.][0-9]+$|^[1][.][0]$'

## xhmm matrix and pca variables
if ! [[ ${minTargetSize} =~ ${integer_check} ]] ; then
        echo -e "\n## CNV Pipeline ##\nERROR: Non-integer supplied to -minTS / --minTS - Please use -h / --help for documentation\nExiting Now"
        exit
fi


if ! [[ ${maxTargetSize} =~ ${integer_check} ]] ; then
        echo -e "\n## CNV Pipeline ##\nERROR: Non-integer supplied to -maxTS / --maxTS - Please use -h / --help for documentation\nExiting Now"
        exit
fi

if ! [[ ${minMeanTargetRD} =~ ${integer_check} ]] ; then
        echo -e "\n## CNV Pipeline ##\nERROR: Non-integer supplied to -minMeanTRD / --minMeanTRD - Please use -h / --help for documentation\nExiting Now"
        exit
fi

if ! [[ ${maxMeanTargetRD} =~ ${integer_check} ]] ; then
        echo -e "\n## CNV Pipeline ##\nERROR: Non-integer supplied to -maxMeanTRD / --maxMeanTRD - Please use -h / --help for documentation\nExiting Now"
        exit
fi

if ! [[ ${minMeanSampleRD} =~ ${integer_check} ]] ; then
        echo -e "\n## CNV Pipeline ##\nERROR: Non-integer supplied to -minMeanSRD / --minMeanSRD - Please use -h / --help for documentation\nExiting Now"
        exit
fi

if ! [[ ${maxMeanSampleRD} =~ ${integer_check} ]] ; then
        echo -e "\n## CNV Pipeline ##\nERROR: Non-integer supplied to -maxMeanSRD / --maxMeanSRD - Please use -h / --help for documentation\nExiting Now"
        exit
fi

if ! [[ ${maxSdSampleRD} =~ ${integer_check} ]] ; then
        echo -e "\n## CNV Pipeline ##\nERROR: Non-integer supplied to -maxSampleSD / --maxSampleSD - Please use -h / --help for documentation\nExiting Now"
        exit
fi

if ! [[ ${PVE_mean_factor} =~ ${float_check} ]] ; then
        echo -e "\n## CNV Pipeline ##\nERROR: Non-float supplied to -pca / --pca - Please use -h / --help for documentation\nExiting Now"
        exit
fi

if ! [[ ${maxSdTargetRD} =~ ${integer_check} ]] ; then
        echo -e "\n## CNV Pipeline ##\nERROR: Non-integer supplied to -maxTargetSD / --maxTargetSD - Please use -h / --help for documentation\nExiting Now"
        exit
fi



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

##INTERVAL
if [[ "$int" == "NULL" ]]; then
	echo -e "\n## CNV Pipeline ##\nERROR: No interval file (exome/targeted panel etc) specified, provided by -v / --interval - Please use -h / --help for documentation\nExiting Now"
	exit
fi

if [[ ! -f ${int} ]]; then
	echo -e "\n## CNV Pipeline ##\nERROR: Specified interval file (exome/targeted panel etc) does not exist, provided by -v / --interval - Please use -h / --help for documentation\nExiting Now"
	exit
fi

##GATK
if [[ ! -f ${gatk} ]]; then
        echo -e "\n## CNV Pipeline ##\nERROR: Specified GATK jar file does not exist, provided by -gatk / --gatk - Please use -h / --help for documentation\nExiting Now"
        exit
fi

##ref
if [[ ! -f ${ref} ]]; then
        echo -e "\n## CNV Pipeline ##\nERROR: Specified reference fasta does not exist, provided by -r / --ref - Please use -h / --help for documentation\nExiting Now"
        exit
fi

##Generating work-environment folder
if [ -d "${outputfolder}cnv_analysis" ]; then
	echo -e "\n## CNV Pipeline ## - Root folder exists - folder not generated\n"
else
	mkdir ${outputfolder}cnv_analysis
fi
cp cnvPCA.R ${outputfolder}cnv_analysis/
cp cnvANNO.R ${outputfolder}cnv_analysis/
cp cnvPLOTS.R ${outputfolder}cnv_analysis/
cp BC1958_CNVs.txt ${outputfolder}cnv_analysis/
cd ${outputfolder}cnv_analysis


###XHMM output process block##
echo -e "## XHMM ANALYSIS ## - Started at: $(date)\n"

###Output directory
if [ -d "xhmm_analysis_${cohort}" ]; then
	echo -e "## XHMM ANALYSIS ## - Analysis folder exists - folder not generated\n"
else
	mkdir xhmm_analysis_${cohort}
fi

###Temp folder generation
if [ -d "xhmm_analysis_${cohort}/temp" ]; then
	echo -e "## XHMM ANALYSIS ## - Temporary folder exists - folder not generated\n"
else 
	mkdir xhmm_analysis_${cohort}/temp
fi
mv cnvPCA.R xhmm_analysis_${cohort}/temp
mv cnvANNO.R xhmm_analysis_${cohort}
mv cnvPLOTS.R xhmm_analysis_${cohort}
mv BC1958_CNVs.txt xhmm_analysis_${cohort}
cd xhmm_analysis_${cohort}/temp

##Variable reporting
echo -e "## CNV Pipeline ## Argument summary:"
echo -2 "## CNV Pipeline ## Version ${version}"
echo -e "## CNV Pipeline ## Genome analysis TK jar file - ${gatk}" | tee -a cnv.log
echo -e "## CNV Pipeline ## Reference fasta file - ${ref}" | tee -a cnv.log
echo -e "## CNV Pipeline ## Input folder - ${inputfolder}" | tee -a cnv.log
echo -e "## CNV Pipeline ## Interval file - ${int}" | tee -a cnv.log
echo -e "## CNV Pipeline ## Parameters file - ${params}" | tee -a cnv.log
echo -e "## CNV Pipeline ## Cohort/Project name - ${cohort}" | tee -a cnv.log
echo -e "## CNV Pipeline ## Output folder - ${outputfolder}cnv_analysis/" | tee -a cnv.log
echo -e "## CNV Pipeline ## Temporary file deletion = ${temp}" | tee -a cnv.log
echo -e "## CNV Pipeline ## Minimum target size = ${minTargetSize}" | tee -a cnv.log
echo -e "## CNV Pipeline ## Maximum target size = ${maxTargetSize}" | tee -a cnv.log
echo -e "## CNV Pipeline ## Minimum mean target read depth = ${minMeanTargetRD}" | tee -a cnv.log
echo -e "## CNV Pipeline ## maximum mean target read depth = ${maxMeanTargetRD}" | tee -a cnv.log
echo -e "## CNV Pipeline ## Minimum mean sample read depth = ${minMeanSampleRD}" | tee -a cnv.log
echo -e "## CNV Pipeline ## Maximum mean sample read depth = ${maxMeanSampleRD}" | tee -a cnv.log
echo -e "## CNV Pipeline ## Maximum standard dev. sample read depth = ${maxSdSampleRD}" | tee -a cnv.log
echo -e "## CNV Pipeline ## Maximum standard dev. target read depth = ${maxSdTargetRD}" | tee -a cnv.log
echo -e "## CNV Pipeline ## PCA variance cut off set to ${PVE_mean_factor}" | tee -a cnv.log
echo -e "## CNV Pipeline ## Call only mode set to ${call}" | tee -a cnv.log

## folder setup
cp ${int} xhmm.intervals
vim -c "%s/\(\S\+\)\t\(\S\+\)\t\(\S\+\)\t\(\S\+\)/\1:\2-\3/g|wq" xhmm.intervals
interval="xhmm.intervals"
ls ${inputfolder}*.bam > bam_list_xhmm

#sleep 120
###XHMM Analysis
if [[ ${call} = "FALSE"  ]]; then
	echo "inside if"
	exit
	echo -e "## XHMM ANALYSIS ## - Bam files split into 6 sets...(Stage 1 of 10)\n"
	split -a 1 --numeric-suffixes=1 --additional-suffix=.list -n l/6 bam_list_xhmm bam_chunk

	echo -e "## XHMM ANALYSIS ## - Performing depth of coverage...(Stage 2 of 10)\n"

	java -Xmx30g -jar ${gatk} -T DepthOfCoverage -I bam_chunk1.list -L ${interval} -R ${ref} -dt BY_SAMPLE -dcov 5000 \
	-l INFO \
	--omitDepthOutputAtEachBase \
	--omitLocusTable \
	--minBaseQuality 0 \
	--minMappingQuality 20 \
	--start 1 \
	--stop 5000 \
	--nBins 200 \
	--includeRefNSites \
	--countType COUNT_FRAGMENTS \
	-o bam_chunkOUT1 > /dev/null 2>&1 &
	
	java -Xmx30g -jar ${gatk} -T DepthOfCoverage -I bam_chunk2.list -L ${interval} -R ${ref} -dt BY_SAMPLE -dcov 5000 \
	-l INFO \
	--omitDepthOutputAtEachBase \
	--omitLocusTable \
	--minBaseQuality 0 \
	--minMappingQuality 20 \
	--start 1 \
	--stop 5000 \
	--nBins 200 \
	--includeRefNSites \
	--countType COUNT_FRAGMENTS \
	-o bam_chunkOUT2 > /dev/null 2>&1 &
	
	java -Xmx30g -jar ${gatk} -T DepthOfCoverage -I bam_chunk3.list -L ${interval} -R ${ref} -dt BY_SAMPLE -dcov 5000 \
	-l INFO \
	--omitDepthOutputAtEachBase \
	--omitLocusTable \
	--minBaseQuality 0 \
	--minMappingQuality 20 \
	--start 1 \
	--stop 5000 \
	--nBins 200 \
	--includeRefNSites \
	--countType COUNT_FRAGMENTS \
	-o bam_chunkOUT3 &
	
	java -Xmx30g -jar ${gatk} -T DepthOfCoverage -I bam_chunk4.list -L ${interval} -R ${ref} -dt BY_SAMPLE -dcov 5000 \
	-l INFO \
	--omitDepthOutputAtEachBase \
	--omitLocusTable \
	--minBaseQuality 0 \
	--minMappingQuality 20 \
	--start 1 \
	--stop 5000 \
	--nBins 200 \
	--includeRefNSites \
	--countType COUNT_FRAGMENTS \
	-o bam_chunkOUT4 > /dev/null 2>&1 &
	
	java -Xmx30g -jar ${gatk} -T DepthOfCoverage -I bam_chunk5.list -L ${interval} -R ${ref} -dt BY_SAMPLE -dcov 5000 \
	-l INFO \
	--omitDepthOutputAtEachBase \
	--omitLocusTable \
	--minBaseQuality 0 \
	--minMappingQuality 20 \
	--start 1 \
	--stop 5000 \
	--nBins 200 \
	--includeRefNSites \
	--countType COUNT_FRAGMENTS \
	-o bam_chunkOUT5 > /dev/null 2>&1 &
	
	java -Xmx30g -jar ${gatk} -T DepthOfCoverage -I bam_chunk6.list -L ${interval} -R ${ref} -dt BY_SAMPLE -dcov 5000 \
	-l INFO \
	--omitDepthOutputAtEachBase \
	--omitLocusTable \
	--minBaseQuality 0 \
	--minMappingQuality 20 \
	--start 1 \
	--stop 5000 \
	--nBins 200 \
	--includeRefNSites \
	--countType COUNT_FRAGMENTS \
	-o bam_chunkOUT6 > /dev/null 2>&1 &
	

	###Allow for all child processes in parallel to complete
	wait
	sleep 5
	
	echo -e "## XHMM ANALYSIS ## - Merging depth of coverage files & Calculating GC content...(Stage 3 of 10)\n"
		
	###Combines GATK Depth-of-Coverage outputs for multiple samples (at same loci):
	xhmm --mergeGATKdepths -o xhmmCNV.mergeDepths.txt \
	--GATKdepths bam_chunkOUT1.sample_interval_summary \
	--GATKdepths bam_chunkOUT2.sample_interval_summary \
	--GATKdepths bam_chunkOUT3.sample_interval_summary \
	--GATKdepths bam_chunkOUT4.sample_interval_summary \
	--GATKdepths bam_chunkOUT5.sample_interval_summary \
	--GATKdepths bam_chunkOUT6.sample_interval_summary > /dev/null 2>&1
	
	else 
	echo -e "## XHMM ANALYSIS ## - Call and annotate ONLY mode - Depth of coverage and merge steps skipped (Steps 1-3)"
fi

if [[ ${call} == "TRUE" ]] && [[ `ls xhmmCNV.mergeDepths.txt | wc -l` < 1 ]]; then
	echo -e "## CNV Pipeline ## No merged depth of coverage file found whilst using call-only mode - You may need to re-run DepthofCoverage"
	exit
fi

exit
###calculates the GC Content of the exome intervals
java -Xmx30g -jar ${gatk} -T GCContentByInterval -L ${interval} -R ${ref} -o DATA_GC_percent.txt > /dev/null 2>&1
echo -e "## XHMM ANALYSIS ## - Removing extreme GC regions and centering to mean read depth...(Stage 4 of 10)\n"
	
###Concatonates and asseses GC content (if less than 0.1 or more than 0.9 -> print to new file
cat DATA_GC_percent.txt | awk '{if ($2 < 0.1 || $2 > 0.9) print $1}' > extreme_gc_targets.txt
###Centers the data about the mean and filters high/low GC intervals out of analysis
xhmm --matrix -r xhmmCNV.mergeDepths.txt --centerData --centerType target -o xhmmCNV.filtered_centered.RD.txt --outputExcludedTargets xhmmCNV.filtered_centered.RD.txt.filtered_targets.txt --outputExcludedSamples xhmmCNV.filtered_centered.RD.txt.filtered_samples.txt --excludeTargets extreme_gc_targets.txt --minTargetSize ${minTargetSize} --maxTargetSize ${maxTargetSize} --minMeanTargetRD ${minMeanTargetRD} --maxMeanTargetRD ${maxMeanTargetRD} --minMeanSampleRD ${minMeanSampleRD} --maxMeanSampleRD ${maxMeanSampleRD} --maxSdSampleRD ${maxSdSampleRD} > /dev/null 2>&1

##version of line before command line control was added - not used
#xhmm --matrix -r xhmmCNV.mergeDepths.txt --centerData --centerType target -o xhmmCNV.filtered_centered.RD.txt --outputExcludedTargets xhmmCNV.filtered_centered.RD.txt.filtered_targets.txt --outputExcludedSamples xhmmCNV.filtered_centered.RD.txt.filtered_samples.txt --excludeTargets extreme_gc_targets.txt --minTargetSize 1 --maxTargetSize 20000 --minMeanTargetRD 10 --maxMeanTargetRD 6000 --minMeanSampleRD 30 --maxMeanSampleRD 2000 --maxSdSampleRD 500 > /dev/null 2>&1

echo -e "## XHMM ANALYSIS ## - Analysing PCA plot...(Stage 5 of 10)\n"
###Performs PCA to generate component variation - decreases data variability due to 1st-nth priciple components
xhmm --PCA -r xhmmCNV.filtered_centered.RD.txt --PCAfiles xhmmCNV.mergeDepths_PCA > /dev/null 2>&1

wd=`pwd`
Rscript cnvPCA.R ${wd}

vim -c '%s/\(Variance threshold at Principle component: \)\(\S\+\)/\r\1\2/|wq' PCA_summary.txt
vim -c '%s/ Standard deviation/\rStandard deviation/g|wq' PCA_summary.txt
vim -c '%s/ Proportion of Variance/\rProportion of Variance/g|wq' PCA_summary.txt
vim -c '%s/ Cumulative Proportion/\rCumulative Proportion/g|wq' PCA_summary.txt
vim -c '%s/\s\+\(\PC[0-9]\+\)\s\+\(\PC[0-9]\+\)\s\+\(\PC[0-9]\+\)\s\+\(\PC[0-9]\+\)\s\+\(\PC[0-9]\+\)/\rComponent\t\1\t\2\t\3\t\4\t\5/g|wq' PCA_summary.txt
vim -c '%s/\(Standard deviation\)\s\+\(\S\+\)\s\+\(\S\+\)\s\+\(\S\+\)\s\+\(\S\+\)\s\+\(\S\+\)/\1\t\2\t\3\t\4\t\5\t\6/g|wq' PCA_summary.txt
vim -c '%s/\(Proportion of Variance\)\s\+\(\S\+\)\s\+\(\S\+\)\s\+\(\S\+\)\s\+\(\S\+\)\s\+\(\S\+\)/\1\t\2\t\3\t\4\t\5\t\6/g|wq' PCA_summary.txt
vim -c '%s/\(Cumulative Proportion\)\s\+\(\S\+\)\s\+\(\S\+\)\s\+\(\S\+\)\s\+\(\S\+\)\s\+\(\S\+\)/\1\t\2\t\3\t\4\t\5\t\6/g|wq' PCA_summary.txt
vim -c '1d|wq' PCA_summary.txt
vim -c '$m 0|wq' PCA_summary.txt

echo -e "## XHMM ANALYSIS ## - Normalising RD data on PCA results...(Stage 6 of 10)\n"
###Normalises the mean centered data using the PCA data
xhmm --normalize -r xhmmCNV.filtered_centered.RD.txt --PCAfiles xhmmCNV.mergeDepths_PCA --normalizeOutput xhmmCNV.PCA_normalized.txt --PCnormalizeMethod PVE_mean --PVE_mean_factor ${PVE_mean_factor} > /dev/null 2>&1

###Generates and asseses z-score distribution of mean centered-normalised read depth data and filters inappropriate intervals
xhmm --matrix -r xhmmCNV.PCA_normalized.txt --centerData --centerType sample --zScoreData -o xhmmCNV.PCA_normalized.filtered.sample_zscores.RD.txt --outputExcludedTargets xhmmCNV.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt --outputExcludedSamples xhmmCNV.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt --maxSdTargetRD ${maxSdTargetRD} > /dev/null 2>&1

###applies the normalisation and z-scoring to the standard non-normalised and centered data set
xhmm --matrix -r xhmmCNV.mergeDepths.txt --excludeTargets xhmmCNV.filtered_centered.RD.txt.filtered_targets.txt --excludeTargets xhmmCNV.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt --excludeSamples xhmmCNV.filtered_centered.RD.txt.filtered_samples.txt --excludeSamples xhmmCNV.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt -o xhmmCNV.same_filtered.RD.txt > /dev/null 2>&1

echo -e "## XHMM ANALYSIS ## - Discovering CNVs...(Stage 7 of 10)\n"	
###performs assessment of the z-score to identify high levels of statistcal deviation in interval regions
xhmm --discover -p ${params} -r xhmmCNV.PCA_normalized.filtered.sample_zscores.RD.txt -R xhmmCNV.same_filtered.RD.txt -c xhmmCNV.xcnv -a xhmmCNV.aux_xcnv -s xhmmCNV > /dev/null 2>&1

echo -e "## XHMM ANALYSIS ## - Genotyping called CNVs...(Stage 7 of 10)\n"

###genotypes indentified CNV during prior discovery steps
xhmm --genotype -p ${params} -r xhmmCNV.PCA_normalized.filtered.sample_zscores.RD.txt -R xhmmCNV.same_filtered.RD.txt -g xhmmCNV.xcnv -F ${ref} -v xhmmCNV.vcf > /dev/null 2>&1

###Results annotation & formatting##
echo -e "## XHMM ANALYSIS ## - Moving output from TEMP...(Stage 8 of 10)\n"

###xcnv to bed format conversion
if (( $(cat xhmmCNV.xcnv | wc -l) < '2' )); then
        echo -e "## XHMM ANALYSIS ## - ERROR: No CNVs called - Likely too few samples\n"
        echo -e "## XHMM ANALYSIS ## - XHMM analysis exiting"
        exit
fi

mv xhmmCNV.xcnv ../xhmmCNV.xcnv
mv bam_list_xhmm ../xhmm_samplelist.txt
mv PCA_Scree.png ../PCA_Scree.png
mv PCA_summary.txt ../PCA_summary.txt
mv xhmmCNV.aux_xcnv ../xhmmCNV.aux_xcnv
mv cnv.log ../cnv.log

cd ../

echo -e "## XHMM ANALYSIS ## - Annotating output...(Stage 9 of 10)\n"

Rscript cnvANNO.R ${int} > /dev/null 2>&1

echo -e "## XHMM ANALYSIS ## - Plotting Graphs...(Stage 10 of 10)\n"
Rscript cnvPLOTS.R > /dev/null 2>&1

mkdir cnv_plots
mv cnv_calls_*.png cnv_plots/

if [[ "$temp" == "FALSE" ]]; then
              rm -r temp
fi

echo -e "## XHMM ANALYSIS ## - COMPLETE!"
date

