# CNV Pipeline (Black Gold Radian)
###### Script for running xhmm CNV algorithms on a set of bam files using xhmm - [Fromer et al. 2014](http://onlinelibrary.wiley.com/doi/10.1002/0471142905.hg0723s81/abstract;jsessionid=FDFFC8DB52043AE0D6DDE0DEEF2C9F9B.f03t02)
###### _NB: This will only work on Whisperwind (Dept of Medical Genetics) and will require major revisions on other machines_

##### Quick Start:
In order to run a basic analysis with the default arguments run the code below with modified variables for the data of interest:

###### _Check that all [dependecies](https://github.com/Phil9S/CNV_pipeline#dependencies) are installed then:_
######  _clone the git_
```
git clone https://github.com/Phil9S/CNV_pipeline.git
```
###### _move to the git directory & run the following command_
```Bash
cd CNV_pipeline/
./CNV_analysis.sh -c My_project -i /data/BAMS/ -o /data/CNV_REULTS/ -p params.txt -v ref_exome.bed
```
_Remember to correctly specify the params.txt, interval file, and change the ref (-r) to one matching the interval file_
##### _Known issues_ 
* ###### Failed runs can occasionally result in non-terminated GATK processes - Clear these manually to prevent memory loading
## Arguments:
#### Required arguments

|ARGUMENT|TYPE|DESCRIPTION|
|---|---|---|
|**-i** --input         |Folder             |A directory path containing the .bam files to be used in the CNV analysis.Must end in .bam & have index files present in the same folder|
|**-c** --cohort			|String				|A string or identifying label to mark the appropraite output folder. e.g "RCC" for the	kidney cancer dataset|		
|**-o** --output			|Folder				|A directory path for a preferrably empty folder in which all resulting files will be deposited|
|**-p** --params		|File				|Path for the params.txt file specifically required for the xhmm analysis. Dile is included in the xhmm download|
|**-v** --interval			|File				|Interval file provided by illumina, typically a file named Nextera_Exome_hg38, in tsv bed format, BED4 without header.|
#### Mode tags
* These are provided at any point within the command to specify running modes - e.g. temporary file retention
* They do not require an argument to be provided after their declaration

|ARGUMENT|TYPE|DESCRIPTION|
|---|---|---|
|**-tr** --temp-remove		|Argument			|Providing the option -tr or --temp-remove, given without a proceeding file/folder/string results in the DELETION of all temporary files generated	during the CNV analysis, inlcuding Depth of coverage.|
|**-call** --call-only		|Argument			|Providing this option allows the pipeline to run without performing the depth of coverage analysis. This is particularly useful for re-calling data with new parameters or cutoffs|

#### Optional Arguments
* These options are not required to run as defaults have been set internally
* Alterations to GATK and REF are encouraged based on the data set, genome build and server
* All other values pertain to the data normalisation and CNV calling stringency - See [XHMM documentation](http://atgu.mgh.harvard.edu/xhmm/tutorial.shtml)

|ARGUMENT|TYPE|DESCRIPTION|
|---|---|---|
|**-gatk** --gatk			|File				|Location of the GenomeAnalysisTK.jar file used for depth of coverage analysis|
|**-r** --ref			|File				|Reference fasta for the genome used in the analysis. Reference must match contigs used in bed file provided - GATK will error if reference does not match the bed.|
|**-minTS**			|Integer				|Set minimum allowed target size for target filtering in xhmm --matrix step. Higher values are more stringent|
|**-maxTS**			|Integer				|Set maximum allowed target size for target filtering in xhmm --matrix step. Lower values are more stringent|
|**-minMeanTRD**	|Integer				|Set minimum mean value read depth for each target in the interval file. Targets with mean Read depth below this value are removed from the analysis|
|**-maxMeanTRD**	|Integer				|Set maximum mean value read depth for each target in the interval file - targets with mean Read depth above this value are removed from the analysis|
|**-minMeanSRD**    |Integer                         |Set minimum mean value read depth for each sample in the interval file - samples with mean Read depth below this value are removed from the analysis|
|**-maxMeanSRD**    |Integer                         |Set maximum mean value read depth for each sample in the interval file - samples with mean Read depth above this value are removed from the analysis|
|**-maxSampleSD**	|Integer				|Set the maximum allowed stardard deviation for each sample - Samples are removed if their SD exceeds this value|
|**-maxTargetSD** 	|Integer				|Set the maximum allowed stardard deviation for each target - Targets are removed if their SD exceeds this value|
|**-pca**			|Float				|PCA cut-off used for removing variance - Principle components exceeding this value as a percentage are normalised and removed from the data to reduce noise within the data|
---
## Dependencies:
Currently only tested and implemented on **Whisperwind**  _without modification_


| Softwares:                                               | R Packages:                                                                 | Core utilities:            |
|----------------------------------------------------------|-----------------------------------------------------------------------------|----------------------------|
|[xHMM](http://atgu.mgh.harvard.edu/xhmm/citations.shtml) | [ggplot2](http://ggplot2.tidyverse.org/reference/)                          | vim |
| [GATK](https://software.broadinstitute.org/gatk/)        | [dplyr](http://dplyr.tidyverse.org/reference/)                              | cat                        |
| [R](https://www.rproject.org/)                           | [tidyr](http://tidyr.tidyverse.org/)                                        | find                       |
|                                                          | [data.table](https://github.com/Rdatatable/data.table)                      | sort                       |
|                                                          | [stringr](https://github.com/tidyverse/stringr)                             | sed                        |
|                                                          | [ggrepel](https://github.com/slowkow/ggrepel)                               | awk                        |
|                                                          | [reshape2](https://www.rdocumentation.org/packages/reshape2/versions/1.4.2) |                            |

---
## Examples:
##### General syntax - Run from the git directory
`./CNV_analysis.sh [Required arguments] [optional arguments] [Mode tags]`
##### Specific cases:
###### Running xhmm analysis on a folder containing bams & remove temp files:
```
./CNV_analysis.sh -c My_project -i /data/BAMS/ -o /data/CNV_REULTS/ -p params.txt -v ref_exome.bed -tr
```
###### Running xhmm call-only mode with new values for minimum Target size:
```
./CNV_analysis.sh -c Existing_project -i /data/BAMS/ -o /data/CNV_REULTS/ -p params.txt -v ref_exome.bed -minTS 10 -maxTS 1000 -call
```
