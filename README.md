
##Pipeline to effectively assign OTU from barcoded amplicon sequence reads adapted from Mathon et al 2021  https://doi.org/10.1111/1755-0998.13430 
## 
## 
##Dr. Christoph Deeg 2025##In brief, the pipeline does the following in the 
##1.	Read merging
##2.	Demultiplexing based on integrated barcode sequence (in primer)
##3.	Primer removal
##4.	Annotate for tracking
##5.	Dereplicate to amplicon sequence variants (ASV)
##6.	Quality filter
##7.	Cluster reads
##8.	Annotate clusters
##9.	BLAST for tax ID
##10.	Run MEGAN on BLAST output to identify and inspect OTUs.
##11.	Run R script on output to assign read numbers to ASV with OTU assignment and combine for OTU
################################################################################
##requirements:
#OBItools
#vserach
#cutadapt
#swarm
#python2
#python3
#BLAST
#MEGAN
#R
################################################################################

## Usage:
## Copy the directory architecture to destination and prepare the following
##1. Unzipped raw read data in "02_raw_data" and rename
##2. Add tag file for demultiplexing in "00_archive"
##3. Edit a version of this script with the appropriate primer sequences
##4. Edit this file appropriately
##5. Go to pipeline main directory and run the following
##    ./01_scripts/optimized_pipeline_PSSI_example.sh
