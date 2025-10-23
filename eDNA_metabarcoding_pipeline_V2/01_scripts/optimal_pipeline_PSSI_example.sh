	#!/bin/bash
###############################################################################
##Pipeline to effectively assign OTU from barcoded amplicon sequence reads adapted from Mathon et al 2021  https://doi.org/10.1111/1755-0998.13430 
## 
##Dr. Christoph Deeg 2025
##
##In brief, the pipeline does the following in the 
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

##
##
##
###############################################################################

## load config global variables


# Global variables:
## adjust these to your system capacity and the amplicon you are working with
INTERP_FOLDER="00_archive"
SCRIPT_FOLDER="01_scripts"
RAW_FOLDER="02_raw_data"
MERGED_FOLDER="03_merged"
SAMPLE_FOLDER="04_samples"
LMIN=55 # lower size limit
LMAX=310 # upper size limit
CORES=75 # cores to use for blastn, clustering, readmerging, etc

# The variable for your output names for the results and the artificial positive controls (APC; optional)
Blast_NAME=$libname"_YOUR_OUTPUT_short_x"
Blast_NAME_APC=$libname"_YOUR_OUTPUT_APC_short_x"

#YOUR primers used for this amplicon (edit these!)
ADAPT1="CGAGAAGACCCTRTGGAGCT"
ADAPT2="GGATTGCGCTGTTATCCCT"


#
# Promt for file name, this should be the name of your input file for each read, e.g. "YOUR_FILE_R1.fastq" and "YOUR_FILE_R2.fastq"
echo What chip?
read libname




mkdir 10_log_files


	
## Step 1: Merge the read files
ls -1 $RAW_FOLDER/*.fastq | \
 # Remove path
	  perl -pe 's/R[12]\.fastq//' | \
    sort -u | \
    while read i 
    do
        echo $i
        vsearch --fastq_mergepairs $i"R1.fastq" --fastq_maxdiffpct 25 --reverse $i"R2.fastq" --fastq_allowmergestagger --threads $CORES --fastaout $i"ali.fq"
    done
	

mv $RAW_FOLDER/*ali.fq $MERGED_FOLDER

##Step 2: Demultiplex by assign each sequence to a sample based on barcode on primer
ls -1 $MERGED_FOLDER/*ali.fq | \
    # Remove path
    awk -F"/" '{ print $2 }' | \
    
    # Remove end of name   
    perl -pe 's/ali\.fq//' | \
    sort -u | \
    # Run cutadapt on each input file
    while read i
    do
        echo $i
		cutadapt $MERGED_FOLDER/$i"ali.fq" -g file:$INTERP_FOLDER/$libname"_tags.fasta" -y '; sample={name};' -e 0 -O 8 -j $CORES --revcomp -o ./$MERGED_FOLDER/$i"ali_assi.fasta" \
		--untrimmed-output $MERGED_FOLDER/$i"unidentified.fasta"
    done
	
	# Merging multiple files consecutively
ls -1 $MERGED_FOLDER/*assi.fasta | \
    # Remove path
    awk -F"/" '{ print $2 }' | \
    
    # Remove end of name   
      perl -pe 's/assi\.fasta//' | \
    sort -u | \
    
    # Step3 : Remove the sequencing primers sequences
    # Note that primers are being removed from the 5' end of forward and reverse reads
    while read i
    do
        echo $i
        cutadapt -g $ADAPT1 -a $ADAPT2 -o $SAMPLE_FOLDER/$i"assi_noprime.fasta" -j $CORES --revcomp $MERGED_FOLDER/$i"assi.fasta" 
    done
	
	
	
## Step 4 annotate each sequence with a unique ID
	
ls -1 $SAMPLE_FOLDER/*.fasta | \
    perl -pe 's/\.fasta//' | \
    
    # Remove the directory part of name 
    awk -F/ '{ print $2 }' - | \

    sort -u | \
    while read i 
    do
        echo $i

        # Run obiannotate 
        obiannotate -S sample:"$i" $SAMPLE_FOLDER/$i".fasta" > $SAMPLE_FOLDER/$i"_annot.fa"
    done
## Step 5: Dereplicate identical sequnces with vsearch
	
ls -1 $SAMPLE_FOLDER/*.fa | \
    perl -pe 's/\.fa//' | \
    
    # Remove the directory part of name 
    awk -F/ '{ print $2 }' - | \

    sort -u | \
    while read i 
    do
        echo $i

        # Run vsearch 
        vsearch --derep_fulllength $SAMPLE_FOLDER/$i".fa" --sizeout --fasta_width 0 --notrunclabels --threads $CORES --relabel_keep --minseqlength 1 --output $SAMPLE_FOLDER/$i"_derep.fa" 
    done
	
	
	
	# convert to obi forat
	
ls -1 $SAMPLE_FOLDER/*derep.fa | \
    perl -pe 's/\.fa//' | \
    
    # Remove the directory part of name 
    awk -F/ '{ print $2 }' - | \

    sort -u | \
    while read i 
    do
        echo $i

        # Format the output for OBI tools 
        python2 $SCRIPT_FOLDER/vsearch_to_obifasta.py -f $SAMPLE_FOLDER/$i".fa" -o $SAMPLE_FOLDER/$i"_format.fa"
    done
	
	# Keep sequences longuer than 55bp without ambiguous bases
	
	ls -1 $SAMPLE_FOLDER/*_format.fa | \
    perl -pe 's/\.fa//' | \
    
    # Remove the directory part of name 
    awk -F/ '{ print $2 }' - | \
 ## Step 6: Quality and size filter
    sort -u | \
    while read i 
    do
        echo $i

        # Run vsearch 
        vsearch --fastx_filter $SAMPLE_FOLDER/$i".fa" --notrunclabels --threads $CORES --fastq_maxns 0 --fastq_minlen $LMIN --fastaout $SAMPLE_FOLDER/$i"_55.fa"
    done
	
	
# Length annotate and format the header
	
	ls -1 $SAMPLE_FOLDER/*_55.fa | \
    perl -pe 's/\.fa//' | \
    
    # Remove the directory part of name 
    awk -F/ '{ print $2 }' - | \
    # annotate with the OBI IDS
    sort -u | \
    while read i 
    do
        echo $i
        obiannotate -S 'size:count' $SAMPLE_FOLDER/$i".fa" | python3 $SCRIPT_FOLDER/formate_header.py > $SAMPLE_FOLDER/$i"_header.fa"
    done	
	
## Step 7 Cluster read variants around centromer using swarm and record the abundance
	
	ls -1 $SAMPLE_FOLDER/*_header.fa | \
    perl -pe 's/\.fa//' | \
    
        awk -F/ '{ print $2 }' - | \

    sort -u | \
    while read i 
    do
        echo $i

        # Run swarm 
        swarm $SAMPLE_FOLDER/$i".fa" -z -f -t $CORES -w $SAMPLE_FOLDER/$i"_clean.fa" -u $SAMPLE_FOLDER/$i"u_clusters"   
		done	 
	# Format the output
	sed -i 's/;/; /g' $SAMPLE_FOLDER//*_clean.fa
	sed -i 's/:/: /g' $SAMPLE_FOLDER//*_clean.fa
	sed -i 's/SUB;/SUB/g' $SAMPLE_FOLDER//*_clean.fa
	sed -i 's/}/}; /g' $SAMPLE_FOLDER//*_clean.fa
 
 # clean up the output file
 
	ls -1 $SAMPLE_FOLDER/*_clean.fa | \
    perl -pe 's/\.fa//' | \
    
        awk -F/ '{ print $2 }' - | \

    sort -u | \
    while read i 
    do
        echo $i

        # clean that up 
	 obiannotate  --delete-tag=scientific_name_by_db --delete-tag=obiclean_samplecount \
	 --delete-tag=obiclean_count --delete-tag=obiclean_singletoncount \
	 --delete-tag=obiclean_cluster --delete-tag=obiclean_internalcount \
	 --delete-tag=obiclean_head  --delete-tag=obiclean_headcount \
	 --delete-tag=id_status --delete-tag=rank_by_db --delete-tag=obiclean_status \
	 --delete-tag=seq_length_ori --delete-tag=sminL --delete-tag=sminR \
	 --delete-tag=reverse_score --delete-tag=reverse_primer --delete-tag=reverse_match --delete-tag=reverse_tag \
	 --delete-tag=forward_tag --delete-tag=forward_score --delete-tag=forward_primer --delete-tag=forward_match \
	 --delete-tag=tail_quality --delete-tag=mode --delete-tag=seq_a_single $SAMPLE_FOLDER/$i".fa" > $SAMPLE_FOLDER/$i"_cf.fa"
 		done	 
 
# # Step 8: Annotate the clusters with abundance information
# sort by abundance
 
	ls -1 $SAMPLE_FOLDER/*_cf.fa | \
    perl -pe 's/\.fa//' | \
    
        awk -F/ '{ print $2 }' - | \

    sort -u | \
    while read i 
    do
        echo $i

        # clean that up and count the occurances
	 obisort -k count -r $SAMPLE_FOLDER/$i".fa" > $SAMPLE_FOLDER/$i"_sort.fa"
 		done	 
		
  # make unique ID for obifasta
 
	ls -1 $SAMPLE_FOLDER/*_sort.fa | \
    perl -pe 's/\.fa//' | \
    
        awk -F/ '{ print $2 }' - | \

    sort -u | \
    while read i 
    do
        echo $i

        # clean that up 
	 python3 $SCRIPT_FOLDER/unique_id_obifasta.py $SAMPLE_FOLDER/$i".fa" > $SAMPLE_FOLDER/$i"_id.fa"
 		done	 
  
  
  # make the OBI tab format file
	ls -1 $SAMPLE_FOLDER/*sort_id.fa | \
    perl -pe 's/\.fa//' | \
    sort -u | \
    while read p
    do
        echo $p
        #export
        obitab --output-seq $p".fa" > $p"_tab.txt"
    done
  
## Step 9: Run a local blast for OTU assignment 

blastn -db /data1/genomics/Programs/nt_database/nt -query $SAMPLE_FOLDER/*id.fa -num_threads $CORES -max_target_seqs 20 -outfmt 5 -out $SAMPLE_FOLDER/$Blast_NAME 
 # repeat this for the APC database (if applicable)
blastn -db /data1/cdeeg/DB/IYS_APC/IYS_APC -query $SAMPLE_FOLDER/*id.fa -num_threads $CORES  -max_target_seqs 5 -outfmt 5 -out $SAMPLE_FOLDER/$Blast_NAME_APC 

## Step 10: Run Megan on BLAST output
#A:Import from blast with the following LCA parameters use the following LCA parameters:
#min score 100
#max expected 0.00000001
#min % ID 97
#top % 5
#min support % 0 (off)
#min support 1
#weighted
#B:Inspect all OTUs of interest, disable suspicious assignments, if satisfied with results, select all, export to csv, select "readName_to_taxonPath", "assigned", and "tab".

##Step 11: Your key output files are the cluster file ("YOUR_FILE_ali_assi_noprime_annot_derep_format_55_headeru_clusters"), your sequence ID tab file 
#("YOUR_FILE_ali_assi_noprime_annot_derep_format_55_header_clean_cf_sort_id_tab.txt"), your MEGAN output ("YOUR_FILE_PSSI_YOUR_OUTPUT_short_x-ex.txt"), 
#your APC MEGAN output if applicable ("YOUR_FILE_PSSI_YOUR_OUTPUT_APC_short_x-ex.txt"). 
#
# Providin the input data above, run the R script "PSSI_count_to_anno_example.R" to get the ASV table and the OTU table for each amplicon. 
