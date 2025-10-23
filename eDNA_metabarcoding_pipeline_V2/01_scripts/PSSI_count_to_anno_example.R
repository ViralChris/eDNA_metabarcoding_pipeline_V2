# Connect read counts to the annotations

# Input: output of obitab ID counts (read counts per cluster per sample), MEGAN annotation (taxonomy ID per cluster), and cluster file (cluster ID)


#### 0. Setup (no changes reqd) ####
# Install Packages
#install.packages("RColorBrewer")
library("RColorBrewer")
library("dplyr")
library("stringr")
library("tidyverse")


options(mc.cores = parallel::detectCores())

# Set working directory depending on the dataset

setwd("C:/01_Data/Chris_Deeg/projects/02_eDNA/WCVI_PSSI/metabarcoding/2024_Q1_Q2_metabarcoding/large_stuff/all_combined")



# Choose your dataset

libraries <- c("12S_1A_1","12S_1B_2","12S_2A_1","12S_2B_2","12S_3A_1","12S_3B_2",
               "16S_1A_1","16S_1B_2","16S_2A_1","16S_2B_2","16S_3A_1","16S_3B_2",
               "12S_1A_2","12S_1B_2","12S_2A_2","12S_2B_2","12S_3A_2","12S_3B_2",
               "16S_1A_2","16S_1B_2","16S_2A_2","16S_2B_2","16S_3A_2","16S_3B_2",
               "COI_1A","COI_1B","COI_2A","COI_2B","COI_3A","COI_3B",
                "Seb_2A","Seb_2B","Seb_3A","Seb_3B"
                )


for(i in libraries) {
 

 datatype <- paste(i)
#datatype <-"16S_1B_1"
  paste("You are analyzing ", datatype, sep = "")



#### 1.0 Import input data and merge #####
#paste("You are analyzing ", datatype, sep = "")

counts <- read.delim2(paste(datatype,"_ali_assi_noprime_annot_derep_format_55_header_clean_cf_sort_id_tab.txt", sep = ""))
annot <- read.delim2(paste(datatype,"_PSSI_cum_2024_Q12_short_x-ex.txt", sep = ""), header = F
                     , col.names = c("id","taxon"))
clusters <- read.delim2(paste(datatype,"_ali_assi_noprime_annot_derep_format_55_headeru_clusters", sep = ""), header = F, col.names = c("type","cluster","length_hits","similarity","orientation","NA1",
                                                                                    "NA2","score","definition","cluster_def"))



head(counts)
head(annot)
head(clusters)

names(counts) 
names(annot)
names(clusters)


### now extract the sample ID and the counts ###

counts$definition <- gsub(" ","", counts$definition)


counts <- counts %>% 
  mutate("sample" = str_extract(definition,  "(?<=\')[^\']+"))


counts <- counts %>% 
  mutate("count" = str_extract(definition,  "(?<=\':)[^}]+"))

counts$count <- sub(" ","", counts$count)



# edit the cluster file

clusters <- subset(clusters,type !="C")
clusters<- 
    clusters %>% select(type, cluster, definition,cluster_def) %>%  mutate(cluster_def = ifelse(type=="S", definition, cluster_def))

# here you will need to edit the name of the seuquencing chip to remove it from the analysos (in this case LH00207R:) and 
# get the read length info for the run (in this case 155). Look at your data and edit accordingly
clusters$definition <- sub("LH00207R:","", clusters$definition)
clusters$cluster_def <- sub("LH00207R:","", clusters$cluster_def)
clusters <-   clusters %>% 
  mutate("clust_def" = str_extract(cluster_def,  "155[^;]+"))
clusters <-   clusters %>% 
  mutate("sample" = str_extract(definition,  "(?<=\')[^\']+"))
clusters <-   clusters %>% 
  mutate("count" = str_extract(definition,  "(?<=\':)[^}]+"))

#get the cluster_def from the counts file

counts <- counts %>% 
  mutate("clust_def" = str_extract(definition,  "155[^;]+"))

# merge the cluster and count file to get ID for all sequence variants in all samples

cluster_count <- merge(clusters,counts,by="clust_def")

cluster_count <- cluster_count[ , -which(names(cluster_count) %in% c( "type", "sample.y", "count.y", "definition.y","type","cluster"))]

colnames(cluster_count) <- sub("*.x","", colnames(cluster_count))

# now merge that with the annotation file to get to OTU for each of the sequence variants

cluster_count$id <- as.character(cluster_count$id)
annot$id <- gsub( ";","", annot$id)

cluster_count_otu <- merge(cluster_count, annot, by="id")

# turn that into an OTU by site table in wide format

cluster_count_otu$count <- as.numeric(cluster_count_otu$count)


#print the ASV (amplicon sequence variant based data
write.csv(cluster_count_otu, paste(datatype, "_count_ASV.csv", sep = ""))

# now narrow this down to OTU level

wide_cluster_count_otu <- cluster_count_otu %>% 
  group_by(taxon, sample) %>% 
  summarize(sum = sum(count)) %>% 
  spread(sample, sum)

#save that as OTU counts per site

otu.site.filename <- paste( datatype, "_OTU_sites.csv", sep = "")
write.csv(x = wide_cluster_count_otu, file = otu.site.filename)
}
