#gene viewer
#https://nvelden.github.io/geneviewer/articles/geneviewer.html

#The aim of this code is to load in the data about our opsins and produce some plots best representing their respective gene / cluster.
#Installing the package
  install.packages("geneviewer")

#Step 2
#Loading the package and setting the working directory
  setwd("C:/Users/danie/Desktop/Opsins genes")
  library(geneviewer)
  library(stringr)
  library(GenomicRanges)

#step 3
#here i load in my relevant gff file and seperate our target values as we want to look at the specific ones. (i filter out mRNA as we are only interested in CDs)
  Oce_gene_data <- read.table("C:/Users/danie/Desktop/Opsins genes/SynOce1.gff",
                          sep = "\t", header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
  Splen_gene_data <- read.table("C:/Users/danie/Desktop/Opsins genes/GCA_019776025_1.gff",
                              sep = "\t", header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
  Pic_gene_data <- read.table("C:/Users/danie/Desktop/Opsins genes/SynPic1.gff",
                                sep = "\t", header = FALSE, comment.char = "#", stringsAsFactors = FALSE)

#Here i am defining the column names so its more legible
  colnames(Oce_gene_data) <- c("Seqid", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes")
  colnames(Splen_gene_data) <- c("Seqid", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes")
  colnames(Pic_gene_data) <- c("Seqid", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes")

#Here i go into our "attributes" variable and extract the target and parent, so i know what the query seq was and which chromosome / scaffold they are located in
  Oce_gene_data$Target <- str_extract(Oce_gene_data$Attributes, "Target=[^;]+")
  Oce_gene_data$Parent <- str_extract(Oce_gene_data$Attributes, "Parent=[^;]+")
  Splen_gene_data$Target <- str_extract(Splen_gene_data$Attributes, "Target=[^;]+")
  Splen_gene_data$Parent <- str_extract(Splen_gene_data$Attributes, "Parent=[^;]+")
  Pic_gene_data$Target <- str_extract(Pic_gene_data$Attributes, "Target=[^;]+")
  Pic_gene_data$Parent <- str_extract(Pic_gene_data$Attributes, "Parent=[^;]+")

#here im simply removing the "target" and "parent" part of the varaible text so i only have the information that is actually needed
  Oce_gene_data$Target <- sub("Target=", "", Oce_gene_data$Target)
  Oce_gene_data$Parent <- sub("Parent=", "", Oce_gene_data$Parent)
  Splen_gene_data$Target <- sub("Target=", "", Splen_gene_data$Target)
  Splen_gene_data$Parent <- sub("Parent=", "", Splen_gene_data$Parent)
  Pic_gene_data$Target <- sub("Target=", "", Pic_gene_data$Target)
  Pic_gene_data$Parent <- sub("Parent=", "", Pic_gene_data$Parent)

#Step 4
#Here i filter out everything that isn't our target of interest, specifically here i am filtering for mRNA as we want the whole gene
  Oce_gene_data_mrna <- subset(Oce_gene_data, Type == "mRNA")
  Splen_gene_data_mrna <- subset(Splen_gene_data, Type == "mRNA")
  Pic_gene_data_mrna <- subset(Pic_gene_data, Type == "mRNA")

#Here i filter for specific opsins (sws1, sws2, rh1 or rh2), after filtering in cases of single genes the obs. with highest score is selected and filtered for. In the case of gene clusters alle individual genes are selected based on score and then filtered for. 
  Oce_gene_cluster <- subset(Oce_gene_data_mrna, grepl("rh1_exorh-Oreochromis_niloticus---rna-XM_003438995.5_1|rh1_exorh-Hippocampus_zosterae---rna-XM_052057555.1_1", Target))
  Splen_gene_cluster <- subset(Splen_gene_data_mrna, grepl("rh1_exorh-Oreochromis_niloticus---rna-XM_003438995.5_1|rh1_exorh-Hippocampus_zosterae---rna-XM_052057555.1_1 ", Target))
  Pic_gene_cluster <- subset(Pic_gene_data_mrna, grepl("rh1_exorh-Oreochromis_niloticus---rna-XM_003438995.5_1|rh1_exorh-Hippocampus_zosterae---rna-XM_052057555.1_1 1 354|rh1_exorh-Hippocampus_zosterae---rna-XM_052057555.1_1", Target))

#Here i am making another variable called label, which is defined at the start and end combined.
  Oce_gene_cluster$Label <- paste0(Oce_gene_cluster$Start, "-", Oce_gene_cluster$End)
  Splen_gene_cluster$Label <- paste0(Splen_gene_cluster$Start, "-", Splen_gene_cluster$End)
  Pic_gene_cluster$Label <- paste0(Pic_gene_cluster$Start, "-", Pic_gene_cluster$End)


#Step 5
#Here i make the graph itself
  GC_chart(Oce_gene_cluster, 
           start = "Start",
           end = "End",
           group = "Label", 
           height = "150px",
           cluster = "Seqid") %>%
           GC_genes(marker="rbox", marker_size = "small") %>%
           GC_title(title = "Synchiropus Ocellatus rh1") %>%
           GC_color(colorScheme = "schemePastel2")

  GC_chart(Splen_gene_cluster, 
           start = "Start",
           end = "End",
           group = "Label", 
           height = "150px",
           cluster = "Seqid") %>%
           GC_genes(marker="rbox", marker_size = "small") %>%
           GC_title(title = "Synchiropus Splendidus rh1") %>%
           GC_color(colorScheme = "schemePastel2")

  GC_chart(Pic_gene_cluster, 
           start = "Start",
           end = "End",
           group = "Label", 
           height = "150px",
           cluster = "Seqid") %>%
           GC_genes(marker="rbox", marker_size = "small") %>%
           GC_title(title = "Synchiropus Picturatus rh1") %>%
           GC_color(colorScheme = "schemePastel2")

