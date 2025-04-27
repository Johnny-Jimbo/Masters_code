#The aim of this code is to load in braker, eggnogg and stringtie merge files and merge them all together so we have on big informative document.
#Here i am loading in the required packages and setting my working directive
  setwd("C:/Users/danie/Desktop/Eggnogg_Braker")
  library(readr)
  library(rtracklayer)
  library(dplyr)

#Here i am loading in my ocellatus files
  OceEgg = read_tsv("C:/Users/danie/Desktop/Eggnogg_Braker/SynSpl_eggnog.tsv")
  OceBraker = readGFF("C:/Users/danie/Desktop/Eggnogg_Braker/SynOceBraker.gtf")
  OceMerge = readGFF("C:/Users/danie/Desktop/Eggnogg_Braker/OcellatusMerge.gtf")

#Here i am load in my Picturatus files
  PicEgg = read_tsv("C:/Users/danie/Desktop/Eggnogg_Braker/SynPic1_eggnog.tsv")
  PicBraker = readGFF("C:/Users/danie/Desktop/Eggnogg_Braker/SynPicBraker.gtf")
  PicMerge = readGFF("C:/Users/danie/Desktop/Eggnogg_Braker/PictuMerge.gtf")

#The reason we have no splendidus, is because there was / is no braker file for it

#Here i just wanted to note down what variables were in the different files
#eggnogg file: Query, seed ortholog, evalue, score, eggnog OGs, max annotated level, COG category, description, prefered name, GOs, EC, KEGG ko, KEGG pathway, KEGG module, KEGG reaction, KEGG rclass, BRITE, KEGG TC, CAZy, BiGG reaction, and PFAMs.
#Merge GTF file: Sequence id, source, type, start, end, score, strand, phase, gene id, transcritp id, and exon number.
#Braker GTF file: sequence id, Source, type, start, end, score, strand, phase, gene id, transcritp id, and exon number.

#Here i am using the full join function from dplyr to filter our 2 data files by seqid, start, end and strand. so if both datafiles have all 4 things match, they are noted.
  Oce_exact_matches <- full_join(
    OceBraker,
    OceMerge,
    by = c("seqid", "start", "end", "strand")
  )

  Pic_exact_matches <- full_join(
    PicBraker,
    PicMerge,
    by = c("seqid", "start", "end", "strand")
  )

#Here i filter out anything that isn't coding sequences (could also be done for exons)
#Initial idea was to filter by transcript, but there is no start and end matches, likely meaning this entire methode doesn't work
  Picturatus <- Pic_exact_matches[Pic_exact_matches$type.x =="CDS", ]
  Ocellatus <- Oce_exact_matches[Oce_exact_matches$type.x == "CDS", ]

