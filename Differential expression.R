#The aim of this code is to load in our data from stringtie and using ballgown and limma i'll performe a differenital expression analysis, for genes upregulated in skin.
setwd("C:/Users/danie/Desktop/Stringtie_results")
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(limma)
library(ggplot2)

# Here i'm defining my dataset
data_directory_Splen = c("C:/Users/danie/Desktop/Stringtie_results/Splendidus/ET_557", 
                         "C:/Users/danie/Desktop/Stringtie_results/Splendidus/ET_570", 
                         "C:/Users/danie/Desktop/Stringtie_results/Splendidus/ET_603", 
                         "C:/Users/danie/Desktop/Stringtie_results/Splendidus/ET_648",
                         "C:/Users/danie/Desktop/Stringtie_results/Splendidus/ET_NonSkinSplen")
data_directory_Pic = c("C:/Users/danie/Desktop/Stringtie_results/Picturatus/ET_612", 
                       "C:/Users/danie/Desktop/Stringtie_results/Picturatus/ET_622",
                       "C:/Users/danie/Desktop/Stringtie_results/Picturatus/ET_YG020-05N1",
                       "C:/Users/danie/Desktop/Stringtie_results/Picturatus/ET_YG020-06N1",
                       "C:/Users/danie/Desktop/Stringtie_results/Picturatus/ET_YG020-07N1")
data_directory_Oce = c("C:/Users/danie/Desktop/Stringtie_results/Ocellatus/ET_630", 
                       "C:/Users/danie/Desktop/Stringtie_results/Ocellatus/ET_641",
                       "C:/Users/danie/Desktop/Stringtie_results/Ocellatus/ET_YG021-05N1",
                       "C:/Users/danie/Desktop/Stringtie_results/Ocellatus/ET_YG021-06N1",
                       "C:/Users/danie/Desktop/Stringtie_results/Ocellatus/ET_YG021-07N1",
                       "C:/Users/danie/Desktop/Stringtie_results/Ocellatus/ET_YG021-08N1")

#For fin to Non-skin comparison
data_directory_Pic = c("C:/Users/danie/Desktop/Stringtie_results/Picturatus/ET_612",
                       "C:/Users/danie/Desktop/Stringtie_results/Picturatus/ET_YG020-05N1",
                       "C:/Users/danie/Desktop/Stringtie_results/Picturatus/ET_YG020-06N1",
                       "C:/Users/danie/Desktop/Stringtie_results/Picturatus/ET_YG020-07N1")
data_directory_Oce = c("C:/Users/danie/Desktop/Stringtie_results/Ocellatus/ET_YG021-05N1",
                       "C:/Users/danie/Desktop/Stringtie_results/Ocellatus/ET_YG021-06N1",
                       "C:/Users/danie/Desktop/Stringtie_results/Ocellatus/ET_YG021-07N1",
                       "C:/Users/danie/Desktop/Stringtie_results/Ocellatus/ET_YG021-08N1")


#For fin to skin comparison
data_directory_Oce = c("C:/Users/danie/Desktop/Stringtie_results/Ocellatus/ET_630", 
                       "C:/Users/danie/Desktop/Stringtie_results/Ocellatus/ET_641",
                       "C:/Users/danie/Desktop/Stringtie_results/Ocellatus/ET_YG021-06N1",
                       "C:/Users/danie/Desktop/Stringtie_results/Ocellatus/ET_YG021-07N1")



# Here i'm loading in the document that defines what is relevant in our samples
Splen_info = read.csv("C:/Users/danie/Desktop/Stringtie_results/Splen_disc.csv", sep = ",")
Pic_info = read.csv("C:/Users/danie/Desktop/Stringtie_results/PictuDesc.csv", sep = ",")
Oce_info = read.csv("C:/Users/danie/Desktop/Stringtie_results/OceDesc.csv", sep = ",")

#For fin to Non-skin comparison
Oce_info = read.csv("C:/Users/danie/Desktop/Stringtie_results/OceDescFin.csv", sep = ",")
Pic_info = read.csv("C:/Users/danie/Desktop/Stringtie_results/PictuDescFin.csv", sep = ",")

#For fin to skin comparison
Oce_info = read.csv("C:/Users/danie/Desktop/Stringtie_results/OceDescSKinFin.csv", sep = ",")


#Here i am making sure that Skin is compared to not skin (it works alphabetiacally, but better to be sure)
Splen_info$Tissue <- factor(Splen_info$Tissue, levels = c("NotSkin", "Skin"))
Pic_info$Tissue <- factor(Pic_info$Tissue, levels = c("NotSkin", "Skin"))
Oce_info$Tissue <- factor(Oce_info$Tissue, levels = c("NotSkin", "Skin"))

#Here i'm creating the ballgown object
bg_splendidus = ballgown(samples = data_directory_Splen, samplePattern = "all", pData = Splen_info)
bg_Picturatus = ballgown(samples = data_directory_Pic, samplePattern = "all", pData = Pic_info)
bg_Ocellatus = ballgown(samples = data_directory_Oce, samplePattern = "all", pData = Oce_info)

#Here i am filtering out low-abundance genes, by removing all transcripts that have a varaince across samples of less than one
bg_splendidus_filt = subset(bg_splendidus,"rowVars(texpr(bg_splendidus)) >1",genomesubset=TRUE)
bg_Picturatus_filt = subset(bg_Picturatus,"rowVars(texpr(bg_Picturatus)) >1",genomesubset=TRUE)
bg_Ocellatus_filt = subset(bg_Ocellatus,"rowVars(texpr(bg_Ocellatus)) >1",genomesubset=TRUE)

#Here i am extracting FPKM data from our ballgown object
fpkm_matrix_Splen <- gexpr(bg_splendidus_filt)
fpkm_matrix_Pic <- gexpr(bg_Picturatus_filt)
fpkm_matrix_Oce <- gexpr(bg_Ocellatus_filt)

#Here i am log2 transforming my FPKM data
log_fpkm_Splen <- log2(fpkm_matrix_Splen + 0.1)
log_fpkm_Pic <- log2(fpkm_matrix_Pic + 0.1)
log_fpkm_Oce <- log2(fpkm_matrix_Oce + 0.1)

#Here i am defining my design matrix used for fitting the linear model
design_Splen <- model.matrix(~ Tissue, data = Splen_info)
design_Pic <- model.matrix(~ Tissue, data = Pic_info)
design_Oce <- model.matrix(~ Tissue, data = Oce_info)

#Here i am fitting the linear model
fit_Splen <- lmFit(log_fpkm_Splen, design_Splen)
fit_Pic <- lmFit(log_fpkm_Pic, design_Pic)
fit_Oce <- lmFit(log_fpkm_Oce, design_Oce)

#Here i am applying the empirical bayes moderation to the linear model
fit_Splen <- eBayes(fit_Splen, trend = TRUE)
fit_Pic <- eBayes(fit_Pic, trend = TRUE)
fit_Oce <- eBayes(fit_Oce, trend = TRUE)

#Here i am getting the results from the linear model, specifically for the skin
results_limma_Splen <- topTable(fit_Splen, coef = "TissueSkin", number = Inf, adjust.method = "BH")
results_limma_Pic <- topTable(fit_Pic, coef = "TissueSkin", number = Inf, adjust.method = "BH")
results_limma_Oce <- topTable(fit_Oce, coef = "TissueSkin", number = Inf, adjust.method = "BH")

#Here i am making the column names a bit more legible
colnames(results_limma_Splen) <- c("Log2FoldChange", "MeanExpression", "T_statistic", "P_value", "Q_value", "LogOdds_DE")
colnames(results_limma_Pic) <- c("Log2FoldChange", "MeanExpression", "T_statistic", "P_value", "Q_value", "LogOdds_DE")
colnames(results_limma_Oce) <- c("Log2FoldChange", "MeanExpression", "T_statistic", "P_value", "Q_value", "LogOdds_DE")

#Here i am filtering the results from the linear model, by both the BH adjusted p value and log2 fold change value.
Splen_sig_genes <- results_limma_Splen[results_limma_Splen$Q_value < 0.001 & abs(results_limma_Splen$Log2FoldChange) > 1, ]
Pic_sig_genes <- results_limma_Pic[results_limma_Pic$Q_value < 0.001 & abs(results_limma_Pic$Log2FoldChange) > 1, ]
Oce_sig_genes <- results_limma_Oce[results_limma_Oce$Q_value < 0.001 & abs(results_limma_Oce$Log2FoldChange) > 1, ]
#Remember the row names are gene name and then the transcript

#Here i am printing out the full results:
write.csv(Splen_sig_genes, "Splendidus_gene_results.csv", row.names=FALSE)
write.csv(Pic_sig_genes, "Picturatus_gene_fin_results.csv", row.names=FALSE)
write.csv(Oce_sig_genes, "Ocellatus_gene_Skin_fin_results.csv", row.names=FALSE)

#Here i am printing only the gene ids
splen_gene_ids <- rownames(Splen_sig_genes)
pic_gene_ids <- rownames(Pic_sig_genes)
oce_gene_ids <- rownames(Oce_sig_genes)

write.csv(splen_gene_ids, "Splen_sig_gene_ids.csv", row.names = FALSE)
write.csv(pic_gene_ids, "Pic_sig_gene_ids.csv", row.names = FALSE)
write.csv(oce_gene_ids, "Oce_sig_gene_ids.csv", row.names = FALSE)

#Here i am plotting the fpkm values for visualization
Color_splen <- ifelse(Splen_info$Tissue == "Skin", "blue", "red")
Color_pic <- ifelse(Pic_info$Tissue == "Skin", "blue", "red")
Color_oce <- ifelse(Oce_info$Tissue == "Skin", "blue", "red")

boxplot(log_fpkm_Splen,
        col=Color_splen,
        las=2,
        ylab="log2(FPKM+0.1)",
        main = "FPKM distribution across Splendidus species")

boxplot(log_fpkm_Pic,
        col=Color_pic,
        las=2,
        ylab="log2(FPKM+0.1)",
        main = "FPKM distribution across Picturatus species")

boxplot(log_fpkm_Oce,
        col=Color_oce,
        las=2,
        ylab='log2(FPKM+0.1)',
        main = "FPKM distribution across Ocellatus species")

#Bargraph for differentially expressed transcripts.
det_counts <- data.frame(
  Species = c("Splendidus", "Picturatus", "Ocellatus"),
  DEG_Count = c(nrow(Splen_sig_genes), nrow(Pic_sig_genes), nrow(Oce_sig_genes))
)

ggplot(det_counts, 
       aes(x = Species, y = DEG_Count, fill = Species)) +
       geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Differentially expressed genes per species",
       x = "Species",
       y = "Number of DEGs") +
  theme(legend.position="none")

#For fin to nonskin comparison
det2_counts <- data.frame(
  Species = c("Picturatus", "Ocellatus"),
  DEG_Count = c(nrow(Pic_sig_genes), nrow(Oce_sig_genes))
)

ggplot(det2_counts, 
       aes(x = Species, y = DEG_Count, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Differentially expressed genes per species",
       x = "Species",
       y = "Number of DEGs") +
  theme(legend.position="none")

