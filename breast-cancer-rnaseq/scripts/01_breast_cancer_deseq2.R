# ==============================================================================
# Project: Breast Cancer RNA-Seq Differential Expression Analysis
# Description: Identification of DEGs between Breast Cancer subtypes/treatments
# Author: Bioinfo Portfolio
# ==============================================================================

# 1. Load Required Libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)

# 2. Load / Prepare Counts and Metadata
# (Nous allons charger la matrice de comptage et les métadonnées ici)
# Option A : Chargement d'un jeu de données d'exemple intégré pour le cancer du sein
# (Nous simulons ici une expérience d'expression différentielle ER+ vs ER-)

set.seed(123)

# Génération d'une matrice de comptage synthétique (1000 gènes, 6 échantillons)
# 3 échantillons ER_positive (Contrôle), 3 échantillons ER_negative (Traité)
genes <- paste0("GENE_", 1:1000)
samples <- c("ER_pos_1", "ER_pos_2", "ER_pos_3", "ER_neg_1", "ER_neg_2", "ER_neg_3")

counts <- matrix(
  rpois(6000, lambda = 100), 
  nrow = 1000, 
  ncol = 6,
  dimnames = list(genes, samples)
)

# Ajout de gènes sur-exprimés artificiels pour simuler un vrai signal biologique
counts[1:50, 4:6] <- counts[1:50, 4:6] * 4

# Métadonnées des échantillons
coldata <- data.frame(
  condition = factor(c("ER_positive", "ER_positive", "ER_positive", 
                       "ER_negative", "ER_negative", "ER_negative")),
  row.names = samples
)

# 3. Création de l'objet DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

print(dds)
# 4. Filtrage des gènes peu exprimés
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# 5. Exécution du modèle DESeq2
dds <- DESeq(dds)

# 6. Extraction des résultats de différentiation (ER_negative vs ER_positive)
res <- results(dds)

# 7. Structuration du tableau de résultats
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "gene_id") %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

# Afficher les 10 gènes les plus significatifs
head(res_df, 10)

# 8. Sauvegarde du fichier CSV dans le dossier results
write.csv(res_df, "../results/breast_cancer_degs.csv", row.names = FALSE)

print("Analyse terminée ! Résultats exportés dans results/breast_cancer_degs.csv")
# ==============================================================================
# 9. Data Visualization (Volcano Plot & Heatmap)
# ==============================================================================

# A. Volcano Plot
png("../results/volcano_breast_cancer.png", width = 800, height = 600)
EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    title = 'Volcano Plot - Breast Cancer ER- vs ER+',
    pCutoff = 0.05,
    FCcutoff = 1.0,
    pointSize = 2.0,
    labSize = 4.0)
dev.off()

# B. Heatmap des 20 gènes les plus significatifs
top20_genes <- head(order(res$padj), 20)
mat <- assay(rlog(dds))[top20_genes, ]
mat <- mat - rowMeans(mat)

png("../results/heatmap_breast_cancer.png", width = 800, height = 600)
pheatmap(mat, annotation_col = coldata["condition"], 
         main = "Top 20 DEGs - Breast Cancer ER Status")
dev.off()

print("Graphiques générés et sauvegardés dans le dossier results !")