library(tidyverse)
library(vegan)
library(pheatmap)

abundance <- read.csv("/home/rachel/kreport/abundancegen5k.csv", row.names = 1)
meta <- read.csv("/home/rachel/kreport/meta.csv")

keywords <- c(
  "Methanosarcina", "Methanosphaera", "Methanobrevibacter",
  "Ruminococcus", "Fibrobacter", "Selenomonas",
  "Butyrivibrio", "Streptococcus bovis", "Prevotella",
  "Clostridium", "Wolinella", "Denitrobacterium", "Treponema"
)

abundance_filtered <- abundance %>%
  rownames_to_column("taxon") %>%
  filter(str_detect(taxon, paste(keywords, collapse = "|"))) %>%
  column_to_rownames("taxon")


otu_matrix <- t(as.matrix(abundance_filtered))
otu_matrix <- otu_matrix[, apply(otu_matrix, 2, var) > 0]


shannon <- diversity(otu_matrix, index = "shannon")
inverse_simpson <- 1 / diversity(otu_matrix, index = "simpson")
pielou <- shannon / log(rowSums(otu_matrix > 0))

diversity_df <- data.frame(
  sample = rownames(otu_matrix),
  Shannon = shannon,
  Inverse_Simpson = inverse_simpson,
  Pielou = pielou
) %>%
  left_join(meta, by = "sample") %>%
  mutate(label = paste(vache_id, condition, sep = "_"))


bray_dist <- vegdist(otu_matrix, method = "bray")
bray_matrix <- as.matrix(bray_dist)

sample_labels <- meta %>%
  mutate(label = paste(vache_id, condition, sep = "_")) %>%
  column_to_rownames("sample")

if (all(rownames(bray_matrix) %in% rownames(sample_labels))) {
  new_labels <- sample_labels[rownames(bray_matrix), "label"]
  rownames(bray_matrix) <- new_labels
  colnames(bray_matrix) <- new_labels
}

pheatmap(bray_matrix,
         clustering_distance_rows = "euclidean",  # méthode de clustering (optionnelle)
         clustering_distance_cols = "euclidean",
         clustering_method = "average",  # méthode de linkage (optionnelle),
         main = "Beta Diversity (Bray-Curtis) of selected taxa",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA)

