library(tidyverse)
library(vegan)
library(pheatmap)

abundance <- read.csv("/home/rachel/kreport/abundancegen5k.csv", row.names = 1)
meta <- read.csv("/home/rachel/kreport/meta.csv")

keywords <- c(
  "Methanosarcina", "Methanosphaera", "Methanobrevibacter",
  "Ruminococcus", "Fibrobacter", "Selenomonas",
  "Butyrivibrio", "Streptococcus bovis", "Prevotella",
  "Clostridium", "Wolinella"
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

abundance <- abundance %>%
  rownames_to_column(var = "taxon") %>%
  pivot_longer(-taxon, names_to = "sample", values_to = "abundance")

df <- abundance %>%
  left_join(meta, by = "sample")


df_interest <- df %>%
  filter(str_detect(taxon, paste(keywords, collapse = "|")))

df_interest <- df_interest %>%
  mutate(sample_label = paste(vache_id, condition, sep = "_"))

heatmap_mat <- df_interest %>%
  group_by(taxon, sample_label) %>%
  summarise(abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = sample_label, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("taxon") %>%
  as.matrix()

ordered_cols <- meta %>%
  mutate(sample_label = paste(vache_id, condition, sep = "_")) %>%
  distinct(sample_label, vache_id, condition) %>%
  arrange(vache_id, factor(condition, levels = c("before", "after"))) %>%
  pull(sample_label)

ordered_cols <- ordered_cols[ordered_cols %in% colnames(heatmap_mat)]
heatmap_mat <- heatmap_mat[, ordered_cols]

taxon_order <- tibble(taxon = rownames(heatmap_mat)) %>%
  mutate(type = case_when(
    str_detect(taxon, "^Genus:") ~ "Genus",
    str_detect(taxon, "^Species:") ~ "Species",
    TRUE ~ "Other"
  )) %>%
  arrange(factor(type, levels = c("Genus", "Species", "Other"))) %>%
  pull(taxon)

heatmap_mat <- heatmap_mat[taxon_order, ]
rownames(heatmap_mat) <- str_replace(rownames(heatmap_mat), "^Genus:", "")


pheatmap(heatmap_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "row",
         fontsize_row = 8,
         fontsize_col = 10,
         main = "Relative abundances of selected genera")


pcoa_result <- cmdscale(bray_dist, k = 2, eig = TRUE)
pcoa_df <- as.data.frame(pcoa_result$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$sample <- rownames(pcoa_df)

pcoa_df <- pcoa_df %>%
  left_join(meta, by = "sample") %>%
  mutate(label = paste(vache_id, condition, sep = "_"))


var_expl <- round(100 * pcoa_result$eig[1:2] / sum(pcoa_result$eig), 1)

ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = condition, label = label)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 3) +
  theme_minimal() +
  labs(title = "PCoA (Bray-Curtis) – Selected Taxa",
       x = paste0("PCoA1 (", var_expl[1], "%)"),
       y = paste0("PCoA2 (", var_expl[2], "%)"))
