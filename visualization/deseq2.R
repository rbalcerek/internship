library(DESeq2)
library(tidyverse)
library(pheatmap)

counts <- read.csv("/home/rachel/kreport/raw_counts.csv", row.names = 1)
meta <- read.csv("/home/rachel/kreport/meta.csv")

counts <- counts[, meta$sample]
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ condition)

dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "after", "before"))
res_df <- as.data.frame(res)
res_df$taxon <- rownames(res_df)


res_genus <- res_df %>% filter(str_detect(taxon, "^Genus:"))
res_species <- res_df %>% filter(str_detect(taxon, "^Species:"))
write.csv(res_genus, "/home/rachel/kreport/deseq2_results_genus.csv", row.names = FALSE)
write.csv(res_species, "/home/rachel/kreport/deseq2_results_species.csv", row.names = FALSE)

res_df <- res_df %>%
  mutate(sig = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "Not significant"))

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

sig_taxa_names <- rownames(res)[which(res$padj < 0.05)]
vsd_mat <- assay(vsd)[sig_taxa_names, ]
up_taxa <- res_df %>%
  filter(padj < 0.05, log2FoldChange > 0) %>%
  arrange(padj)

down_taxa <- res_df %>%
  filter(padj < 0.05, log2FoldChange < 0) %>%
  arrange(padj)


res_df <- res_df %>%
  mutate(classification = case_when(
    padj < 0.05 ~ "Significant",
    padj >= 0.05 & padj < 0.1 ~ "Near-significant",
    padj >= 0.1 ~ "Not significant",
    TRUE ~ "NA"
  ))

write.csv(res_df, "/home/rachel/kreport/deseq2_results_classified.csv", row.names = FALSE)

res_sig_or_near <- res_df %>%
  filter(classification != "Not significant") %>%
  arrange(padj) %>%
  slice_min(padj, n = 40) %>%
  mutate(taxon = factor(taxon, levels = taxon)) 

ggplot(res_sig_or_near, aes(x = taxon, y = log2FoldChange, fill = classification)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Significant" = "red", "Near-significant" = "orange")) +
  theme_minimal() +
  labs(title = "Top 40 Differentially Abundant Taxa (with Significance Classes)",
       x = "Genus or Species",
       y = "log2 Fold Change (After/Before)",
       fill = "Significance") +
  theme(axis.text.y = element_text(size = 6))

