library(tidyverse)
library(readr)
library(ggpubr)

kreport_path <- "/home/rachel/kreport"
meta <- read_csv(file.path(kreport_path, "meta.csv"))
kreport_files <- list.files(path = kreport_path, pattern = "_kreport$", full.names = TRUE)

parse_kreport_bacteria <- function(file) {
  df <- read_tsv(file, col_names = FALSE, col_types = cols(.default = "c"))
  colnames(df) <- c("percent", "reads_clade", "reads_direct", "rank_code", "taxid", "name")
  
  df <- df %>%
    mutate(
      across(c(percent, reads_clade, reads_direct), as.numeric),
      sample = str_remove(basename(file), "_kreport$")
    )

  bacteria_row <- which(df$name == "Bacteria" & df$rank_code == "D")
  next_D_row <- which(df$rank_code == "D" & seq_along(df$name) > bacteria_row)[1]
  
  if (length(bacteria_row) == 0 || is.na(bacteria_row)) return(tibble())
  
  end_row <- if (!is.na(next_D_row)) next_D_row - 1 else nrow(df)
  bacteria_df <- df[(bacteria_row + 1):end_row, ]
  
  if (nrow(bacteria_df) == 0) return(tibble())
  
  bacteria_df <- bacteria_df %>%
    filter(rank_code == "P", reads_clade > 5000)
  
  if (nrow(bacteria_df) == 0) return(tibble())
  
  total_reads <- sum(bacteria_df$reads_clade)
  
  bacteria_df <- bacteria_df %>%
    mutate(rel_abundance = (reads_clade / total_reads) * 100)
  
  return(bacteria_df)
}


bacteria_data <- map_dfr(kreport_files, parse_kreport_bacteria)

bacteria_data <- bacteria_data %>%
  left_join(meta, by = "sample") %>%
  mutate(clean_taxon = str_trim(name))

plot_data_bacteria <- bacteria_data %>%
  group_by(condition, clean_taxon) %>%
  summarise(rel_abundance = sum(rel_abundance), .groups = "drop")

total_abund_bact <- plot_data_bacteria %>%
  group_by(condition) %>%
  summarise(total = sum(rel_abundance))

plot_data_bacteria <- plot_data_bacteria %>%
  left_join(total_abund_bact, by = "condition") %>%
  mutate(global_rel_abundance = (rel_abundance / total) * 100)

plot_data_bacteria <- plot_data_bacteria %>%
  mutate(grouped_phylum = ifelse(global_rel_abundance < 1, "Other genera", clean_taxon)) %>%
  group_by(condition, grouped_phylum) %>%
  summarise(global_rel_abundance = sum(global_rel_abundance), .groups = "drop")

plot_data_bacteria$condition <- factor(
  plot_data_bacteria$condition,
  levels = c("before", "after")
)

ggplot(plot_data_bacteria, aes(x = condition, y = global_rel_abundance, fill = grouped_phylum)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Relative abundance of bacterial phyla",
    subtitle = "Summed across all cows (Filtered: >5000 reads, Phylum level)",
    x = "Condition",
    y = "Relative abundance (%)",
    fill = "Phylum"
  ) 




