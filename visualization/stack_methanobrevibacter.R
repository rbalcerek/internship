library(tidyverse)
library(readr)

kreport_path <- "/home/rachel/kreport"
meta <- read_csv(file.path(kreport_path, "meta.csv"))
kreport_files <- list.files(path = kreport_path, pattern = "_kreport$", full.names = TRUE)

parse_kreport_archaea_species <- function(file) {
  df <- read_tsv(file, col_names = FALSE, col_types = cols(.default = "c"))
  colnames(df) <- c("percent", "reads_clade", "reads_direct", "rank_code", "ncbi_taxid", "name")
  
  df <- df %>%
    mutate(
      across(c(percent, reads_clade, reads_direct), as.numeric),
      sample = str_remove(basename(file), "_kreport$")
    )
  
  archaea_row <- which(df$name == "Archaea" & df$rank_code == "D")
  next_D_row <- which(df$rank_code == "D" & seq_along(df$name) > archaea_row)[1]
  
  if (length(archaea_row) == 0 || is.na(archaea_row)) return(tibble())
  
  end_row <- if (!is.na(next_D_row)) next_D_row - 1 else nrow(df)
  archaea_df <- df[(archaea_row + 1):end_row, ]
  
  if (nrow(archaea_df) == 0) return(tibble())
  
  archaea_df <- archaea_df %>%
    filter(rank_code == "S", reads_clade > 5000)
  
  if (nrow(archaea_df) == 0) return(tibble())
  
  total_reads <- sum(archaea_df$reads_clade)
  
  archaea_df <- archaea_df %>%
    mutate(rel_abundance = (reads_clade / total_reads) * 100)
  
  return(archaea_df)
}


archaea_species_data <- map_dfr(kreport_files, parse_kreport_archaea_species)

archaea_species_data <- archaea_species_data %>%
  left_join(meta, by = "sample") %>%
  mutate(clean_taxon = str_trim(name))
methanobrevibacter_species <- archaea_species_data %>%
  filter(str_detect(clean_taxon, "Methanobrevibacter"))

methanobrevibacter_species <- methanobrevibacter_species %>%
  mutate(condition = factor(condition, levels = c("before", "after")))

plot_methanobrevibacter_species <- methanobrevibacter_species %>%
  group_by(condition, clean_taxon) %>%
  summarise(rel_abundance = sum(rel_abundance), .groups = "drop")

total_methano <- plot_methanobrevibacter_species %>%
  group_by(condition) %>%
  summarise(total = sum(rel_abundance))

plot_methanobrevibacter_species <- plot_methanobrevibacter_species %>%
  left_join(total_methano, by = "condition") %>%
  mutate(global_rel_abundance = (rel_abundance / total) * 100)

ggplot(plot_methanobrevibacter_species, aes(x = condition, y = global_rel_abundance, fill = clean_taxon)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Relative abundance of Methanobrevibacter species",
    subtitle = "Summed across all cows (Filtered >5000 reads)",
    x = "Condition",
    y = "Relative abundance (%)",
    fill = "Species"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.grid = element_blank()
  )


plot_methanobrevibacter_species_vache <- methanobrevibacter_species %>%
  group_by(vache_id, condition, clean_taxon) %>%
  summarise(rel_abundance = sum(rel_abundance), .groups = "drop")

ggplot(plot_methanobrevibacter_species_vache, aes(x = condition, y = rel_abundance, fill = clean_taxon)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~ vache_id) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Relative abundance of Methanobrevibacter species per cow",
    subtitle = "Species level (Filtered >5000 reads)",
    x = "Condition",
    y = "Relative abundance (%)",
    fill = "Species"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.grid = element_blank()
  )
