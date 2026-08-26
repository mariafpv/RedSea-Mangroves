##### Functional comparison of Mangrove and Seagrass MAGs Jaccard #####

#1. Libraries
#2. Read data
#3. Build complete KO matrix
#4. GLOBAL ANALYSIS
# Jaccard
# PCoA
# Plots
# PERMANOVA
# PERMDISP
#5. SELECTED GENE ANALYSIS
# Jaccard
# PCoA
# Plots
# PERMANOVA
# PERMDISP

### 1. Libraries ###
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

theme_set(theme_classic(base_family = "ArialMT"))

## Working directory
setwd("/path/.../") # nolint: line_length_linter.

## Annotation table
ko <- read.csv(
  "tables/combined_annotations.csv",
  sep = ";",
  check.names = FALSE
)
## Metadata
metadata <- read.csv(
  "tables/rsxchina_metadata.csv",
  sep = ";",
  check.names = FALSE
)
## Pathway table
pathway <- read.delim(
  "tables/kos_interest_SN.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

### 2. Build base counts matrix ###
##Extract ko codes values to make the count
ko_table <- ko |>
  # Keep only rows containing at least one valid KO
  filter(str_detect(kegg_ko, "ko:K[0-9]{5}")) |>
  # Split multiple KO annotations
  separate_rows(kegg_ko, sep = ",") |>
  # Remove "ko:" and keep only the KO ID
  mutate(
    kegg_ko = str_extract(kegg_ko, "K[0-9]{5}")
  ) |>
  # Remove rows where no valid KO was extracted
  filter(!is.na(kegg_ko))
## Count the number of KOs assigned to each CDS
ko_table <- ko_table |>
  group_by(query) |>
  mutate(
    n_KO = n(),
    weight = 1 / n_KO
  ) |>
  ungroup()
## Sum weighted counts per MAG
ko_counts <- ko_table |>
  group_by(genome, kegg_ko) |>
  summarise(
    abundance = sum(weight),
    .groups = "drop"
  )
## Correct by completeness
completeness <- metadata |>
  select(genome, completeness) |>
  distinct()

ko_counts <- ko_counts |>
  left_join(
    completeness,
    by = "genome"
  ) |>
  mutate(
    abundance = abundance / (completeness / 100)
  )

## verify
sum(is.na(ko_counts$completeness))

## Check values
ggplot(ko_counts, aes(x = abundance)) +
  geom_histogram(
    bins = 50,
    color = "black",
    fill = "#1f77b4"
  ) +
  labs(
    title = "Distribution of non-zero KO counts",
    x = "Weighted KO count",
    y = "Frequency"
  ) +
  theme_classic(base_family = "ArialMT")
## same plot but log
ggplot(ko_counts, aes(x = abundance)) +
  geom_histogram(
    bins = 50,
    color = "black",
    fill = "#1f77b4"
  ) +
  scale_y_log10() +
  labs(
    x = "Weighted KO count",
    y = "Frequency (log scale)"
  ) +
  theme_classic(base_family = "ArialMT")
## genomes and KOs responsible for the highest corrected abundances
ko_counts |>
  arrange(desc(abundance)) |>
  head(20)
## summary statistics
summary(ko_counts$abundance)
quantile(
  ko_counts$abundance,
  probs = c(0.5, 0.9, 0.95, 0.99, 0.999)
)
##Another plot to check the correction by completeness
ko_counts_raw <- ko_counts
ggplot(ko_counts, aes(x = ko_counts_raw$abundance, y = abundance)) +
  geom_point(alpha = 0.3, size = 0.7) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "red"
  ) +
  labs(
    x = "Raw weighted abundance",
    y = "Completeness-corrected abundance"
  ) +
  theme_classic(base_family = "ArialMT")

### 3. Build complete KO matrix ###
ko_matrix_global <- ko_counts |>
  select(
    genome,
    kegg_ko,
    abundance
  ) |>
  pivot_wider(
    names_from = kegg_ko,
    values_from = abundance,
    values_fill = 0
  )

ko_matrix_global <- as.data.frame(
  ko_matrix_global
)

rownames(ko_matrix_global) <- ko_matrix_global$genome

ko_matrix_global$genome <- NULL

## Metadata
meta_global <- metadata |>
  select(
    genome,
    location,
    compartment,
    domain,
    phylum,
    class,
    order,
    family,
    genus,
    specie
  ) |>
  distinct()

meta_global <- meta_global[
  match(
    rownames(ko_matrix_global),
    meta_global$genome
  ),
]

stopifnot(
  all(
    meta_global$genome ==
      rownames(ko_matrix_global)
  )
)

## Rename locations
meta_global$location <- factor(
  meta_global$location,
  levels = c(
    "red_sea",
    "China"
  ),
  labels = c(
    "Mangroves",
    "Seagrass"
  )
)
## Highlight variable
meta_global <- meta_global |>
  mutate(
    Highlight = case_when(
      order == "Desulfobacterales" &
        location == "Mangroves" ~
        "Mangrove_Desulfobacterales",

      order == "Chromatiales" &
        location == "Mangroves" ~
        "Mangrove_Chromatiales",

      order == "Desulfobacterales" &
        location == "Seagrass" ~
        "Seagrass_Desulfobacterales",

      order == "Chromatiales" &
        location == "Seagrass" ~
        "Seagrass_Chromatiales",

      TRUE ~ "Other"
    )
  )

## SELECTED GENE MATRIX (another matrix for specific genes of interest)
selected_kos <- unique(
  pathway$ko_code
)

ko_matrix_selected <- ko_matrix_global[
  ,
  colnames(ko_matrix_global) %in% selected_kos
]
## Remove genomes without selected genes
ko_matrix_selected <- ko_matrix_selected[
  rowSums(ko_matrix_selected) > 0,
]
## Remove genes absent from all genomes
ko_matrix_selected <- ko_matrix_selected[
  ,
  colSums(ko_matrix_selected) > 0
]
## Metadata for selected matrix
meta_selected <- meta_global[
  match(
    rownames(ko_matrix_selected),
    meta_global$genome
  ),
]

stopifnot(
  all(
    meta_selected$genome ==
      rownames(ko_matrix_selected)
  )
)

### 4. GLOBAL ANALYSIS (JACCARD)
## Presence / Absence matrix
ko_matrix_global_pa <- ko_matrix_global

ko_matrix_global_pa[
  ko_matrix_global_pa > 0
] <- 1

## Jaccard distance
dist_jaccard_global <- vegdist(
  ko_matrix_global_pa,
  method = "jaccard",
  binary = TRUE
)

## PCoA
pcoa_jaccard_global <- cmdscale(
  dist_jaccard_global,
  eig = TRUE,
  k = 2
)

coords_pcoa_jaccard_global <- as.data.frame(
  pcoa_jaccard_global$points
)

colnames(coords_pcoa_jaccard_global) <- c(
  "PC1",
  "PC2"
)

coords_pcoa_jaccard_global <- cbind(
  coords_pcoa_jaccard_global,
  meta_global
)


## Percentage of variance explained
eig_pcoa_jaccard_global <- pcoa_jaccard_global$eig

variance_pcoa_jaccard_global <- round(
  eig_pcoa_jaccard_global[eig_pcoa_jaccard_global > 0] /
    sum(eig_pcoa_jaccard_global[eig_pcoa_jaccard_global > 0]) * 100,
  1
)

## Global PCoA
plot_pcoa_jaccard_global <- ggplot(
  coords_pcoa_jaccard_global,
  aes(PC1, PC2)
) +

  stat_ellipse(
    aes(
      fill = location,
      colour = location
    ),
    geom = "polygon",
    alpha = 0.15,
    linewidth = 0.8
  ) +

  geom_point(
    aes(fill = location),
    shape = 24,
    size = 4.2,
    colour = "black",
    stroke = 0.35
  ) +

  scale_fill_manual(
    values = c(
      Mangroves = "#8C6BB1",
      Seagrass = "#F2E55C"
    )
  ) +

  scale_colour_manual(
    values = c(
      Mangroves = "#8C6BB1",
      Seagrass = "#F2E55C"
    )
  ) +

  xlab(
    paste0(
      "PC1 (",
      variance_pcoa_jaccard_global[1],
      "%)"
    )
  ) +

  ylab(
    paste0(
      "PC2 (",
      variance_pcoa_jaccard_global[2],
      "%)"
    )
  )

plot_pcoa_jaccard_global

## Highlight Desulfobacterales and Chromatiales
plot_pcoa_jaccard_global_h <- ggplot(
  coords_pcoa_jaccard_global,
  aes(PC1, PC2)
) +

  stat_ellipse(
    aes(
      fill = location,
      colour = location
    ),
    geom = "polygon",
    alpha = 0.15,
    linewidth = 0.8
  ) +

  geom_point(
    data = subset(
      coords_pcoa_jaccard_global,
      Highlight == "Other"
    ),
    shape = 24,
    size = 4.2,
    fill = "grey80",
    colour = "black",
    stroke = 0.25
  ) +

  geom_point(
    data = subset(
      coords_pcoa_jaccard_global,
      Highlight != "Other"
    ),
    aes(fill = Highlight),
    shape = 24,
    size = 4.2,
    colour = "black",
    stroke = 0.25
  ) +

  scale_fill_manual(
    values = c(

      Mangroves = "#8C6BB1",
      Seagrass = "#F2E55C",

      Mangrove_Desulfobacterales = "#5E3C99",
      Mangrove_Chromatiales = "#C2A5CF",

      Seagrass_Desulfobacterales = "#E6AB02",
      Seagrass_Chromatiales = "#FFF176"
    )
  ) +

  scale_colour_manual(
    values = c(
      Mangroves = "#8C6BB1",
      Seagrass = "#F2E55C"
    )
  ) +

  xlab(
    paste0(
      "PC1 (",
      variance_pcoa_jaccard_global[1],
      "%)"
    )
  ) +

  ylab(
    paste0(
      "PC2 (",
      variance_pcoa_jaccard_global[2],
      "%)"
    )
  )

plot_pcoa_jaccard_global_h

## PERMANOVA
permanova_pcoa_jaccard_global <- adonis2(
  dist_jaccard_global ~ location,
  data = meta_global,
  permutations = 9999
)

permanova_pcoa_jaccard_global

## PERMDISP
betadisper_pcoa_jaccard_global <- betadisper(
  dist_jaccard_global,
  meta_global$location
)

anova(
  betadisper_pcoa_jaccard_global
)

permutest(
  betadisper_pcoa_jaccard_global,
  permutations = 9999
)

## Save figures
ggsave(
  "results/PCoA_jaccard_global.png",
  plot = plot_pcoa_jaccard_global,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  "results/PCoA_jaccard_global_highlight.png",
  plot = plot_pcoa_jaccard_global_h,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

### 5. SELECTED GENE ANALYSIS (JACCARD)
## Presence / Absence matrix
ko_matrix_selected_pa <- ko_matrix_selected

ko_matrix_selected_pa[
  ko_matrix_selected_pa > 0
] <- 1

## Jaccard distance
dist_jaccard_selected <- vegdist(
  ko_matrix_selected_pa,
  method = "jaccard",
  binary = TRUE
)

## PCoA
pcoa_jaccard_selected <- cmdscale(
  dist_jaccard_selected,
  eig = TRUE,
  k = 2
)

coords_pcoa_jaccard_selected <- as.data.frame(
  pcoa_jaccard_selected$points
)

colnames(coords_pcoa_jaccard_selected) <- c(
  "PC1",
  "PC2"
)

coords_pcoa_jaccard_selected <- cbind(
  coords_pcoa_jaccard_selected,
  meta_selected
)

## Percentage of variance explained
eig_pcoa_jaccard_selected <- pcoa_jaccard_selected$eig

variance_pcoa_jaccard_selected <- round(
  eig_pcoa_jaccard_selected[eig_pcoa_jaccard_selected > 0] /
    sum(eig_pcoa_jaccard_selected[eig_pcoa_jaccard_selected > 0]) * 100,
  1
)

## PCoA plot
plot_pcoa_jaccard_selected <- ggplot(
  coords_pcoa_jaccard_selected,
  aes(PC1, PC2)
) +

  stat_ellipse(
    aes(
      fill = location,
      colour = location
    ),
    geom = "polygon",
    alpha = 0.15,
    linewidth = 0.8
  ) +

  geom_point(
    aes(fill = location),
    shape = 24,
    size = 4.2,
    colour = "black",
    stroke = 0.35
  ) +

  scale_fill_manual(
    values = c(
      Mangroves = "#8C6BB1",
      Seagrass = "#F2E55C"
    )
  ) +

  scale_colour_manual(
    values = c(
      Mangroves = "#8C6BB1",
      Seagrass = "#F2E55C"
    )
  ) +

  xlab(
    paste0(
      "PC1 (",
      variance_pcoa_jaccard_selected[1],
      "%)"
    )
  ) +

  ylab(
    paste0(
      "PC2 (",
      variance_pcoa_jaccard_selected[2],
      "%)"
    )
  )

plot_pcoa_jaccard_selected

## Highlight Desulfobacterales and Chromatiales
plot_pcoa_jaccard_selected_h <- ggplot(
  coords_pcoa_jaccard_selected,
  aes(PC1, PC2)
) +

  stat_ellipse(
    aes(
      fill = location,
      colour = location
    ),
    geom = "polygon",
    alpha = 0.15,
    linewidth = 0.8
  ) +

  geom_point(
    data = subset(
      coords_pcoa_jaccard_selected,
      Highlight == "Other"
    ),
    shape = 24,
    size = 4.2,
    fill = "grey80",
    colour = "black",
    stroke = 0.25
  ) +

  geom_point(
    data = subset(
      coords_pcoa_jaccard_selected,
      Highlight != "Other"
    ),
    aes(fill = Highlight),
    shape = 24,
    size = 4.2,
    colour = "black",
    stroke = 0.25
  ) +

  scale_fill_manual(
    values = c(

      Mangroves = "#8C6BB1",
      Seagrass = "#F2E55C",

      Mangrove_Desulfobacterales = "#5E3C99",
      Mangrove_Chromatiales = "#C2A5CF",

      Seagrass_Desulfobacterales = "#E6AB02",
      Seagrass_Chromatiales = "#FFF176"
    )
  ) +

  scale_colour_manual(
    values = c(
      Mangroves = "#8C6BB1",
      Seagrass = "#F2E55C"
    )
  ) +

  xlab(
    paste0(
      "PC1 (",
      variance_pcoa_jaccard_selected[1],
      "%)"
    )
  ) +

  ylab(
    paste0(
      "PC2 (",
      variance_pcoa_jaccard_selected[2],
      "%)"
    )
  )

plot_pcoa_jaccard_selected_h

## PERMANOVA
perm_pcoa_jaccard_selected <- adonis2(
  dist_jaccard_selected ~ location,
  data = meta_selected,
  permutations = 9999
)

perm_pcoa_jaccard_selected

## PERMDISP
betadis_pcoa_jaccard_selected <- betadisper(
  dist_jaccard_selected,
  meta_selected$location
)

anova(
  betadis_pcoa_jaccard_selected
)

permutest(
  betadis_pcoa_jaccard_selected,
  permutations = 9999
)

## Save figures
ggsave(
  "PCoA_jaccard_selected.png",
  plot = plot_pcoa_jaccard_selected,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

ggsave(
  "PCoA_jaccard_selected_highlight.png",
  plot = plot_pcoa_jaccard_selected_h,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)