##### Functional comparison of Sedimenticolaceae MAGs #####
#####              Mangroves vs Cordgrass             #####

#1. Libraries
#2. Read data
#3. Build KO matrix

#4. GLOBAL ANALYSIS
# Bray–Curtis
# PCoA
# Plots
# PERMANOVA
# PERMDISP
# Save figures

#5. SELECTED GENE ANALYSIS
# Bray–Curtis
# PCoA
# Plots
# PERMANOVA
# PERMDISP
# Save figures

# 6. ENVFIT


### 1. Libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

theme_set(theme_classic(base_family = "ArialMT"))

## 2. Read data
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

###3. Build KO matrix
## Parse KO annotations
ko_table <- ko |>
  filter(str_detect(kegg_ko, "K[0-9]{5}")) |>
  separate_rows(kegg_ko, sep = ",") |>
  mutate(
    kegg_ko = str_extract(kegg_ko, "K[0-9]{5}")
  ) |>
  filter(!is.na(kegg_ko))
## Count number of KO annotations per CDS
ko_table <- ko_table |>
  group_by(query) |>
  mutate(
    n_KO = n(),
    weight = 1 / n_KO
  ) |>
  ungroup()
## Sum weighted counts per genome
ko_counts <- ko_table |>
  group_by(genome, kegg_ko) |>
  summarise(
    abundance = sum(weight),
    .groups = "drop"
  )
## Correct abundances by genome completeness
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
## Sedimenticolaceae genomes
sed_genomes <- metadata |>
  filter(family == "Sedimenticolaceae") |>
  pull(genome)

## GLOBAL KO MATRIX (all KO codes)
ko_matrix_global <- ko_counts |>
  filter(genome %in% sed_genomes) |>
  select(genome, kegg_ko, abundance) |>
  pivot_wider(
    names_from = kegg_ko,
    values_from = abundance,
    values_fill = 0
  ) |>
  as.data.frame()

rownames(ko_matrix_global) <- ko_matrix_global$genome
ko_matrix_global$genome <- NULL

## Metadata for global matrix
meta_global <- metadata |>
  filter(genome %in% rownames(ko_matrix_global)) |>
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
  distinct() |>
  arrange(match(genome, rownames(ko_matrix_global)))

stopifnot(
  all(meta_global$genome == rownames(ko_matrix_global))
)

## SELECTED KO MATRIX (only pathways of interest)
selected_kos <- unique(pathway$ko_code)

ko_matrix_selected <- ko_matrix_global[
  ,
  colnames(ko_matrix_global) %in% selected_kos,
  drop = FALSE
]

## Remove genomes with no selected KOs
ko_matrix_selected <- ko_matrix_selected[
  rowSums(ko_matrix_selected) > 0,
  ,
  drop = FALSE
]

## Remove KOs absent from all genomes
ko_matrix_selected <- ko_matrix_selected[
  ,
  colSums(ko_matrix_selected) > 0,
  drop = FALSE
]

## Metadata for selected matrix
meta_selected <- meta_global |>
  filter(genome %in% rownames(ko_matrix_selected)) |>
  arrange(match(genome, rownames(ko_matrix_selected)))

stopifnot(
  all(meta_selected$genome == rownames(ko_matrix_selected))
)

### 4. GLOBAL ANALYSIS
## Bray–Curtis distance
bc_global <- vegdist(
  ko_matrix_global,
  method = "bray"
)

## PCoA
pcoa_global <- cmdscale(
  bc_global,
  eig = TRUE,
  k = 2
)

coords_global <- as.data.frame(
  pcoa_global$points
)

colnames(coords_global) <- c(
  "PC1",
  "PC2"
)

coords_global <- cbind(
  coords_global,
  meta_global
)

## Percentage of variance explained
eig_global <- pcoa_global$eig

variance_global <- round(
  eig_global[eig_global > 0] /
    sum(eig_global[eig_global > 0]) * 100,
  1
)

## PCoA plot
plot_global <- ggplot(
  coords_global,
  aes(
    PC1,
    PC2
  )
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
    aes(
      fill = location
    ),
    shape = 24,
    size = 4.2,
    colour = "black",
    stroke = 0.35
  ) +

  scale_fill_manual(
    values = c(
      red_sea = "#8C6BB1",
      China = "#F2E55C"
    )
  ) +

  scale_colour_manual(
    values = c(
      red_sea = "#8C6BB1",
      China = "#F2E55C"
    )
  ) +

  xlab(
    paste0(
      "PC1 (",
      variance_global[1],
      "%)"
    )
  ) +

  ylab(
    paste0(
      "PC2 (",
      variance_global[2],
      "%)"
    )
  )

plot_global

## PERMANOVA
permanova_global <- adonis2(
  bc_global ~ location,
  data = meta_global,
  permutations = 9999
)

permanova_global

## PERMDISP
betadisper_global <- betadisper(
  bc_global,
  meta_global$location
)

anova(
  betadisper_global
)

permutest(
  betadisper_global,
  permutations = 9999
)

## Save figures
ggsave(
  "tables/PCoA_Sedimenticolaceae_global.png",
  plot = plot_global,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

### 5. SELECTED GENE ANALYSIS
## Bray–Curtis distance
bc_selected <- vegdist(
  ko_matrix_selected,
  method = "bray"
)

## PCoA
pcoa_selected <- cmdscale(
  bc_selected,
  eig = TRUE,
  k = 2
)

coords_selected <- as.data.frame(
  pcoa_selected$points
)

colnames(coords_selected) <- c(
  "PC1",
  "PC2"
)

coords_selected <- cbind(
  coords_selected,
  meta_selected
)

## Percentage of variance explained
eig_selected <- pcoa_selected$eig

variance_selected <- round(
  eig_selected[eig_selected > 0] /
    sum(eig_selected[eig_selected > 0]) * 100,
  1
)

## PCoA plot
plot_selected <- ggplot(
  coords_selected,
  aes(
    PC1,
    PC2
  )
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
    aes(
      fill = location
    ),
    shape = 24,
    size = 4.2,
    colour = "black",
    stroke = 0.35
  ) +

  scale_fill_manual(
    values = c(
      red_sea = "#8C6BB1",
      China = "#F2E55C"
    )
  ) +

  scale_colour_manual(
    values = c(
      red_sea = "#8C6BB1",
      China = "#F2E55C"
    )
  ) +

  xlab(
    paste0(
      "PC1 (",
      variance_selected[1],
      "%)"
    )
  ) +

  ylab(
    paste0(
      "PC2 (",
      variance_selected[2],
      "%)"
    )
  )

plot_selected

## PERMANOVA
permanova_selected <- adonis2(
  bc_selected ~ location,
  data = meta_selected,
  permutations = 9999
)

permanova_selected

## PERMDISP
betadisper_selected <- betadisper(
  bc_selected,
  meta_selected$location
)

anova(
  betadisper_selected
)

permutest(
  betadisper_selected,
  permutations = 9999
)

## Save figures
ggsave(
  "results/PCoA_Sedimenticolaceae_selected.png",
  plot = plot_selected,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

### 6. ENVFIT
library(ggrepel)

## Fit KO vectors
envfit_selected <- envfit(
  pcoa_selected$points,
  ko_matrix_selected,
  permutations = 9999
)

## Extract vectors
vectors <- scores(
  envfit_selected,
  display = "vectors"
)

envfit_results <- data.frame(
  KO  = rownames(vectors),
  PC1 = vectors[, 1],
  PC2 = vectors[, 2],
  r2  = envfit_selected$vectors$r,
  p   = envfit_selected$vectors$pvals
)

## Add pathway annotations
envfit_results <- envfit_results |>
  left_join(
    pathway,
    by = c("KO" = "ko_code")
  ) |>
  arrange(desc(r2))

## Significant KOs only
envfit_plot <- envfit_results |>
  filter(p < 0.05) |>
  slice_head(n = 10)


## Scale arrows
ord_max <- max(
  abs(coords_selected$PC1),
  abs(coords_selected$PC2)
)

vec_max <- max(
  abs(envfit_plot$PC1),
  abs(envfit_plot$PC2)
)

arrow_scale <- 0.7 * ord_max / vec_max

envfit_plot$PC1 <- envfit_plot$PC1 * arrow_scale
envfit_plot$PC2 <- envfit_plot$PC2 * arrow_scale

## PCoA and ENVFIT
plot_selected_envfit <- ggplot(
  coords_selected,
  aes(
    PC1,
    PC2
  )
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
    aes(
      fill = location
    ),
    shape = 24,
    size = 4,
    colour = "black",
    stroke = 0.35
  ) +

  geom_segment(
    data = envfit_plot,
    aes(
      x = 0,
      y = 0,
      xend = PC1,
      yend = PC2
    ),
    inherit.aes = FALSE,
    linewidth = 0.7,
    colour = "black",
    arrow = arrow(length = unit(0.18, "cm"))
  ) +

  geom_text_repel(
    data = envfit_plot,
    aes(
      x = PC1,
      y = PC2,
      label = gene
    ),
    inherit.aes = FALSE,
    nudge_x = 0.25,
    direction = "y",
    hjust = 0,
    segment.color = "grey40",
    family = "sans",
    fontface = "italic",
    size = 3.5
  ) +

  scale_fill_manual(
    values = c(
      red_sea = "#8C6BB1",
      China = "#F2E55C"
    )
  ) +

  scale_colour_manual(
    values = c(
      red_sea = "#8C6BB1",
      China = "#F2E55C"
    )
  ) +

  labs(
    x = paste0(
      "PC1 (",
      variance_selected[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      variance_selected[2],
      "%)"
    )
  ) +

  theme_classic(base_family = "Arial") +

  theme(
    legend.title = element_blank()
  )

plot_selected_envfit

## Export results
write.csv(
  envfit_results,
  "results/Envfit_significant_KOs.csv",
  row.names = FALSE
)

ggsave(
  "results/PCoA_envfit_top10_KOs.png",
  plot = plot_selected_envfit,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)
