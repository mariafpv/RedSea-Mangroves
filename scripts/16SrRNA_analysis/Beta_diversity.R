library(vegan)
library(phyloseq)
library(tidyverse)

rm(list = ls())
load("Data/mangrove_phyloseq_dada.RData")

# Estimate alpha diversity
ps_filtered <- subset_samples(mangrove_ps, !(Compartment %in% c("Leaves", "Root") & PNA == "NonPNA"))

# Remove taxa that became zero after filtering
ps_filtered <- prune_taxa(taxa_sums(ps_filtered) > 0, ps_filtered)

table(sample_data(ps_filtered)$Compartment, sample_data(ps_filtered)$PNA)
#=============================================================
## Beta Diversity
### Drop control samples (soil + water)
#-----------------------------------------------
# 1. Remove control compartments and zero-abundance taxa/samples
#-----------------------------------------------
controls <- c("Control Soil", "Control Water")

ps_nc <- ps_filtered %>%
  subset_samples(!Compartment %in% controls) %>%           # remove controls
  prune_taxa(taxa_sums(.) > 0, .) %>%                     # drop zero-abundance taxa
  prune_samples(sample_sums(.) > 0, .) %>%
  transform_sample_counts(., function(x) x / sum(x))

#-----------------------------------------------
# 2. Prepare OTU matrix and metadata
#-----------------------------------------------
meta_df<- sample_data(ps_nc) %>%
  data.frame() %>%
  mutate(
    Compartment = factor(Compartment,
                         levels = c("Leaves", "Root", "Rhizo", "Sediment"))
  )

#-----------------------------------------------
# 3. Compute Bray-Curtis distance
#-----------------------------------------------
bray_tax <- phyloseq::distance(ps_nc, method = "bray")

#-----------------------------------------------
# 4. PERMANOVA and dispersion test
#-----------------------------------------------
# Run PERMANOVA
adonis_res <- adonis2(
  bray_tax ~ Compartment,       # formula: distance ~ grouping factor
  data = meta_df,
  permutations = 999            # number of random permutations
)

adonis_tab <- adonis_res %>% as.data.frame() 
adonis_tab$Term <- rownames(adonis_tab)

adonis_tab <- adonis_tab %>% 
  dplyr::select(Term, Df, SumOfSqs, R2, F,`Pr(>F)`)

rownames(adonis_tab) <- NULL

adonis_tab

adonis_main <-adonis_tab %>%
  filter(Term == "Model")

adonis_main
#=============================================
# Compute pairwise_permanova
#==============================================
pairwise_permanova <- function(dist_mat, meta_df, group_var,
                               permutations = 999,
                               p_adjust_method = "BH") {
  
  # Check if grouping variable exists
  if (!group_var %in% colnames(meta_df)) {
    stop("Grouping variable not found in metadata.")
  }
  
  # Convert to factor
  meta_df[[group_var]] <- as.factor(meta_df[[group_var]])
  
  groups <- levels(meta_df[[group_var]])
  
  if (length(groups) < 2) {
    stop("Need at least two groups for pairwise comparisons.")
  }
  
  combs <- combn(groups, 2)
  
  results <- apply(combs, 2, function(x) {
    
    # Subset distance matrix and metadata
    subset_idx <- meta_df[[group_var]] %in% x
    dist_sub <- as.dist(as.matrix(dist_mat)[subset_idx, subset_idx])
    meta_sub <- meta_df[subset_idx, , drop = FALSE]
    
    # Run PERMANOVA
    ad <- adonis2(dist_sub ~ meta_sub[[group_var]], permutations = permutations)
    
    data.frame(
      Group1 = x[1],
      Group2 = x[2],
      R2 = ad$R2[1],
      F = ad$F[1],
      p = ad$`Pr(>F)`[1]
    )
  })
  
  results_df <- bind_rows(results) %>%
    mutate(
      R2 = round(R2, 3),
      F = round(F, 3),
      p = signif(p, 3),
      p_adj = signif(p.adjust(p, method = p_adjust_method), 3)
    ) %>%
    arrange(p_adj)
  
  return(results_df)
}

pairwise_df_permanova <- pairwise_permanova(bray_tax, meta_df, "Compartment")
pairwise_df_permanova
knitr::kable(
  pairwise_df_permanova,
  caption = "Pairwise PERMANOVA comparisons among compartments (Bray–Curtis, 999 permutations, BH-corrected p-values)"
)

#==============================================
#Dispersion
#=============================================

disp <- betadisper(bray_tax, meta_df$Compartment)

# Test significance
anova(disp)        # F-test for dispersion differences
permutest(disp, permutations = 99999)

# Post-hoc pairwise comparisons 
# pairwise_disp only if overall dispersion significant
pairwise_disp <- function(dist_mat, meta_df, group_var,
                          permutations = 99999,
                          p_adjust_method = "BH") {
  
  # Ensure grouping variable exists
  if (!group_var %in% colnames(meta_df)) {
    stop("Grouping variable not found in metadata.")
  }
  
  # Ensure factor
  meta_df[[group_var]] <- as.factor(meta_df[[group_var]])
  
  groups <- levels(meta_df[[group_var]])
  
  if (length(groups) < 2) {
    stop("Need at least two groups.")
  }
  
  combs <- combn(groups, 2)
  
  results <- apply(combs, 2, function(x) {
    
    subset_idx <- meta_df[[group_var]] %in% x
    
    dist_sub <- as.dist(as.matrix(dist_mat)[subset_idx, subset_idx])
    meta_sub <- meta_df[subset_idx, , drop = FALSE]
    
    disp_sub <- betadisper(dist_sub, meta_sub[[group_var]])
    test <- permutest(disp_sub, permutations = permutations)
    
    data.frame(
      Group1 = x[1],
      Group2 = x[2],
      Df = test$tab[1, "Df"],
      SS= test$tab[1, "Sum Sq"],
      MS = test$tab[1, "Mean Sq"],
      F = test$tab[1, "F"],
      NPerm = test$tab[1, "N.Perm"],
      p = test$tab[1, "Pr(>F)"]
    )
  })
  
  results_df <- bind_rows(results) %>%
    mutate(
      p_adj = p.adjust(p, method = p_adjust_method),
      SS = round(SS, 3),
      MS = round(MS, 3),
      F = round(F, 3),
      p = signif(p, 3),
      p_adj = signif(p_adj, 3)
    ) %>%
    arrange(p_adj)
  
  return(results_df)
}

pairwise_disp(bray_tax, meta_df, "Compartment", permutations = 99999)

#-----------------------------------------------
# 5. PCoA
#-----------------------------------------------
ord_tax <- ordinate(ps_nc, method = "PCoA", distance = bray_tax, trymax = 100)
#tax_dist_nc <- vegdist(otu_nc, method = "bray")

# Extract sample coordinates
pcoa_df <- as.data.frame(ord_tax$vectors)

pcoa_df$Sample <- rownames(pcoa_df)
pcoa_df$Compartment <- sample_data(ps_nc)$Compartment

# Rename axes for clarity
colnames(pcoa_df)[1:2] <- c("PCoA1", "PCoA2")

# Fix order of compartments
pcoa_df$Compartment <- factor(
  pcoa_df$Compartment,
  levels = c("Leaves", "Root", "Rhizo", "Sediment")
)

# Percent variance explained
eig_vals <- ord_tax$values$Relative_eig * 100

x_lab <- paste0("PCoA1 (", round(eig_vals[1], 1), "%)")
y_lab <- paste0("PCoA2 (", round(eig_vals[2], 1), "%)")

#-----------------------------------------------
# 6. PCoA Plot
#-----------------------------------------------
comp_colors <- c(
  "#33A02C",  # Leaves
  "#6A3D9A",  # Root
  "#FDAE61",  # Rhizo
  "#A73030FF" # Sediment
)

p_pcoa <- ggplot(pcoa_df,
                 aes(x = PCoA1,
                     y = PCoA2,
                     color = Compartment,
                     shape = Compartment)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(fill = Compartment),
               geom = "polygon",
               level = 0.95,
               alpha = 0.25,
               color = NA) +
  scale_color_manual(values = comp_colors) +
  scale_fill_manual(values = comp_colors) +
  labs(x = x_lab, y = y_lab) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

p_pcoa

# Save plot
# ggsave("Figures_DECIPHER/PCOA_Compartment.pdf", plot = p_pcoa, width = 6, height = 4.5)

#=============================================================
## NMDS
### Drop control samples (soil + water)
#-----------------------------------------------
# 1. NMDS (2 dimensions)
#-----------------------------------------------

nmds_tax <- ordinate(ps_nc, method = "NMDS",
                     distance = bray_tax,
                     trymax = 100)

#-----------------------------------------------
# 2. Prepare NMDS coordinates dataframe
#-----------------------------------------------
nmds_df <- data.frame(
  Sample = rownames(nmds_tax$points),
  Compartment = sample_data(ps_nc)$Compartment,
  NMDS1 = nmds_tax$points[, 1],
  NMDS2 = nmds_tax$points[, 2]
)

# Ensure compartments follow a fixed order
nmds_df$Compartment <- factor(
  nmds_df$Compartment,
  levels = c("Leaves", "Root", "Rhizo", "Sediment")
)

#-----------------------------------------------
# 3. Colors
#-----------------------------------------------
comp_colors <- c(
  "#33A02C",  # Leaves
  "#6A3D9A",  # Root
  "#FDAE61",  # Rhizo
  "#A73030FF" # Sediment
)

#-----------------------------------------------
# 4. NMDS Plot
#-----------------------------------------------
p_nmds <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Compartment, shape = Compartment)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(fill = Compartment), geom = "polygon", level = 0.95, alpha = 0.25, color = NA) +
  scale_color_manual(values = comp_colors) +
  scale_fill_manual(values = comp_colors) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

p_nmds
# Save plot
# ggsave("Figures_DECIPHER/NMDS_Compartment.pdf", plot = p_nmds, width = 5.5, height = 5)

#-----------------------------------------------
# 5. NMDS stress
#-----------------------------------------------
nmds_tax$stress

#=============================================================================
# 3. Beta diversity Leaves and Roots only
# 3.1 Bray–Curtis (relative abundance)
#=============================================================================
# Subset Root and Leaves
ps_pna_rel <- mangrove_ps %>%
  subset_samples(Compartment %in% c("Root","Leaves")) %>%           # remove controls
  prune_taxa(taxa_sums(.) > 0, .) %>%                     # drop zero-abundance taxa
  prune_samples(sample_sums(.) > 0, .) %>%
  transform_sample_counts(., function(x) x / sum(x)) 

meta_df_pna <- data.frame(sample_data(ps_pna_rel))

bray_dist_pna <- phyloseq::distance(ps_pna_rel, method = "bray")

adonis2(bray_dist_pna ~ Compartment * PNA, data = meta_df_pna)


meta_df_pna$Group <- interaction(meta_df_pna$Compartment,
                                 meta_df_pna$PNA)
pairwise_permanova(bray_dist_pna, meta_df_pna, "Group")

#=====================================================
# PERMDISP
#=====================================================
disp <- betadisper(bray_dist_pna, interaction(meta_df_pna$Compartment, meta_df_pna$PNA))
anova(disp)
permutest(disp, permutations = 99999)

pairwise_disp(bray_dist_pna, meta_df_pna, "Group", permutations = 99999)

#-----------------------------------------------------

ord_tax <- ordinate(ps_pna_rel, method = "PCoA", distance = bray_dist_pna, trymax = 100)

# Extract sample coordinates
pcoa_df <- as.data.frame(ord_tax$vectors)

pcoa_df$Sample <- rownames(pcoa_df)
pcoa_df$Compartment <- sample_data(ps_pna_rel)$Compartment
pcoa_df$PNA <- sample_data(ps_pna_rel)$PNA
# Rename axes for clarity
colnames(pcoa_df)[1:2] <- c("PCoA1", "PCoA2")

# Fix order of compartments
pcoa_df$Compartment <- factor( pcoa_df$Compartment, levels = c("Leaves", "Root"))

pcoa_df$PNA <- factor(pcoa_df$PNA, levels = c( "PNA", "NonPNA"))

# Percent variance explained
eig_vals <- ord_tax$values$Relative_eig * 100

x_lab <- paste0("PCoA1 (", round(eig_vals[1], 1), "%)")
y_lab <- paste0("PCoA2 (", round(eig_vals[2], 1), "%)")

#-----------------------------------------------
# 6. PCoA Plot
#-----------------------------------------------

comp_colors <- c("#1A8612", "#6A3D9A", "#B2DF8A", "#FCCDE5")
p_pcoa <- ggplot(pcoa_df,
                 aes(x = PCoA1,
                     y = PCoA2,
                     color = Compartment,
                     shape = PNA)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(fill = interaction(Compartment, PNA)),
               geom = "polygon",
               level = 0.95,
               alpha = 0.25,
               color = NA) +
  scale_color_manual(values = comp_colors) +
  scale_fill_manual(values = comp_colors) +
  labs(x = x_lab, y = y_lab) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

p_pcoa
# Save plot
#ggsave("Figures_DECIPHER/PCOA_Compartment_PNA_ord.pdf", plot = p_pcoa, width = 6, height = 3.5)

#------------------------------------------------------
ord <- ordinate(ps_pna_rel, method = "NMDS", distance = bray_dist_pna, trymax = 100)
# Extract sample scores (coordinates)
nmds_scores <- scores(ord, display = "sites")  # 'sites' gives sample points

# NMDS1 coordinates
nmds <- nmds_scores[, "NMDS1"]

# Optionally, combine with sample metadata
nmds_df <- as.data.frame(nmds_scores)

# View NMDS1
head(nmds_df$NMDS1)

#ps_pna_rel1 <- prune_samples(!sample_names(ps_pna_rel) %in% c("111_16SPNA"), ps_pna_rel)
ps_pna_rel1 <- ps_pna_rel
bray_dist_pna1 <- phyloseq::distance(ps_pna_rel1, method = "bray")

ord1 <- ordinate(ps_pna_rel1, method = "NMDS", distance = bray_dist_pna1 , trymax = 100)

nmds_scores1 <- scores(ord1, display = "sites") 

nmds_df1 <- as.data.frame(nmds_scores1)
nmds_df1$Sample <- rownames(nmds_df1)
nmds_df1$Compartment <- sample_data(ps_pna_rel1)$Compartment
nmds_df1$PNA <- sample_data(ps_pna_rel1)$PNA

col1 <- c("#1A8612", "#6A3D9A", "#B2DF8A", "#FCCDE5")

ggplot(nmds_df1, aes(
  x = NMDS1,
  y = NMDS2,
  color = Compartment,
  shape = PNA
)
) +
  geom_point(size = 2, alpha = 0.85) +
  #stat_chull(aes(fill = interaction(Compartment, PNA)), geom = "polygon", alpha = 0.25, color = NA)+
  stat_ellipse(aes(fill = interaction(Compartment, PNA)), geom = "polygon", level = 0.95, alpha = 0.5, color = NA) +
  scale_color_manual(values = c("#1A8612", "#6A3D9A")) +
  scale_fill_manual(values = col1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

#ggsave("Figures_DECIPHER/NMDS_Compartment_PNA_ord.pdf", width = 6, height = 4)

#===============================================================
# Subset Leaves 
#--------------------------------------------------------
ps_leaf_rel <- mangrove_ps %>%
  subset_samples(Compartment == "Leaves") %>%           # remove controls
  prune_taxa(taxa_sums(.) > 0, .) %>%                     # drop zero-abundance taxa
  prune_samples(sample_sums(.) > 0, .) %>%
  transform_sample_counts(., function(x) x / sum(x)) 

# Remove taxa with zero counts after subsetting

#ps_leaf_rel <- prune_samples(!sample_names(ps_leaf_rel) %in% c("111_16SPNA"), ps_leaf_rel)

# Bray-Curtis distance
bray_leaf <- phyloseq::distance(ps_leaf_rel, method = "bray")

# NMDS
ord_leaf <- ordinate(ps_leaf_rel, method = "NMDS",
                     distance = bray_leaf,
                     trymax = 100)

nmds_scores_leaf <- scores(ord_leaf, display = "sites") 
# Extract coordinates

leaf_df <- as.data.frame(nmds_scores_leaf)
leaf_df$Sample <- rownames(nmds_scores_leaf)
leaf_df$Compartment <- sample_data(ps_leaf_rel)$Compartment
leaf_df$PNA <- sample_data(ps_leaf_rel)$PNA

ggplot(leaf_df, aes(NMDS1, NMDS2, color = PNA, shape = PNA)) +
  geom_point(size = 2.5) +
  stat_ellipse(aes(fill = PNA),
             geom = "polygon",
             alpha = 0.25,
             color = NA) +
  scale_color_manual(values = c("#C0FF3E", "#1A8612")) +
  scale_fill_manual(values = c("#C0FF3E", "#1A8612")) +
  theme_bw() +
  ggtitle("NMDS – Leaf ")

#ggsave("Figures_DECIPHER/NMDS_Compartment_Leaf.pdf", width = 5, height = 3.5)
#===============================================================
# Subset Roots
#--------------------------------------------------------

ps_root_rel <- mangrove_ps %>%
  subset_samples(Compartment == "Root") %>%           # remove controls
  prune_taxa(taxa_sums(.) > 0, .) %>%                     # drop zero-abundance taxa
  prune_samples(sample_sums(.) > 0, .) %>%
  transform_sample_counts(., function(x) x / sum(x)) 

# Bray-Curtis distance
bray_root <- phyloseq::distance(ps_root_rel, method = "bray")

# NMDS
ord_root <- ordinate(ps_root_rel, method = "NMDS", distance = bray_root, trymax = 100)

nmds_scores_root <- scores(ord_root, display = "sites") 
# Extract coordinates

root_df <- as.data.frame(nmds_scores_root)
root_df$Sample <- rownames(nmds_scores_root)
root_df$Compartment <- sample_data(ps_root_rel)$Compartment
root_df$PNA <- sample_data(ps_root_rel)$PNA

ggplot(root_df, aes(NMDS1, NMDS2, color = PNA, shape = PNA)) +
  geom_point(size = 2.5) +
  stat_ellipse(aes(fill = PNA),
               geom = "polygon",
               alpha = 0.25,
               color = NA) +
  scale_color_manual(values = c("#E066FF", "#68228B")) +
  scale_fill_manual(values = c("#E066FF", "#68228B")) +
  theme_bw() +
  ggtitle("NMDS – Root")

#ggsave("Figures_DECIPHER/NMDS_Compartment_root.pdf", width = 5, height = 3.5)

meta_df_root <- sample_data(ps_root_rel)
permanova_root <- adonis2(
  bray_root ~ PNA,
  data = meta_df_root,
  permutations = 99999,
  by = "margin",
  method = "bray"
)
#--------------------------------------------------

meta_df_leaf <- sample_data(ps_leaf_rel) %>% data.frame()
permanova_leaf <- adonis2(
  bray_leaf ~ PNA,
  data = meta_df_leaf,
  permutations = 99999,
  by = "margin",
  method = "bray"
)

ord_leaf_pco <- ordinate(ps_leaf_rel, method = "PCoA", distance = bray_leaf, trymax = 100)

# Extract sample coordinates
pcoa_df_leaf <- as.data.frame(ord_leaf_pco$vectors)

pcoa_df_leaf$Sample <- rownames(pcoa_df_leaf)
pcoa_df_leaf$Compartment <- sample_data(ps_leaf_rel)$Compartment
pcoa_df_leaf$PNA <- sample_data(ps_leaf_rel)$PNA
# Rename axes for clarity
colnames(pcoa_df_leaf)[1:2] <- c("PCoA1", "PCoA2")

pcoa_df_leaf$PNA <- factor(pcoa_df_leaf$PNA, levels = c( "PNA", "NonPNA"))

# Percent variance explained
eig_vals <- pcoa_df_leaf$values$Relative_eig * 100

x_lab <- paste0("PCoA1 (", round(eig_vals[1], 1), "%)")
y_lab <- paste0("PCoA2 (", round(eig_vals[2], 1), "%)")

#-----------------------------------------------
# 6. PCoA Plot Leaves
#-----------------------------------------------

comp_colors <- c("#C0FF3E", "#1A8612")
p_pcoa_leaf <- ggplot(pcoa_df_leaf,
                 aes(x = PCoA1,
                     y = PCoA2,
                     color = PNA,
                     shape = PNA)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(fill =  PNA),
               geom = "polygon",
               level = 0.95,
               alpha = 0.25,
               color = NA) +
  scale_color_manual(values = comp_colors) +
  scale_fill_manual(values = comp_colors) +
  labs(x = x_lab, y = y_lab) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

p_pcoa_leaf

#--------------------------------------------------
ord_root_pco <- ordinate(ps_root_rel, method = "PCoA", distance = bray_root, trymax = 100)

# Extract sample coordinates
pcoa_df_root <- as.data.frame(ord_root_pco$vectors)

pcoa_df_root$Sample <- rownames(pcoa_df_root)
pcoa_df_root$Compartment <- sample_data(ps_root_rel)$Compartment
pcoa_df_root$PNA <- sample_data(ps_root_rel)$PNA
# Rename axes for clarity
colnames(pcoa_df_root)[1:2] <- c("PCoA1", "PCoA2")

# Fix order of compartments

pcoa_df_root$PNA <- factor(pcoa_df_root$PNA, levels = c( "PNA", "NonPNA"))

# Percent variance explained
eig_vals <- pcoa_df_leaf$values$Relative_eig * 100

x_lab <- paste0("PCoA1 (", round(eig_vals[1], 1), "%)")
y_lab <- paste0("PCoA2 (", round(eig_vals[2], 1), "%)")

#-----------------------------------------------
# 6. PCoA Plot Root
#-----------------------------------------------

comp_colors <- c("#E066FF", "#68228B")
p_pcoa_root <- ggplot(pcoa_df_root,
                 aes(x = PCoA1,
                     y = PCoA2,
                     color = PNA,
                     shape = PNA)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(fill =  PNA),
               geom = "polygon",
               level = 0.95,
               alpha = 0.25,
               color = NA) +
  scale_color_manual(values = comp_colors) +
  scale_fill_manual(values = comp_colors) +
  labs(x = x_lab, y = y_lab) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

p_pcoa_root

bray_roots <- phyloseq::distance(ps_root_rel, method = "bray")

