rm(list = ls())
library(phyloseq)
library(tidyverse)
library(ggpubr)
library("rstatix")  # add_xy_position

load("Data/mangrove_phyloseq_dada.RData")

# Estimate alpha diversity
ps_filtered <- subset_samples(mangrove_ps, !(Compartment %in% c("Leaves", "Root") & PNA == "NonPNA"))

# Remove taxa that became zero after filtering
ps_filtered <- prune_taxa(taxa_sums(ps_filtered) > 0, ps_filtered)

table(sample_data(ps_filtered)$Compartment, sample_data(ps_filtered)$PNA)

theme0 <- theme(
  axis.title = element_text(size = 12, colour = "black"),
  axis.text = element_text(size = 10, colour = "black"),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  axis.text.x = element_text(angle = 90, hjust = 1),
  panel.grid = element_blank()
)

plot_alpha_diversity <- function(ps_obj,
                                 metric = "Shannon",
                                 compartments = NULL,      # order of all compartments for plotting
                                 present_comps = NULL,     # which compartments to use for pairwise tests
                                 colors = NULL,
                                 outfile = NULL) {
  
  # Compute alpha diversity
  alpha_div <- estimate_richness(
    ps_obj,
    measures = c("Observed", "Chao1", "Shannon", "Simpson")
  )
  
  # Add metadata
  metadata <- sample_data(ps_obj)
  alpha_div <- cbind(Compartment = metadata$Compartment, alpha_div) %>%
    as.data.frame()
  
  # Set compartment order for plotting
  if(is.null(compartments)) {
    compartments <- unique(alpha_div$Compartment)  # default order from metadata
  }
  alpha_div$Compartment <- factor(alpha_div$Compartment, levels = compartments)
  
  # Compute pairs only for the specified present compartments
  if(!is.null(present_comps)) {
    pair_list <- combn(present_comps, 2, simplify = FALSE)
    
    # Wilcoxon tests for each pair
    sig_pairs <- lapply(pair_list, function(pair) {
      vals1 <- alpha_div[[metric]][alpha_div$Compartment == pair[1]]
      vals2 <- alpha_div[[metric]][alpha_div$Compartment == pair[2]]
      p <- wilcox.test(vals1, vals2)$p.value
      if(p < 0.05) return(pair) else return(NULL)
    }) %>% Filter(Negate(is.null), .)
  } else {
    sig_pairs <- NULL
  }
  
  # Default colors
  if(is.null(colors)) {
    colors <- rainbow(length(compartments))
  }
  
  # Plot
  p <- ggplot(alpha_div, aes(x = Compartment, y = .data[[metric]], fill = Compartment)) +
    geom_boxplot(outlier.colour = "transparent", width = 0.6) +
    geom_jitter(aes(colour = Compartment), height = 0, width = 0.15, alpha = 0.5, size = 1) +
    scale_fill_manual(values = colors) +
    scale_colour_manual(values = colors) +
    stat_compare_means(
      comparisons = sig_pairs,  # only for specified present_comps
      label = "p.signif",
      tip.length = 0.015,
      step.increase = 0.05,
      vjust = 1
    ) +
    labs(x = '', y = paste0(metric, " diversity"), fill = 'Compartment') +
    theme_bw() +
    theme0
  print(p)
}

# Example usage
fixed_order <- c("Control Soil","Control Water", "Leaves", "Root", "Rhizo", "Sediment")
present_comps <- c("Leaves", "Root","Rhizo", "Sediment")
col_comps <- c("#3B9AB2", "gray", "#33A02C", "#6A3D9A", "#FDAE61", "#A73030FF")

p_shannon <- plot_alpha_diversity(ps_filtered,
                                  metric = "Shannon",
                                  compartments = fixed_order,
                                  present_comps = present_comps,
                                  colors = col_comps,
                                  outfile = "Figures_DECIPHER/Shannon_Compartment.pdf")



#ggsave("Figures_DECIPHER/Alpha_Shannon_Compartment.pdf", plot = p_shannon, width = 5.5, height = 6)

p_observed <- plot_alpha_diversity(ps_filtered,
                                   metric = "Observed",
                                   compartments = fixed_order,
                                   present_comps = present_comps,
                                   colors = col_comps)

p_observed


#ggsave("Figures_DECIPHER/Alpha_Observed_Compartment.pdf", plot = p_observed, width = 5.5, height = 6)
#================================================================
#Alpha diversity statistics
#================================================================
alpha_df <- estimate_richness(
  ps_filtered,
  measures = c("Observed", "Chao1", "Shannon", "Simpson")
)
rownames(alpha_df) <- sub("^X", "", rownames(alpha_df))
head(alpha_df)

sampledat <- data.frame(sample_data(ps_filtered)) %>%
  select(PNA, Compartment, Site)%>%
  rownames_to_column("SampleID") 

alpha_df <- alpha_df %>% rownames_to_column("SampleID") 

alpha_df <- alpha_df %>% full_join(sampledat)

#Alpha diversity statistics
shannon_lm <- lm(Shannon ~ Compartment , data = alpha_df)
anova(shannon_lm)

# Post-hoc (only if interaction is significant):
library(emmeans)
emmeans(shannon_lm, pairwise ~ Compartment)

alpha_df$logObserved <- log1p(alpha_df$Observed)
#Richness (Observed ASVs – log-transform)
rich_lm <- lm(logObserved ~ Compartment , data = alpha_df)
anova(rich_lm)

#================================================================
### Alpha Diversity PNA vs NON-PNA is leaves and Root
#================================================================
# Subset Root and Leaves

mangrove_ps_pna <- mangrove_ps %>%
  subset_samples(Compartment %in% c("Root","Leaves")) %>%           # remove controls
  prune_taxa(taxa_sums(.) > 0, .) %>%                     # drop zero-abundance taxa
  prune_samples(sample_sums(.) > 0, .) 

alpha_df <- estimate_richness(
  mangrove_ps_pna,
  measures = c("Observed", "Chao1", "Shannon", "Simpson")
)
rownames(alpha_df) <- sub("^X", "", rownames(alpha_df))
head(alpha_df)

sampledat <- data.frame(sample_data(mangrove_ps_pna)) %>%
  select(PNA, Compartment, Site)%>%
  rownames_to_column("SampleID") 

alpha_df <- alpha_df %>% rownames_to_column("SampleID") 

alpha_df <- alpha_df %>% full_join(sampledat)

# Make sure factors are set correctly
alpha_df$Compartment <- factor(alpha_df$Compartment,
                               levels = c("Leaves", "Root"))

alpha_df$PNA <- factor(alpha_df$PNA,
                       levels = c("NonPNA", "PNA"))

plot_alpha <- function(alpha_df,
                       metric,
                       plot_type = c("violin", "boxplot"),
                       scale_type = c("normal", "log10")) {
  
  plot_type  <- match.arg(plot_type)
  scale_type <- match.arg(scale_type)
  
  df <- alpha_df
  
  # -----------------------------
  # Apply transformation if needed
  # -----------------------------
  if (scale_type == "log10") {
    
    if (any(df[[metric]] <= 0, na.rm = TRUE)) {
      stop("Log10 transformation requires strictly positive values.")
    }
    
    df[[metric]] <- log10(df[[metric]])
    y_label <- paste0("log10(", metric, ")")
    
  } else {
    y_label <- metric
  }
  
  # -----------------------------
  # Statistical test
  # -----------------------------
  stat.test <- df %>%
    group_by(Compartment) %>%
    wilcox_test(formula = as.formula(paste(metric, "~ PNA"))) %>%
    add_significance("p") %>%
    mutate(
      p.format = format.pval(p, digits = 3),
      label = paste0(p.signif, "\n", "p = ", p.format)
    ) %>%
    add_xy_position(x = "PNA") %>%
    mutate(
      group1_inter = paste0(group1, ".", Compartment),
      group2_inter = paste0(group2, ".", Compartment)
    )
  
  # -----------------------------
  # Base plot
  # -----------------------------
  p <- ggplot(df,
              aes(x = PNA,
                  y = .data[[metric]],
                  fill = Compartment)) +
    facet_wrap(~Compartment) +
    scale_fill_manual(values = c(Root = "#6A3D9A",
                                 Leaves = "#33A02C")) +
    scale_color_manual(values = c(Root = "#6A3D9A",
                                  Leaves = "#33A02C")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme_bw() +
    xlab("") +
    ylab(y_label)
  
  # -----------------------------
  # Add geometry
  # -----------------------------
  if (plot_type == "violin") {
    
    p <- p +
      geom_violin(adjust = 1.5, alpha =0.7) +
      geom_jitter(aes(color = Compartment),
                  width = 0.15,
                  size = 1.5,
                  show.legend = FALSE)
    
  } else {
    
    p <- p +
      geom_boxplot(outlier.shape = NA, width = 0.7, alpha =0.7) +
      geom_jitter(aes(color = Compartment),
                  width = 0.15,
                  size = 1.5,
                  show.legend = FALSE)
  }
  
  # -----------------------------
  # Add p-values
  # -----------------------------
  p +
    stat_pvalue_manual(stat.test,
                       label = "label",
                       size = 4)
}
plot_alpha(alpha_df, "Shannon", "violin", "normal")

#ggsave("Figures_DECIPHER/Alpha_Shannon_Compartment_PNA_Violin.pdf", width = 5, height = 4.5)

plot_alpha(alpha_df, "Shannon", "boxplot", "normal")
# ggsave("Figures_DECIPHER/Alpha_Shannon_Compartment_PNA_boxplot.pdf", width = 5, height = 4.5)

plot_alpha(alpha_df, "Observed", "boxplot", "log10")
# ggsave("Figures_DECIPHER/Alpha_Observed_Compartment_PNA_boxplot_log.pdf", width = 5, height = 4.5)

plot_alpha(alpha_df, "Observed", "boxplot", "normal")
# ggsave("Figures_DECIPHER/Alpha_Observed_Compartment_PNA_boxplot.pdf", width = 5, height = 4.5)

plot_alpha(alpha_df, "Observed", "violin", "normal")
# ggsave("Figures_DECIPHER/Alpha_Observed_Compartment_PNA_Violin.pdf", width = 5, height = 4.5)

#================================================

#==============================================================
# Check
unique(sample_data(mangrove_ps_pna)$PNA)
#================================================================
#Alpha diversity statistics
#================================================================
alpha_df <- estimate_richness(
  mangrove_ps_pna,
  measures = c("Observed", "Chao1", "Shannon", "Simpson")
)
rownames(alpha_df) <- sub("^X", "", rownames(alpha_df))
head(alpha_df)

sampledat <- data.frame(sample_data(mangrove_ps_pna)) %>%
  select(PNA, Compartment, Site)%>%
  rownames_to_column("SampleID") 

alpha_df <- alpha_df %>% rownames_to_column("SampleID") 

alpha_df <- alpha_df %>% full_join(sampledat)

#Alpha diversity statistics
shannon_lm <- lm(Shannon ~ Compartment * PNA, data = alpha_df)
anova(shannon_lm)

# Post-hoc (only if interaction is significant):
library(emmeans)
emmeans(shannon_lm, pairwise ~ PNA | Compartment)

alpha_df$logObserved <- log1p(alpha_df$Observed)
#Richness (Observed ASVs – log-transform)
rich_lm <- lm(logObserved ~ Compartment * PNA, data = alpha_df)
anova(rich_lm)

#-----------------------------------------------
# 1. Prepare sequencing depth dataframe
#-----------------------------------------------
depth_df <- data.frame(
  Depth = sample_sums(mangrove_ps_pna),
  Compartment = sample_data(mangrove_ps_pna)$Compartment,
  PNA = sample_data(mangrove_ps_pna)$PNA
)

#-----------------------------------------------
# 2. Plot sequencing depth (log-scale) by PNA and Compartment
#-----------------------------------------------
cols <- c(
  "#33A02C",  # e.g., Control or PNA=0
  "#6A3D9A"   # e.g., Treatment or PNA=1
)

p_depth <- depth_df %>%
  filter(Depth > 0) %>%
  ggplot(aes(x = PNA, y = Depth, fill = Compartment)) +
  geom_boxplot(outlier.colour = "transparent", width = 0.6) +
  geom_jitter(aes(colour = Compartment), width = 0.1, height = 0.1, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  scale_y_log10(breaks = c( 100, 1000, 10000)) +
# log10 scale to handle sequencing depth differences
  xlab("") +
  facet_grid(~Compartment)+
  ylab("Sequencing depth") +
  theme_bw() 

p_depth
#-----------------------------------------------
# 3. Save plot
#-----------------------------------------------
#ggsave("Figures_DECIPHER/SamplingDepth_Compartment_PNA.pdf", plot = p_depth, width = 5, height = 3.5)
