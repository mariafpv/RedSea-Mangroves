library (tidyr)
library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(pheatmap)

load("Data/mangrove_phyloseq_with_chloroplast.RData")

ps_sub1 <- subset_samples(ps_sub, Compartment %in% c("Root", "Leaves"))

# 2. Replace Domain with "Chloroplast" if Order == "Chloroplast"
tax_df <- as.data.frame(tax_table(ps_sub1))
tax_df$Domain <- ifelse(tax_df$Order == "Chloroplast", "Chloroplast", as.character(tax_df$Domain))
tax_df$Domain <- ifelse(tax_df$Family == "Mitochondria", "Mitochondria", as.character(tax_df$Domain))
tax_table(ps_sub1) <- tax_table(as.matrix(tax_df))

# 3. Function to process each facet independently
process_facet <- function(ps_obj, compartment_name, pna_status){
  ps_facet <- prune_samples(
    sample_data(ps_obj)$Compartment == compartment_name &
      sample_data(ps_obj)$PNA == pna_status,
    ps_obj
  )
  
  # Remove samples with zero total abundance
  ps_facet <- prune_samples(sample_sums(ps_facet) > 0, ps_facet)
  
  if(nsamples(ps_facet) == 0) return(NULL)
  
  # Aggregate at Domain level first
  ps_domain <- tax_glom(ps_facet, taxrank = "Domain")
  
  # Re-normalize counts per sample to sum to 1
  ps_rel <- transform_sample_counts(ps_domain, function(x) x / sum(x))
  
  return(ps_rel)
}

# 4. Generate all facets
facet_list <- list(
  Root_NonPNA = process_facet(ps_sub1, "Root", "NonPNA"),
  Root_PNA    = process_facet(ps_sub1, "Root", "PNA"),
  Leaves_NonPNA = process_facet(ps_sub1, "Leaves", "NonPNA"),
  Leaves_PNA    = process_facet(ps_sub1, "Leaves", "PNA")
)

# Remove empty facets
facet_list <- facet_list[!sapply(facet_list, is.null)]

df_long <- bind_rows(lapply(names(facet_list), function(nm){
  df <- psmelt(facet_list[[nm]])
  df$Facet <- nm
  df
}))

# 6. Factorize for plotting
df_long$PNA <- factor(df_long$PNA, levels = c("NonPNA", "PNA"))
df_long$Compartment <- factor(df_long$Compartment, levels = c("Root", "Leaves"))

# 1. Create combined facet label
df_long <- df_long %>%
  mutate(
    FacetRow = PNA,                # Rows = NonPNA / PNA
    FacetCol = Compartment,        # Columns = Leaves / Root
    FacetLabel = paste(FacetRow, FacetCol, sep = " - ")
  )

# 2. Order facets manually for 2x2 grid
facet_order <- c(
  "NonPNA - Leaves",
  "NonPNA - Root",
  "PNA - Leaves",
  "PNA - Root"
)
df_long$FacetLabel <- factor(df_long$FacetLabel, levels = facet_order)

# 3. Re-level SampleID within each facet to remove empty space
df_long <- df_long %>%
  group_by(FacetLabel) %>%
  mutate(SampleID = factor(SampleID, levels = unique(SampleID))) %>%
  ungroup()

# Remove everything after the first underscore in SampleID
df_long <- df_long %>%
  mutate(SampleID = sub("_.*", "", SampleID)) %>%
  filter(!is.na(Abundance))

p1 <- df_long %>% filter (PNA == "NonPNA") %>%
  ggplot(aes(x = SampleID, y = Abundance, fill = Domain)) +
  geom_bar(stat = "identity") +
  ylab("") + xlab ("") +
  facet_grid(PNA ~ Compartment, scales = "free_x") +
  # theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 10)
  ) +
  scale_fill_manual("Kingdom", values= c("#5869C7", "tomato","#09A39A", "#6A3D9A"))

p2 <- df_long %>% filter (PNA == "PNA") %>%
  ggplot(aes(x = SampleID, y = Abundance, fill = Domain)) +
  geom_bar(stat = "identity") +
  ylab("") + xlab ("") +
  facet_grid(PNA ~ Compartment, scales = "free_x") +
  # theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 10)
  ) +
  scale_fill_manual("Kingdom", values= c("#5869C7", "tomato","#09A39A", "#6A3D9A"))

p3 <- gridExtra::grid.arrange(p1, p2)
p3
#ggsave("Figures_DECIPHER/Relative_Abundance_Kindom1.pdf", plot = p3, width = 7.5, height = 5)
