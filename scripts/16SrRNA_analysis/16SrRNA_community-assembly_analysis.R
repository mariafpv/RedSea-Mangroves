# Mangrove project - R codes used for final analyses and figures
# Generated on: 2026-08-25
# Author: Lucas William Mendes
# This file concatenates the R scripts used during the final analyses.
# Original script paths are indicated at the beginning of each section.
# NOTE: Some scripts use environment variables such as PGPG_OUT_DIR and PGPG_REMOVE_MAGS for corrected figure outputs.


==========================================================================================
# SCRIPT: analysis_assembly/community_assembly_analysis.R
==========================================================================================

library(readxl)
library(openxlsx)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

input_file <- "ASV_table_organizado.xlsx"
out_dir <- "analysis_assembly/results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

treatments <- c(
  Sediment = "Sediments",
  Rhizosphere = "Rhizosphere",
  Root = "Roots",
  Leaves = "Leaves"
)

read_treatment <- function(sheet_name, treatment_name) {
  raw <- readxl::read_excel(input_file, sheet = sheet_name)
  names(raw)[1] <- "ASV"
  mat <- raw |>
    as.data.frame(check.names = FALSE) |>
    column_to_rownames("ASV") |>
    t() |>
    as.data.frame(check.names = FALSE)
  mat[] <- lapply(mat, function(x) as.numeric(replace(x, is.na(x), 0)))
  mat <- as.matrix(mat)
  rownames(mat) <- colnames(raw)[-1]
  metadata <- data.frame(
    sample = rownames(mat),
    treatment = treatment_name,
    stringsAsFactors = FALSE
  )
  list(counts = mat, metadata = metadata)
}

parts <- Map(read_treatment, treatments, names(treatments))
counts <- do.call(rbind, lapply(parts, `[[`, "counts"))
metadata <- bind_rows(lapply(parts, `[[`, "metadata"))

stopifnot(identical(rownames(counts), metadata$sample))
counts <- counts[, colSums(counts) > 0, drop = FALSE]
presence <- ifelse(counts > 0, 1, 0)

sample_summary <- metadata |>
  mutate(
    total_reads = rowSums(counts),
    observed_asvs = rowSums(presence),
    shannon = vegan::diversity(counts, index = "shannon")
  )

treatment_alpha_summary <- sample_summary |>
  group_by(treatment) |>
  summarise(
    n_samples = n(),
    mean_total_reads = mean(total_reads),
    sd_total_reads = sd(total_reads),
    mean_observed_asvs = mean(observed_asvs),
    sd_observed_asvs = sd(observed_asvs),
    mean_shannon = mean(shannon),
    sd_shannon = sd(shannon),
    .groups = "drop"
  )

write.csv(sample_summary, file.path(out_dir, "sample_alpha_summary.csv"), row.names = FALSE)
write.csv(treatment_alpha_summary, file.path(out_dir, "treatment_alpha_summary.csv"), row.names = FALSE)

pair_partition <- function(x, y) {
  x <- x > 0
  y <- y > 0
  a <- sum(x & y)
  b <- sum(x & !y)
  c <- sum(!x & y)
  beta_sor <- if ((2 * a + b + c) == 0) NA_real_ else (b + c) / (2 * a + b + c)
  beta_sim <- if ((a + min(b, c)) == 0) NA_real_ else min(b, c) / (a + min(b, c))
  beta_sne <- beta_sor - beta_sim
  beta_jac <- if ((a + b + c) == 0) NA_real_ else (b + c) / (a + b + c)
  beta_jtu <- if ((a + 2 * min(b, c)) == 0) NA_real_ else 2 * min(b, c) / (a + 2 * min(b, c))
  beta_jne <- beta_jac - beta_jtu
  c(
    shared = a,
    unique_first = b,
    unique_second = c,
    beta_sor = beta_sor,
    beta_sim_turnover = beta_sim,
    beta_sne_nestedness = beta_sne,
    beta_jac = beta_jac,
    beta_jtu_turnover = beta_jtu,
    beta_jne_nestedness = beta_jne
  )
}

pairs <- combn(seq_len(nrow(presence)), 2)
sample_beta <- apply(pairs, 2, function(idx) {
  vals <- pair_partition(presence[idx[1], ], presence[idx[2], ])
  c(
    sample_1 = rownames(presence)[idx[1]],
    sample_2 = rownames(presence)[idx[2]],
    treatment_1 = metadata$treatment[idx[1]],
    treatment_2 = metadata$treatment[idx[2]],
    vals
  )
}) |>
  t() |>
  as.data.frame(stringsAsFactors = FALSE)

num_cols <- setdiff(names(sample_beta), c("sample_1", "sample_2", "treatment_1", "treatment_2"))
sample_beta[num_cols] <- lapply(sample_beta[num_cols], as.numeric)
sample_beta <- sample_beta |>
  mutate(comparison = if_else(treatment_1 == treatment_2, "within_treatment", "between_treatments"))

sample_beta_summary <- sample_beta |>
  group_by(treatment_1, treatment_2, comparison) |>
  summarise(
    n_pairs = n(),
    mean_beta_sor = mean(beta_sor, na.rm = TRUE),
    mean_turnover = mean(beta_sim_turnover, na.rm = TRUE),
    mean_nestedness = mean(beta_sne_nestedness, na.rm = TRUE),
    mean_beta_jac = mean(beta_jac, na.rm = TRUE),
    mean_jaccard_turnover = mean(beta_jtu_turnover, na.rm = TRUE),
    mean_jaccard_nestedness = mean(beta_jne_nestedness, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(sample_beta, file.path(out_dir, "sample_pairwise_beta_partition.csv"), row.names = FALSE)
write.csv(sample_beta_summary, file.path(out_dir, "sample_pairwise_beta_summary.csv"), row.names = FALSE)

pooled_by_treatment <- rowsum(counts, group = metadata$treatment)
pooled_presence <- ifelse(pooled_by_treatment > 0, 1, 0)
pooled_pairs <- combn(rownames(pooled_presence), 2)
pooled_beta <- apply(pooled_pairs, 2, function(pair) {
  vals <- pair_partition(pooled_presence[pair[1], ], pooled_presence[pair[2], ])
  c(treatment_1 = pair[1], treatment_2 = pair[2], vals)
}) |>
  t() |>
  as.data.frame(stringsAsFactors = FALSE)
pooled_beta[setdiff(names(pooled_beta), c("treatment_1", "treatment_2"))] <-
  lapply(pooled_beta[setdiff(names(pooled_beta), c("treatment_1", "treatment_2"))], as.numeric)
write.csv(pooled_beta, file.path(out_dir, "pooled_treatment_beta_partition.csv"), row.names = FALSE)

bray <- vegdist(counts, method = "bray")
jaccard <- vegdist(presence, method = "jaccard", binary = TRUE)

permanova_bray <- adonis2(bray ~ treatment, data = metadata, permutations = 999)
permanova_jaccard <- adonis2(jaccard ~ treatment, data = metadata, permutations = 999)
disp_bray <- anova(betadisper(bray, metadata$treatment), permutations = 999)
disp_jaccard <- anova(betadisper(jaccard, metadata$treatment), permutations = 999)

capture.output(permanova_bray, file = file.path(out_dir, "permanova_bray.txt"))
capture.output(permanova_jaccard, file = file.path(out_dir, "permanova_jaccard.txt"))
capture.output(disp_bray, file = file.path(out_dir, "betadisper_bray.txt"))
capture.output(disp_jaccard, file = file.path(out_dir, "betadisper_jaccard.txt"))

set.seed(42)
nmds <- metaMDS(counts, distance = "bray", k = 2, trymax = 100, autotransform = FALSE, trace = 0)
nmds_scores <- as.data.frame(scores(nmds, display = "sites")) |>
  rownames_to_column("sample") |>
  left_join(metadata, by = "sample")
write.csv(nmds_scores, file.path(out_dir, "nmds_bray_scores.csv"), row.names = FALSE)

p_nmds <- ggplot(nmds_scores, aes(NMDS1, NMDS2, color = treatment)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(linewidth = 0.8, level = 0.68, show.legend = FALSE) +
  theme_bw(base_size = 12) +
  labs(
    title = "NMDS - Bray-Curtis",
    subtitle = paste0("Stress = ", round(nmds$stress, 3)),
    x = "NMDS1",
    y = "NMDS2",
    color = "Treatment"
  )
ggsave(file.path(out_dir, "nmds_bray.png"), p_nmds, width = 7, height = 5, dpi = 300)

radfit_one <- function(treatment_name) {
  abund <- colSums(counts[metadata$treatment == treatment_name, , drop = FALSE])
  abund <- sort(abund[abund > 0], decreasing = TRUE)
  invisible(capture.output(
    fit <- suppressWarnings(radfit(abund)),
    type = "message"
  ))
  aic <- AIC(fit)
  if (is.null(dim(aic))) {
    aic_values <- as.numeric(aic)
    model_names <- names(aic)
  } else {
    aic_values <- as.numeric(aic[, "AIC"])
    model_names <- rownames(aic)
  }
  out <- data.frame(
    treatment = treatment_name,
    model = model_names,
    AIC = aic_values,
    delta_AIC = aic_values - min(aic_values, na.rm = TRUE),
    row.names = NULL
  ) |>
    filter(!is.na(AIC)) |>
    arrange(AIC)
  list(abundance = abund, fit = fit, aic = out)
}

radfits <- lapply(unique(metadata$treatment), radfit_one)
names(radfits) <- unique(metadata$treatment)
radfit_aic <- bind_rows(lapply(radfits, `[[`, "aic"))
write.csv(radfit_aic, file.path(out_dir, "radfit_aic_by_treatment.csv"), row.names = FALSE)

rad_abund <- bind_rows(lapply(names(radfits), function(tr) {
  abund <- radfits[[tr]]$abundance
  data.frame(
    treatment = tr,
    rank = seq_along(abund),
    abundance = as.numeric(abund),
    ASV = names(abund),
    stringsAsFactors = FALSE
  )
}))
write.csv(rad_abund, file.path(out_dir, "rank_abundance_by_treatment.csv"), row.names = FALSE)

p_rad <- ggplot(rad_abund, aes(rank, abundance, color = treatment)) +
  geom_line(linewidth = 0.7) +
  scale_y_log10() +
  theme_bw(base_size = 12) +
  labs(
    title = "Rank-abundance por tratamento",
    x = "Rank da ASV",
    y = "Abundancia total agregada (log10)",
    color = "Treatment"
  )
ggsave(file.path(out_dir, "rank_abundance_curves.png"), p_rad, width = 7, height = 5, dpi = 300)

workbook <- createWorkbook()
addWorksheet(workbook, "Alpha summary")
writeData(workbook, "Alpha summary", sample_summary)
addWorksheet(workbook, "Alpha by treatment")
writeData(workbook, "Alpha by treatment", treatment_alpha_summary)
addWorksheet(workbook, "Beta summary")
writeData(workbook, "Beta summary", sample_beta_summary)
addWorksheet(workbook, "Pooled beta")
writeData(workbook, "Pooled beta", pooled_beta)
addWorksheet(workbook, "RAD AIC")
writeData(workbook, "RAD AIC", radfit_aic)
addWorksheet(workbook, "NMDS scores")
writeData(workbook, "NMDS scores", nmds_scores)
saveWorkbook(workbook, file.path(out_dir, "community_assembly_results.xlsx"), overwrite = TRUE)

sink(file.path(out_dir, "analysis_summary.txt"))
cat("Community assembly exploratory analysis\n")
cat("Input:", input_file, "\n")
cat("Samples:", nrow(counts), "\n")
cat("ASVs retained:", ncol(counts), "\n")
cat("NMDS Bray stress:", round(nmds$stress, 4), "\n\n")
cat("PERMANOVA Bray-Curtis\n")
print(permanova_bray)
cat("\nPERMANOVA binary Jaccard\n")
print(permanova_jaccard)
cat("\nBeta dispersion Bray-Curtis\n")
print(disp_bray)
cat("\nBeta dispersion binary Jaccard\n")
print(disp_jaccard)
cat("\nBest RAD model by treatment\n")
print(radfit_aic |> group_by(treatment) |> slice_min(AIC, n = 1, with_ties = FALSE) |> ungroup())
cat("\nAlpha diversity by treatment\n")
print(treatment_alpha_summary)
sink()

message("Done. Results written to: ", normalizePath(out_dir))


==========================================================================================
# SCRIPT: analysis_assembly/taxonomic_null_assembly.R
==========================================================================================

library(openxlsx)
library(vegan)
library(dplyr)
library(tibble)

input_file <- Sys.getenv("ASV_INPUT_FILE", unset = "ASV_table_organizado.xlsx")
out_dir <- Sys.getenv("TAXONOMIC_NULL_OUT_DIR", unset = "analysis_assembly/taxonomic_null_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

treatments <- c(
  Sediment = "Sediments",
  Rhizosphere = "Rhizosphere",
  Root = "Roots",
  Leaves = "Leaves"
)

read_treatment <- function(sheet_name, treatment_name) {
  raw <- openxlsx::read.xlsx(input_file, sheet = sheet_name)
  names(raw)[1] <- "ASV"
  mat <- raw |>
    as.data.frame(check.names = FALSE) |>
    column_to_rownames("ASV") |>
    t() |>
    as.data.frame(check.names = FALSE)
  mat[] <- lapply(mat, function(x) as.numeric(replace(x, is.na(x), 0)))
  mat <- as.matrix(mat)
  rownames(mat) <- colnames(raw)[-1]
  metadata <- data.frame(
    sample = rownames(mat),
    treatment = treatment_name,
    stringsAsFactors = FALSE
  )
  list(counts = mat, metadata = metadata)
}

parts <- Map(read_treatment, treatments, names(treatments))
counts <- do.call(rbind, lapply(parts, `[[`, "counts"))
metadata <- bind_rows(lapply(parts, `[[`, "metadata"))
counts <- counts[, colSums(counts) > 0, drop = FALSE]

make_null_community <- function(counts, taxon_prob) {
  out <- matrix(0, nrow = nrow(counts), ncol = ncol(counts), dimnames = dimnames(counts))
  sample_reads <- rowSums(counts)
  sample_richness <- rowSums(counts > 0)
  regional_abundance <- colSums(counts)

  for (i in seq_len(nrow(counts))) {
    richness <- min(sample_richness[i], ncol(counts))
    selected <- sample(seq_len(ncol(counts)), size = richness, replace = FALSE, prob = taxon_prob)
    abundance_prob <- regional_abundance[selected]
    if (sum(abundance_prob) == 0) {
      abundance_prob <- rep(1, length(selected))
    }
    out[i, selected] <- as.vector(rmultinom(1, sample_reads[i], prob = abundance_prob))
  }
  out
}

set.seed(42)
n_random <- 999
observed_bray <- as.matrix(vegdist(counts, method = "bray"))
pairs <- combn(seq_len(nrow(counts)), 2)
less_than_observed <- numeric(ncol(pairs))
equal_observed <- numeric(ncol(pairs))

taxon_prob <- colSums(counts > 0)
taxon_prob[taxon_prob == 0] <- 0

for (iter in seq_len(n_random)) {
  null_counts <- make_null_community(counts, taxon_prob)
  null_bray <- as.matrix(vegdist(null_counts, method = "bray"))
  null_vals <- null_bray[cbind(pairs[1, ], pairs[2, ])]
  obs_vals <- observed_bray[cbind(pairs[1, ], pairs[2, ])]
  less_than_observed <- less_than_observed + as.numeric(null_vals < obs_vals)
  equal_observed <- equal_observed + as.numeric(null_vals == obs_vals)
}

rc_bray <- ((less_than_observed + 0.5 * equal_observed) / n_random - 0.5) * 2
obs_vals <- observed_bray[cbind(pairs[1, ], pairs[2, ])]

pairwise_rcbray <- data.frame(
  sample_1 = rownames(counts)[pairs[1, ]],
  sample_2 = rownames(counts)[pairs[2, ]],
  treatment_1 = metadata$treatment[pairs[1, ]],
  treatment_2 = metadata$treatment[pairs[2, ]],
  observed_bray = obs_vals,
  rc_bray = rc_bray,
  assembly_category = case_when(
    rc_bray > 0.95 ~ "dispersal_limitation_or_divergent_processes",
    rc_bray < -0.95 ~ "homogenizing_dispersal_or_convergent_processes",
    TRUE ~ "ecological_drift_or_undominated"
  ),
  stringsAsFactors = FALSE
)

rcbray_summary <- pairwise_rcbray |>
  mutate(comparison = if_else(treatment_1 == treatment_2, "within_treatment", "between_treatments")) |>
  group_by(treatment_1, treatment_2, comparison, assembly_category) |>
  summarise(
    n_pairs = n(),
    mean_observed_bray = mean(observed_bray, na.rm = TRUE),
    mean_rc_bray = mean(rc_bray, na.rm = TRUE),
    median_rc_bray = median(rc_bray, na.rm = TRUE),
    .groups = "drop"
  )

rcbray_proportions <- pairwise_rcbray |>
  mutate(comparison = if_else(treatment_1 == treatment_2, "within_treatment", "between_treatments")) |>
  count(treatment_1, treatment_2, comparison, assembly_category, name = "n_pairs") |>
  group_by(treatment_1, treatment_2, comparison) |>
  mutate(proportion = n_pairs / sum(n_pairs)) |>
  ungroup()

write.csv(pairwise_rcbray, file.path(out_dir, "pairwise_rcbray.csv"), row.names = FALSE)
write.csv(rcbray_summary, file.path(out_dir, "rcbray_summary.csv"), row.names = FALSE)
write.csv(rcbray_proportions, file.path(out_dir, "rcbray_assembly_proportions.csv"), row.names = FALSE)

workbook <- createWorkbook()
addWorksheet(workbook, "Pairwise RCbray")
writeData(workbook, "Pairwise RCbray", pairwise_rcbray)
addWorksheet(workbook, "Summary")
writeData(workbook, "Summary", rcbray_summary)
addWorksheet(workbook, "Proportions")
writeData(workbook, "Proportions", rcbray_proportions)
saveWorkbook(workbook, file.path(out_dir, "taxonomic_null_assembly_results.xlsx"), overwrite = TRUE)

message("Done. Results written to: ", normalizePath(out_dir))


==========================================================================================
# SCRIPT: analysis_assembly/clam_pairwise_analysis.R
==========================================================================================

library(openxlsx)
library(vegan)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)

input_file <- Sys.getenv("ASV_INPUT_FILE", unset = "ASV_table_organizado.xlsx")
out_dir <- Sys.getenv("CLAM_OUT_DIR", unset = "analysis_assembly/clam_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

treatments <- c(
  Sediment = "Sediments",
  Rhizosphere = "Rhizosphere",
  Roots = "Roots",
  Leaves = "Leaves"
)

treatment_colors <- c(
  Sediment = "#963C35",
  Rhizosphere = "#EAB378",
  Roots = "#613F8F",
  Leaves = "#549B40"
)

read_treatment <- function(sheet_name, treatment_name) {
  raw <- openxlsx::read.xlsx(input_file, sheet = sheet_name)
  names(raw)[1] <- "ASV"
  mat <- raw |>
    as.data.frame(check.names = FALSE) |>
    column_to_rownames("ASV") |>
    t() |>
    as.data.frame(check.names = FALSE)
  mat[] <- lapply(mat, function(x) as.numeric(replace(x, is.na(x), 0)))
  mat <- as.matrix(mat)
  rownames(mat) <- colnames(raw)[-1]
  metadata <- data.frame(
    sample = rownames(mat),
    treatment = treatment_name,
    stringsAsFactors = FALSE
  )
  list(counts = mat, metadata = metadata)
}

parts <- Map(read_treatment, treatments, names(treatments))
counts <- do.call(rbind, lapply(parts, `[[`, "counts"))
metadata <- bind_rows(lapply(parts, `[[`, "metadata"))
stopifnot(identical(rownames(counts), metadata$sample))
counts <- counts[, colSums(counts) > 0, drop = FALSE]

run_clam_pair <- function(pair) {
  treatment_a <- pair[1]
  treatment_b <- pair[2]
  keep <- metadata$treatment %in% c(treatment_a, treatment_b)
  comm <- counts[keep, , drop = FALSE]
  comm <- comm[, colSums(comm) > 0, drop = FALSE]
  groups <- factor(metadata$treatment[keep], levels = c(treatment_a, treatment_b))

  clam <- clamtest(
    comm,
    groups,
    coverage.limit = 10,
    specialization = 2 / 3,
    npoints = 20,
    alpha = 0.05 / 20
  )

  raw_clam <- as.data.frame(clam)
  names(raw_clam)[2:3] <- c("count_a", "count_b")

  df <- raw_clam |>
    rename(
      ASV = Species,
      clam_class = Classes
    ) |>
    mutate(
      treatment_a = treatment_a,
      treatment_b = treatment_b,
      comparison = paste(treatment_a, treatment_b, sep = " vs "),
      assigned_group = case_when(
        clam_class == paste0("Specialist_", treatment_a) ~ treatment_a,
        clam_class == paste0("Specialist_", treatment_b) ~ treatment_b,
        clam_class == "Specialist_A" ~ treatment_a,
        clam_class == "Specialist_B" ~ treatment_b,
        clam_class == "Generalist" ~ "Generalist",
        TRUE ~ "Too rare"
      ),
      plot_color = case_when(
        assigned_group %in% names(treatment_colors) ~ treatment_colors[assigned_group],
        assigned_group == "Generalist" ~ "#555555",
        TRUE ~ "#BDBDBD"
      )
    )

  attr(df, "coverage") <- attr(clam, "coverage")
  df
}

pairs <- combn(names(treatments), 2, simplify = FALSE)
pair_results <- lapply(pairs, run_clam_pair)
clam_results <- bind_rows(pair_results)

clam_summary <- clam_results |>
  count(comparison, treatment_a, treatment_b, assigned_group, clam_class, name = "n_asvs") |>
  group_by(comparison) |>
  mutate(proportion = n_asvs / sum(n_asvs)) |>
  ungroup()

specialist_summary <- clam_results |>
  filter(assigned_group %in% names(treatment_colors)) |>
  count(comparison, assigned_group, name = "n_specialist_asvs") |>
  arrange(comparison, assigned_group)

coverage_summary <- bind_rows(lapply(seq_along(pairs), function(i) {
  pair <- pairs[[i]]
  comparison <- paste(pair[1], pair[2], sep = " vs ")
  cov <- attr(pair_results[[i]], "coverage")
  data.frame(
    comparison = comparison,
    treatment = names(cov),
    coverage = as.numeric(cov),
    stringsAsFactors = FALSE
  )
}))

write.csv(clam_results, file.path(out_dir, "clam_pairwise_results.csv"), row.names = FALSE)
write.csv(clam_summary, file.path(out_dir, "clam_pairwise_summary.csv"), row.names = FALSE)
write.csv(specialist_summary, file.path(out_dir, "clam_specialist_summary.csv"), row.names = FALSE)
write.csv(coverage_summary, file.path(out_dir, "clam_coverage_summary.csv"), row.names = FALSE)

plot_pair <- function(pair) {
  treatment_a <- pair[1]
  treatment_b <- pair[2]
  comparison <- paste(treatment_a, treatment_b, sep = " vs ")
  group_levels <- c(treatment_a, treatment_b, "Generalist", "Too rare")
  df <- clam_results |>
    filter(.data$comparison == comparison) |>
    mutate(
      x = count_a + 1,
      y = count_b + 1,
      assigned_group = factor(
        assigned_group,
        levels = group_levels
      )
    )

  color_values <- c(
    setNames(treatment_colors[treatment_a], treatment_a),
    setNames(treatment_colors[treatment_b], treatment_b),
    Generalist = "#555555",
    `Too rare` = "#BDBDBD"
  )

  ggplot(df, aes(x = x, y = y, color = assigned_group)) +
    geom_point(alpha = 0.72, size = 1.4) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = color_values, drop = FALSE) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = 10)
    ) +
    labs(
      title = comparison,
      x = paste0(treatment_a, " abundance + 1"),
      y = paste0(treatment_b, " abundance + 1")
    )
}

plots <- lapply(pairs, plot_pair)
panel <- wrap_plots(plots, ncol = 3, guides = "collect") +
  plot_annotation(
    title = "CLAM test pairwise classification of ASVs across treatments",
    subtitle = "Specialists are colored by treatment; generalists and too-rare ASVs are shown in gray."
  ) &
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(out_dir, "clam_pairwise_panel.png"), panel, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "clam_pairwise_panel.pdf"), panel, width = 14, height = 8)

for (i in seq_along(pairs)) {
  comparison_file <- gsub(" ", "_", gsub(" vs ", "_vs_", paste(pairs[[i]], collapse = " vs ")))
  ggsave(
    file.path(out_dir, paste0("clam_", comparison_file, ".png")),
    plots[[i]],
    width = 6,
    height = 5,
    dpi = 300
  )
}

workbook <- createWorkbook()
addWorksheet(workbook, "Pairwise CLAM")
writeData(workbook, "Pairwise CLAM", clam_results)
addWorksheet(workbook, "Summary")
writeData(workbook, "Summary", clam_summary)
addWorksheet(workbook, "Specialists")
writeData(workbook, "Specialists", specialist_summary)
addWorksheet(workbook, "Coverage")
writeData(workbook, "Coverage", coverage_summary)
saveWorkbook(workbook, file.path(out_dir, "clam_pairwise_results.xlsx"), overwrite = TRUE)

caption <- paste(
  "Figure X. Pairwise CLAM test classification of ASVs across sediment, rhizosphere, root, and leaf bacterial communities.",
  "Each panel compares two treatments using pooled ASV abundances and classifies ASVs as specialists of either treatment, generalists, or too rare to classify.",
  "Treatment specialists are shown using the specified treatment colors: Sediment (#963C35), Rhizosphere (#EAB378), Roots (#613F8F), and Leaves (#549B40).",
  "Axes represent total ASV abundance in each treatment plus one, shown on a log10 scale.",
  sep = " "
)
writeLines(caption, file.path(out_dir, "clam_pairwise_figure_caption.txt"))

message("Done. Results written to: ", normalizePath(out_dir))


==========================================================================================
# SCRIPT: analysis_assembly/beta_nti_analysis.R
==========================================================================================

library(readxl)
library(openxlsx)
library(ape)
library(vegan)
library(dplyr)
library(tibble)

input_file <- "ASV_table_organizado.xlsx"
tree_file <- "asv_tree.nwk"
out_dir <- "analysis_assembly/beta_nti_results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(tree_file)) {
  stop(
    "Missing phylogenetic tree: ", tree_file, "\n",
    "Place a Newick tree in the project folder with tip labels matching the ASV IDs."
  )
}

treatments <- c(
  Sediment = "Sediments",
  Rhizosphere = "Rhizosphere",
  Root = "Roots",
  Leaves = "Leaves"
)

read_treatment <- function(sheet_name, treatment_name) {
  raw <- readxl::read_excel(input_file, sheet = sheet_name)
  names(raw)[1] <- "ASV"
  mat <- raw |>
    as.data.frame(check.names = FALSE) |>
    column_to_rownames("ASV") |>
    t() |>
    as.data.frame(check.names = FALSE)
  mat[] <- lapply(mat, function(x) as.numeric(replace(x, is.na(x), 0)))
  mat <- as.matrix(mat)
  rownames(mat) <- colnames(raw)[-1]
  metadata <- data.frame(
    sample = rownames(mat),
    treatment = treatment_name,
    stringsAsFactors = FALSE
  )
  list(counts = mat, metadata = metadata)
}

parts <- Map(read_treatment, treatments, names(treatments))
counts <- do.call(rbind, lapply(parts, `[[`, "counts"))
metadata <- bind_rows(lapply(parts, `[[`, "metadata"))
counts <- counts[, colSums(counts) > 0, drop = FALSE]

tree <- read.tree(tree_file)
shared_asvs <- intersect(colnames(counts), tree$tip.label)
if (length(shared_asvs) < 2) {
  stop("Fewer than 2 ASVs are shared between the count table and tree tip labels.")
}

missing_from_tree <- setdiff(colnames(counts), tree$tip.label)
missing_from_counts <- setdiff(tree$tip.label, colnames(counts))

counts <- counts[, shared_asvs, drop = FALSE]
tree <- keep.tip(tree, shared_asvs)
phylo_dist <- cophenetic.phylo(tree)
phylo_dist <- phylo_dist[colnames(counts), colnames(counts)]

nearest_between <- function(present_a, present_b, dist_mat) {
  taxa_a <- names(present_a)[present_a > 0]
  taxa_b <- names(present_b)[present_b > 0]
  if (length(taxa_a) == 0 || length(taxa_b) == 0) {
    return(NA_real_)
  }
  mean(apply(dist_mat[taxa_a, taxa_b, drop = FALSE], 1, min), na.rm = TRUE)
}

beta_mntd_pair <- function(sample_a, sample_b, dist_mat) {
  pa <- sample_a > 0
  pb <- sample_b > 0
  mntd_ab <- nearest_between(pa, pb, dist_mat)
  mntd_ba <- nearest_between(pb, pa, dist_mat)
  mean(c(mntd_ab, mntd_ba), na.rm = TRUE)
}

set.seed(42)
n_random <- 999
pairs <- combn(seq_len(nrow(counts)), 2)

results <- vector("list", ncol(pairs))
for (i in seq_len(ncol(pairs))) {
  idx <- pairs[, i]
  obs <- beta_mntd_pair(counts[idx[1], ], counts[idx[2], ], phylo_dist)
  null_vals <- numeric(n_random)
  for (j in seq_len(n_random)) {
    shuffled <- phylo_dist
    rn <- sample(rownames(shuffled))
    cn <- sample(colnames(shuffled))
    rownames(shuffled) <- rn
    colnames(shuffled) <- cn
    shuffled <- shuffled[colnames(counts), colnames(counts)]
    null_vals[j] <- beta_mntd_pair(counts[idx[1], ], counts[idx[2], ], shuffled)
  }
  beta_nti <- (obs - mean(null_vals, na.rm = TRUE)) / sd(null_vals, na.rm = TRUE)
  results[[i]] <- data.frame(
    sample_1 = rownames(counts)[idx[1]],
    sample_2 = rownames(counts)[idx[2]],
    treatment_1 = metadata$treatment[idx[1]],
    treatment_2 = metadata$treatment[idx[2]],
    observed_beta_mntd = obs,
    null_mean_beta_mntd = mean(null_vals, na.rm = TRUE),
    null_sd_beta_mntd = sd(null_vals, na.rm = TRUE),
    beta_nti = beta_nti,
    assembly_category = case_when(
      beta_nti > 2 ~ "variable_selection",
      beta_nti < -2 ~ "homogeneous_selection",
      TRUE ~ "stochastic_or_undominated"
    ),
    stringsAsFactors = FALSE
  )
}

beta_nti <- bind_rows(results)
beta_nti_summary <- beta_nti |>
  mutate(comparison = if_else(treatment_1 == treatment_2, "within_treatment", "between_treatments")) |>
  group_by(treatment_1, treatment_2, comparison, assembly_category) |>
  summarise(
    n_pairs = n(),
    mean_beta_nti = mean(beta_nti, na.rm = TRUE),
    median_beta_nti = median(beta_nti, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(beta_nti, file.path(out_dir, "pairwise_beta_nti.csv"), row.names = FALSE)
write.csv(beta_nti_summary, file.path(out_dir, "beta_nti_summary.csv"), row.names = FALSE)
write.csv(
  data.frame(
    missing_from_tree = c(missing_from_tree, rep(NA, max(0, length(missing_from_counts) - length(missing_from_tree)))),
    missing_from_counts = c(missing_from_counts, rep(NA, max(0, length(missing_from_tree) - length(missing_from_counts))))
  ),
  file.path(out_dir, "asv_tree_matching_report.csv"),
  row.names = FALSE
)

workbook <- createWorkbook()
addWorksheet(workbook, "Pairwise betaNTI")
writeData(workbook, "Pairwise betaNTI", beta_nti)
addWorksheet(workbook, "Summary")
writeData(workbook, "Summary", beta_nti_summary)
saveWorkbook(workbook, file.path(out_dir, "beta_nti_results.xlsx"), overwrite = TRUE)

message("Done. Results written to: ", normalizePath(out_dir))


==========================================================================================
# SCRIPT: analysis_assembly/run_bnti_from_asv_sequences.R
==========================================================================================

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(openxlsx)
  library(Biostrings)
  library(DECIPHER)
  library(ape)
  library(phangorn)
  library(ggplot2)
  library(scales)
  library(matrixStats)
  library(parallel)
})

set.seed(123)

base_dir <- "/Users/lwmendes/Documents/Codex/Mangrove"
asv_file <- file.path(base_dir, "ASV_table_organizado_sem_mitochondria.xlsx")
seq_file <- file.path(base_dir, "ASV_sequence_mapping.csv")
out_dir <- file.path(base_dir, "analysis_assembly", "bnti_results_no_mitochondria")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_rand <- as.integer(Sys.getenv("BNTI_NRAND", "999"))
n_cores <- as.integer(Sys.getenv("BNTI_CORES", "4"))

message("Reading ASV table and sequence mapping...")
comm_raw <- read.xlsx(asv_file, sheet = "Frequency", colNames = TRUE, check.names = FALSE)
names(comm_raw)[1] <- "ASV"
comm_raw$ASV <- as.character(comm_raw$ASV)

treatment_sheet_map <- c(
  "Sediments" = "Sediment",
  "Rhizosphere" = "Rhizosphere",
  "Roots" = "Roots",
  "Leaves" = "Leaves"
)
treatment_samples <- do.call(
  rbind,
  lapply(names(treatment_sheet_map), function(sheet_name) {
    sheet_data <- read.xlsx(asv_file, sheet = sheet_name, colNames = TRUE, rows = 1:2, check.names = FALSE)
    data.frame(
      Sample = names(sheet_data)[-1],
      Treatment = unname(treatment_sheet_map[sheet_name]),
      stringsAsFactors = FALSE
    )
  })
)

count_cols <- setdiff(names(comm_raw), "ASV")
count_cols <- intersect(count_cols, treatment_samples$Sample)
comm_counts <- as.matrix(comm_raw[, count_cols, drop = FALSE])
storage.mode(comm_counts) <- "numeric"
rownames(comm_counts) <- comm_raw$ASV
comm_counts <- comm_counts[rowSums(comm_counts, na.rm = TRUE) > 0, , drop = FALSE]
comm_counts[is.na(comm_counts)] <- 0

seq_map <- read.csv(seq_file, stringsAsFactors = FALSE, check.names = FALSE)
names(seq_map)[1:2] <- c("ASV", "Sequence")
seq_map$ASV <- as.character(seq_map$ASV)
seq_map$Sequence <- toupper(as.character(seq_map$Sequence))

common_asvs <- intersect(rownames(comm_counts), seq_map$ASV)
missing_sequences <- setdiff(rownames(comm_counts), seq_map$ASV)
if (length(missing_sequences) > 0) {
  write.csv(
    data.frame(ASV = missing_sequences),
    file.path(out_dir, "asvs_missing_sequences.csv"),
    row.names = FALSE
  )
}

comm_counts <- comm_counts[common_asvs, , drop = FALSE]
seq_map <- seq_map[match(common_asvs, seq_map$ASV), , drop = FALSE]

sample_totals <- colSums(comm_counts)
comm_counts <- comm_counts[, sample_totals > 0, drop = FALSE]

message("ASVs retained: ", nrow(comm_counts))
message("Samples retained: ", ncol(comm_counts))

metadata <- data.frame(
  Sample = colnames(comm_counts),
  Treatment = treatment_samples$Treatment[match(colnames(comm_counts), treatment_samples$Sample)],
  stringsAsFactors = FALSE
)
metadata$Treatment <- factor(metadata$Treatment, levels = c("Sediment", "Rhizosphere", "Roots", "Leaves"))

write.csv(metadata, file.path(out_dir, "bnti_sample_metadata.csv"), row.names = FALSE)

aligned_fasta <- file.path(out_dir, "asv_sequences_aligned.fasta")
tree_file <- file.path(out_dir, "asv_upgma_tree.nwk")

if (!file.exists(tree_file)) {
  if (file.exists(aligned_fasta)) {
    message("Using cached alignment: ", aligned_fasta)
    aligned <- readDNAStringSet(aligned_fasta)
  } else {
    message("Aligning ASV sequences with DECIPHER...")
    dna <- DNAStringSet(seq_map$Sequence)
    names(dna) <- seq_map$ASV
    aligned <- AlignSeqs(dna, processors = 4, verbose = TRUE)
    writeXStringSet(aligned, aligned_fasta)
  }

  message("Building UPGMA tree from aligned-sequence distances...")
  dna_dist <- DistanceMatrix(
    aligned,
    method = "overlap",
    type = "dist",
    includeTerminalGaps = FALSE,
    correction = "none",
    processors = 4,
    verbose = TRUE
  )
  tree <- as.phylo(hclust(dna_dist, method = "average"))
  tree <- ladderize(tree)
  tree$edge.length[is.na(tree$edge.length) | tree$edge.length < 0] <- 0
  write.tree(tree, file = tree_file)
} else {
  message("Using cached tree: ", tree_file)
  tree <- read.tree(tree_file)
}

message("Preparing phylogenetic distance matrix...")
tree <- keep.tip(tree, intersect(tree$tip.label, rownames(comm_counts)))
comm_counts <- comm_counts[tree$tip.label, , drop = FALSE]
phylo_dist <- cophenetic.phylo(tree)
phylo_dist <- phylo_dist[rownames(comm_counts), rownames(comm_counts)]

comm <- t(comm_counts)
sample_totals <- rowSums(comm)
rel_abund <- sweep(comm, 1, sample_totals, "/")
pres_list <- lapply(seq_len(nrow(comm)), function(i) which(comm[i, ] > 0))
weight_list <- lapply(seq_len(nrow(comm)), function(i) rel_abund[i, pres_list[[i]]])
sample_names <- rownames(comm)

pairs <- combn(seq_len(nrow(comm)), 2)
pair_df <- data.frame(
  Sample1 = sample_names[pairs[1, ]],
  Sample2 = sample_names[pairs[2, ]],
  Treatment1 = as.character(metadata$Treatment[match(sample_names[pairs[1, ]], metadata$Sample)]),
  Treatment2 = as.character(metadata$Treatment[match(sample_names[pairs[2, ]], metadata$Sample)]),
  stringsAsFactors = FALSE
)
pair_df$Comparison <- ifelse(
  pair_df$Treatment1 <= pair_df$Treatment2,
  paste(pair_df$Treatment1, pair_df$Treatment2, sep = " vs "),
  paste(pair_df$Treatment2, pair_df$Treatment1, sep = " vs ")
)

calc_beta_mntd <- function(dist_mat, pairs_mat, pres, weights, perm = NULL) {
  out <- numeric(ncol(pairs_mat))
  for (k in seq_len(ncol(pairs_mat))) {
    i <- pairs_mat[1, k]
    j <- pairs_mat[2, k]
    ia <- pres[[i]]
    ib <- pres[[j]]
    if (!is.null(perm)) {
      ia_d <- perm[ia]
      ib_d <- perm[ib]
    } else {
      ia_d <- ia
      ib_d <- ib
    }
    dsub <- dist_mat[ia_d, ib_d, drop = FALSE]
    m1 <- rowMins(dsub)
    m2 <- colMins(dsub)
    out[k] <- (sum(m1 * weights[[i]]) + sum(m2 * weights[[j]])) / 2
  }
  out
}

message("Calculating observed betaMNTD...")
obs_betamntd <- calc_beta_mntd(phylo_dist, pairs, pres_list, weight_list)

message("Running taxa-label null model with ", n_rand, " randomizations...")
message("Using ", n_cores, " core(s) for null randomizations.")
null_mean <- numeric(length(obs_betamntd))
null_m2 <- numeric(length(obs_betamntd))
rand_done <- 0L
while (rand_done < n_rand) {
  batch_size <- min(n_cores, n_rand - rand_done)
  seeds <- sample.int(.Machine$integer.max, batch_size)
  batch <- mclapply(
    seeds,
    function(seed) {
      set.seed(seed)
      perm <- sample.int(ncol(comm))
      calc_beta_mntd(phylo_dist, pairs, pres_list, weight_list, perm = perm)
    },
    mc.cores = batch_size
  )
  for (rand_betamntd in batch) {
    rand_done <- rand_done + 1L
    delta <- rand_betamntd - null_mean
    null_mean <- null_mean + delta / rand_done
    null_m2 <- null_m2 + delta * (rand_betamntd - null_mean)
  }
  if (rand_done %% 5 == 0 || rand_done == n_rand) {
    message("  randomization ", rand_done, "/", n_rand)
  }
}
null_sd <- sqrt(null_m2 / (n_rand - 1))
bnti <- (obs_betamntd - null_mean) / null_sd

pair_df$Observed_betaMNTD <- obs_betamntd
pair_df$Null_mean_betaMNTD <- null_mean
pair_df$Null_sd_betaMNTD <- null_sd
pair_df$betaNTI <- bnti
pair_df$Assembly_process <- ifelse(
  pair_df$betaNTI > 2,
  "Variable selection",
  ifelse(pair_df$betaNTI < -2, "Homogeneous selection", "Stochastic/undominated")
)

write.csv(pair_df, file.path(out_dir, "pairwise_bnti_results.csv"), row.names = FALSE)

summary_df <- do.call(
  rbind,
  lapply(split(pair_df, pair_df$Comparison), function(x) {
    data.frame(
      Comparison = unique(x$Comparison),
      n_pairs = nrow(x),
      mean_betaNTI = mean(x$betaNTI, na.rm = TRUE),
      median_betaNTI = median(x$betaNTI, na.rm = TRUE),
      sd_betaNTI = sd(x$betaNTI, na.rm = TRUE),
      variable_selection_pct = mean(x$betaNTI > 2, na.rm = TRUE) * 100,
      homogeneous_selection_pct = mean(x$betaNTI < -2, na.rm = TRUE) * 100,
      stochastic_undominated_pct = mean(abs(x$betaNTI) <= 2, na.rm = TRUE) * 100
    )
  })
)
summary_df <- summary_df[order(summary_df$Comparison), ]
write.csv(summary_df, file.path(out_dir, "bnti_summary_by_comparison.csv"), row.names = FALSE)

process_summary <- as.data.frame.matrix(table(pair_df$Comparison, pair_df$Assembly_process))
process_summary$Comparison <- rownames(process_summary)
process_long <- reshape(
  process_summary,
  varying = setdiff(names(process_summary), "Comparison"),
  v.names = "n",
  timevar = "Assembly_process",
  times = setdiff(names(process_summary), "Comparison"),
  direction = "long"
)
rownames(process_long) <- NULL
process_long$Percent <- ave(process_long$n, process_long$Comparison, FUN = function(x) x / sum(x) * 100)
write.csv(process_long, file.path(out_dir, "bnti_assembly_process_proportions.csv"), row.names = FALSE)

treatment_colors <- c(
  "Sediment" = "#963C35",
  "Rhizosphere" = "#EAB378",
  "Roots" = "#613F8F",
  "Leaves" = "#549B40"
)
process_colors <- c(
  "Variable selection" = "#8C2D04",
  "Homogeneous selection" = "#2B6CB0",
  "Stochastic/undominated" = "#8A8A8A"
)

p1 <- ggplot(pair_df, aes(x = Comparison, y = betaNTI)) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "grey35", linewidth = 0.4) +
  geom_boxplot(outlier.shape = NA, fill = "grey92", color = "grey25", linewidth = 0.35) +
  geom_jitter(aes(color = Treatment1), width = 0.18, alpha = 0.35, size = 0.8, show.legend = FALSE) +
  scale_color_manual(values = treatment_colors) +
  labs(x = NULL, y = expression(beta*"NTI"), title = expression("Pairwise "*beta*"NTI by treatment comparison")) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold"))

p2 <- ggplot(process_long, aes(x = Comparison, y = Percent, fill = Assembly_process)) +
  geom_col(color = "white", linewidth = 0.2) +
  scale_fill_manual(values = process_colors, name = "Assembly process") +
  scale_y_continuous(labels = percent_format(scale = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Pairwise comparisons (%)", title = "Assembly process classification") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "bnti_boxplot_by_comparison.png"), p1, width = 9.5, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "bnti_boxplot_by_comparison.pdf"), p1, width = 9.5, height = 5.5)
ggsave(file.path(out_dir, "bnti_assembly_process_proportions.png"), p2, width = 9.5, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "bnti_assembly_process_proportions.pdf"), p2, width = 9.5, height = 5.5)

if (requireNamespace("patchwork", quietly = TRUE)) {
  suppressPackageStartupMessages(library(patchwork))
  panel <- p1 / p2 + plot_annotation(tag_levels = "A")
  ggsave(file.path(out_dir, "bnti_combined_panel.png"), panel, width = 10, height = 9.5, dpi = 300)
  ggsave(file.path(out_dir, "bnti_combined_panel.pdf"), panel, width = 10, height = 9.5)
}

wb <- createWorkbook()
addWorksheet(wb, "Pairwise_bNTI")
writeData(wb, "Pairwise_bNTI", pair_df)
addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_df)
addWorksheet(wb, "Process_proportions")
writeData(wb, "Process_proportions", process_long)
addWorksheet(wb, "Metadata")
writeData(wb, "Metadata", metadata)
saveWorkbook(wb, file.path(out_dir, "bnti_results.xlsx"), overwrite = TRUE)

message("Done. Outputs written to: ", out_dir)


==========================================================================================
# SCRIPT: analysis_assembly/plot_taxonomic_assembly.R
==========================================================================================

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

in_dir <- Sys.getenv("TAXONOMIC_NULL_IN_DIR", unset = "analysis_assembly/taxonomic_null_results")
out_dir <- in_dir

pairwise <- read.csv(file.path(in_dir, "pairwise_rcbray.csv"), check.names = FALSE)
props <- read.csv(file.path(in_dir, "rcbray_assembly_proportions.csv"), check.names = FALSE)

category_labels <- c(
  dispersal_limitation_or_divergent_processes = "Divergent / dispersal limitation",
  ecological_drift_or_undominated = "Drift / undominated",
  homogenizing_dispersal_or_convergent_processes = "Homogenizing / convergent"
)

category_colors <- c(
  "Divergent / dispersal limitation" = "#C44E52",
  "Drift / undominated" = "#4C72B0",
  "Homogenizing / convergent" = "#55A868"
)

props <- props |>
  mutate(
    pair = paste(treatment_1, treatment_2, sep = " x "),
    assembly_label = recode(assembly_category, !!!category_labels)
  )

pairwise <- pairwise |>
  mutate(
    pair = paste(treatment_1, treatment_2, sep = " x "),
    comparison = if_else(treatment_1 == treatment_2, "Within treatment", "Between treatments"),
    assembly_label = recode(assembly_category, !!!category_labels)
  )

p_props <- ggplot(props, aes(x = reorder(pair, proportion), y = proportion, fill = assembly_label)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = category_colors) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  labs(
    title = "Community assembly processes inferred by RCbray",
    subtitle = "Taxonomic abundance-based null model; thresholds: RCbray > 0.95 or < -0.95",
    x = NULL,
    y = "Proportion of sample pairs"
  )

ggsave(file.path(out_dir, "rcbray_assembly_proportions.png"), p_props, width = 9, height = 6, dpi = 300)

p_box <- ggplot(pairwise, aes(x = reorder(pair, rc_bray, median), y = rc_bray, fill = comparison)) +
  geom_hline(yintercept = c(-0.95, 0.95), linetype = "dashed", color = "gray35") +
  geom_boxplot(outlier.alpha = 0.35, width = 0.65) +
  coord_flip() +
  scale_fill_manual(values = c("Within treatment" = "#8172B3", "Between treatments" = "#CCB974")) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  labs(
    title = "RCbray distribution by comparison",
    subtitle = "Values near zero suggest stochastic or undominated processes",
    x = NULL,
    y = "RCbray"
  )

ggsave(file.path(out_dir, "rcbray_boxplot_by_comparison.png"), p_box, width = 9, height = 6, dpi = 300)

heat <- pairwise |>
  group_by(treatment_1, treatment_2) |>
  summarise(mean_rc_bray = mean(rc_bray, na.rm = TRUE), .groups = "drop")

p_heat <- ggplot(heat, aes(treatment_1, treatment_2, fill = mean_rc_bray)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = round(mean_rc_bray, 2)), size = 4) +
  scale_fill_gradient2(low = "#55A868", mid = "white", high = "#C44E52", midpoint = 0, limits = c(-1, 1)) +
  coord_equal() +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Mean RCbray between treatments",
    fill = "Mean\nRCbray"
  )

ggsave(file.path(out_dir, "rcbray_mean_heatmap.png"), p_heat, width = 6, height = 5, dpi = 300)

message("Figures written to: ", normalizePath(out_dir))


==========================================================================================
# SCRIPT: analysis_assembly/plot_taxonomic_assembly_panel.R
==========================================================================================

library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)

in_dir <- Sys.getenv("TAXONOMIC_NULL_IN_DIR", unset = "analysis_assembly/taxonomic_null_results")
out_dir <- in_dir

pairwise <- read.csv(file.path(in_dir, "pairwise_rcbray.csv"), check.names = FALSE)
props <- read.csv(file.path(in_dir, "rcbray_assembly_proportions.csv"), check.names = FALSE)

category_labels <- c(
  dispersal_limitation_or_divergent_processes = "Divergent / dispersal limitation",
  ecological_drift_or_undominated = "Drift / undominated",
  homogenizing_dispersal_or_convergent_processes = "Homogenizing / convergent"
)

category_colors <- c(
  "Divergent / dispersal limitation" = "#C44E52",
  "Drift / undominated" = "#4C72B0",
  "Homogenizing / convergent" = "#55A868"
)

props <- props |>
  mutate(
    pair = paste(treatment_1, treatment_2, sep = " x "),
    assembly_label = recode(assembly_category, !!!category_labels)
  )

pairwise <- pairwise |>
  mutate(
    pair = paste(treatment_1, treatment_2, sep = " x "),
    comparison = if_else(treatment_1 == treatment_2, "Within treatment", "Between treatments"),
    assembly_label = recode(assembly_category, !!!category_labels)
  )

pair_order <- props |>
  group_by(pair) |>
  summarise(divergent_prop = sum(proportion[assembly_label == "Divergent / dispersal limitation"]), .groups = "drop") |>
  arrange(divergent_prop) |>
  pull(pair)

props$pair <- factor(props$pair, levels = pair_order)
pairwise$pair <- factor(pairwise$pair, levels = pair_order)

theme_panel <- theme_bw(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

p_props <- ggplot(props, aes(x = pair, y = proportion, fill = assembly_label)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = category_colors) +
  theme_panel +
  labs(
    title = "Assembly process proportions",
    x = NULL,
    y = "Proportion of sample pairs"
  )

p_box <- ggplot(pairwise, aes(x = pair, y = rc_bray, fill = comparison)) +
  geom_hline(yintercept = c(-0.95, 0.95), linetype = "dashed", color = "gray35", linewidth = 0.4) +
  geom_boxplot(outlier.alpha = 0.35, width = 0.65, linewidth = 0.35) +
  coord_flip() +
  scale_fill_manual(values = c("Within treatment" = "#8172B3", "Between treatments" = "#CCB974")) +
  theme_panel +
  labs(
    title = "RCbray distributions",
    x = NULL,
    y = "RCbray"
  )

heat <- pairwise |>
  group_by(treatment_1, treatment_2) |>
  summarise(mean_rc_bray = mean(rc_bray, na.rm = TRUE), .groups = "drop")

p_heat <- ggplot(heat, aes(treatment_1, treatment_2, fill = mean_rc_bray)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = round(mean_rc_bray, 2)), size = 3.2) +
  scale_fill_gradient2(low = "#55A868", mid = "white", high = "#C44E52", midpoint = 0, limits = c(-1, 1)) +
  coord_equal() +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold", size = 11),
    legend.position = "bottom"
  ) +
  labs(
    title = "Mean RCbray",
    fill = "Mean RCbray"
  )

panel <- (p_props | p_box | p_heat) +
  plot_annotation(
    tag_levels = "A",
    title = "Taxonomic null-model inference of bacterial community assembly",
    subtitle = "RCbray was estimated using a taxonomic abundance-based null model with 999 randomizations."
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "rcbray_assembly_single_panel.png"), panel, width = 16, height = 6, dpi = 300)
ggsave(file.path(out_dir, "rcbray_assembly_single_panel.pdf"), panel, width = 16, height = 6)

caption <- paste(
  "Figure X. Taxonomic null-model inference of bacterial community assembly across plant and sediment compartments.",
  "(A) Proportion of sample pairs assigned to divergent/dispersal-limited, stochastic/undominated, or homogenizing/convergent assembly categories based on RCbray thresholds.",
  "(B) Distribution of pairwise RCbray values for each within- and between-treatment comparison; dashed lines indicate the classification thresholds at RCbray = -0.95 and 0.95.",
  "(C) Mean RCbray values among treatments. Positive values indicate communities that are more dissimilar than expected under the null model, whereas negative values indicate communities that are more similar than expected.",
  sep = " "
)
writeLines(caption, file.path(out_dir, "rcbray_assembly_single_panel_caption.txt"))

message("Single-panel figure written to: ", normalizePath(out_dir))


==========================================================================================
# SCRIPT: analysis_assembly/plot_leaves_rcbray_bnti.R
==========================================================================================

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(openxlsx)
})

base_dir <- "/Users/lwmendes/Documents/Codex/Mangrove"
rc_file <- file.path(base_dir, "analysis_assembly", "taxonomic_null_results_no_mitochondria", "pairwise_rcbray.csv")
bnti_file <- file.path(base_dir, "analysis_assembly", "bnti_results_no_mitochondria", "pairwise_bnti_results.csv")
out_dir <- file.path(base_dir, "analysis_assembly", "leaves_rcbray_bnti_panel")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rc <- read.csv(rc_file, stringsAsFactors = FALSE)
bnti <- read.csv(bnti_file, stringsAsFactors = FALSE)

make_pair_id <- function(a, b) {
  ifelse(a <= b, paste(a, b, sep = "__"), paste(b, a, sep = "__"))
}

rc_leaves <- subset(rc, treatment_1 == "Leaves" & treatment_2 == "Leaves")
rc_leaves$Pair_ID <- make_pair_id(rc_leaves$sample_1, rc_leaves$sample_2)
rc_leaves$RCbray_process <- ifelse(
  rc_leaves$rc_bray > 0.95,
  "Dispersal limitation/divergent",
  ifelse(
    rc_leaves$rc_bray < -0.95,
    "Homogenizing dispersal/convergent",
    "Drift/undominated"
  )
)

bnti_leaves <- subset(bnti, Treatment1 == "Leaves" & Treatment2 == "Leaves")
bnti_leaves$Pair_ID <- make_pair_id(bnti_leaves$Sample1, bnti_leaves$Sample2)

merged <- merge(
  rc_leaves[, c("Pair_ID", "sample_1", "sample_2", "observed_bray", "rc_bray", "RCbray_process")],
  bnti_leaves[, c("Pair_ID", "Observed_betaMNTD", "Null_mean_betaMNTD", "Null_sd_betaMNTD", "betaNTI", "Assembly_process")],
  by = "Pair_ID",
  all = FALSE
)

merged$Integrated_process <- ifelse(
  merged$betaNTI > 2,
  "Variable selection",
  ifelse(
    merged$betaNTI < -2,
    "Homogeneous selection",
    ifelse(
      merged$rc_bray > 0.95,
      "Dispersal limitation",
      ifelse(merged$rc_bray < -0.95, "Homogenizing dispersal", "Drift/undominated")
    )
  )
)

process_levels <- c(
  "Homogeneous selection",
  "Variable selection",
  "Homogenizing dispersal",
  "Dispersal limitation",
  "Drift/undominated",
  "RCbray homogenizing/convergent",
  "RCbray dispersal limitation/divergent"
)

scatter_colors <- c(
  "Homogeneous selection" = "#2B6CB0",
  "Variable selection" = "#8C2D04",
  "Homogenizing dispersal" = "#549B40",
  "Dispersal limitation" = "#963C35",
  "Drift/undominated" = "#8A8A8A"
)

bar_colors <- c(
  "Homogeneous selection" = "#2B6CB0",
  "Variable selection" = "#8C2D04",
  "Homogenizing dispersal" = "#549B40",
  "Dispersal limitation" = "#963C35",
  "Drift/undominated" = "#8A8A8A",
  "RCbray homogenizing/convergent" = "#549B40",
  "RCbray dispersal limitation/divergent" = "#963C35"
)

rc_process_display <- ifelse(
  merged$RCbray_process == "Homogenizing dispersal/convergent",
  "RCbray homogenizing/convergent",
  ifelse(
    merged$RCbray_process == "Dispersal limitation/divergent",
    "RCbray dispersal limitation/divergent",
    merged$RCbray_process
  )
)
rc_bar <- data.frame(Method = "RCbray only", Process = rc_process_display, stringsAsFactors = FALSE)
integrated_bar <- data.frame(Method = "βNTI + RCbray", Process = merged$Integrated_process, stringsAsFactors = FALSE)
bar_df <- rbind(rc_bar, integrated_bar)
bar_df$Method <- factor(bar_df$Method, levels = c("RCbray only", "βNTI + RCbray"))
bar_df$Process <- factor(bar_df$Process, levels = process_levels)
bar_summary <- as.data.frame(table(bar_df$Method, bar_df$Process), stringsAsFactors = FALSE)
names(bar_summary) <- c("Method", "Process", "n")
bar_summary <- subset(bar_summary, n > 0)
bar_summary$Percent <- ave(bar_summary$n, bar_summary$Method, FUN = function(x) x / sum(x) * 100)

summary_df <- data.frame(
  Metric = c(
    "Number of leaf pairs",
    "Mean βNTI",
    "Median βNTI",
    "Mean RCbray",
    "Median RCbray",
    "βNTI homogeneous selection (%)",
    "RCbray homogenizing dispersal/convergent (%)",
    "Integrated homogeneous selection (%)",
    "Integrated homogenizing dispersal (%)"
  ),
  Value = c(
    nrow(merged),
    mean(merged$betaNTI),
    median(merged$betaNTI),
    mean(merged$rc_bray),
    median(merged$rc_bray),
    mean(merged$betaNTI < -2) * 100,
    mean(merged$rc_bray < -0.95) * 100,
    mean(merged$Integrated_process == "Homogeneous selection") * 100,
    mean(merged$Integrated_process == "Homogenizing dispersal") * 100
  )
)

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.tag = element_text(face = "bold", size = base_size + 2)
    )
}

p_scatter <- ggplot(merged, aes(x = betaNTI, y = rc_bray, color = Integrated_process)) +
  annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = "#2B6CB0", alpha = 0.06) +
  annotate("rect", xmin = -2, xmax = 2, ymin = -Inf, ymax = -0.95, fill = "#549B40", alpha = 0.08) +
  annotate("rect", xmin = -2, xmax = 2, ymin = 0.95, ymax = Inf, fill = "#963C35", alpha = 0.08) +
  geom_hline(yintercept = c(-0.95, 0.95), linetype = "dashed", color = "grey35", linewidth = 0.35) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey35", linewidth = 0.35) +
  geom_point(size = 2.2, alpha = 0.82) +
  scale_color_manual(values = scatter_colors, name = "Integrated process") +
  scale_x_continuous(breaks = pretty_breaks(6)) +
  scale_y_continuous(breaks = pretty_breaks(6)) +
  coord_cartesian(ylim = c(-1.05, 1.05), xlim = c(min(merged$betaNTI) - 0.2, 2.35)) +
  labs(
    x = "βNTI",
    y = "RCbray",
    title = "Leaf pairwise βNTI and RCbray"
  ) +
  theme_pub(10) +
  theme(legend.position = "none")

p_bar <- ggplot(bar_summary, aes(x = Method, y = Percent, fill = Process)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.25) +
  geom_text(
    aes(label = ifelse(Percent >= 6, paste0(round(Percent, 1), "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "white"
  ) +
  scale_fill_manual(values = bar_colors, name = "Process") +
  guides(fill = guide_legend(ncol = 1, title.position = "top")) +
  scale_y_continuous(labels = percent_format(scale = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = NULL,
    y = "Leaf pairwise comparisons (%)",
    title = "Assembly classification in leaves"
  ) +
  theme_pub(10) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold")
  )

panel <- p_scatter | p_bar
panel <- panel + plot_layout(widths = c(1.12, 1)) + plot_annotation(tag_levels = "A")

ggsave(file.path(out_dir, "leaves_rcbray_bnti_panel.png"), panel, width = 12.5, height = 5.9, dpi = 300)
ggsave(file.path(out_dir, "leaves_rcbray_bnti_panel.pdf"), panel, width = 12.5, height = 5.9)
ggsave(file.path(out_dir, "leaves_rcbray_bnti_scatter.png"), p_scatter, width = 7, height = 5.8, dpi = 300)
ggsave(file.path(out_dir, "leaves_rcbray_bnti_scatter.pdf"), p_scatter, width = 7, height = 5.8)

write.csv(merged, file.path(out_dir, "leaves_pairwise_rcbray_bnti.csv"), row.names = FALSE)
write.csv(bar_summary, file.path(out_dir, "leaves_rcbray_bnti_process_summary.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "leaves_rcbray_bnti_numeric_summary.csv"), row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "Pairwise")
writeData(wb, "Pairwise", merged)
addWorksheet(wb, "Process_summary")
writeData(wb, "Process_summary", bar_summary)
addWorksheet(wb, "Numeric_summary")
writeData(wb, "Numeric_summary", summary_df)
saveWorkbook(wb, file.path(out_dir, "leaves_rcbray_bnti_results.xlsx"), overwrite = TRUE)

cat("Wrote outputs to:", out_dir, "\n")
print(summary_df)
print(bar_summary)


==========================================================================================
# SCRIPT: analysis_assembly/plot_roots_rcbray_bnti.R
==========================================================================================

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(openxlsx)
})

base_dir <- "/Users/lwmendes/Documents/Codex/Mangrove"
rc_file <- file.path(base_dir, "analysis_assembly", "taxonomic_null_results_no_mitochondria", "pairwise_rcbray.csv")
bnti_file <- file.path(base_dir, "analysis_assembly", "bnti_results_no_mitochondria", "pairwise_bnti_results.csv")
out_dir <- file.path(base_dir, "analysis_assembly", "roots_rcbray_bnti_panel")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rc <- read.csv(rc_file, stringsAsFactors = FALSE)
bnti <- read.csv(bnti_file, stringsAsFactors = FALSE)

make_pair_id <- function(a, b) {
  ifelse(a <= b, paste(a, b, sep = "__"), paste(b, a, sep = "__"))
}

rc_roots <- subset(rc, treatment_1 == "Root" & treatment_2 == "Root")
rc_roots$Pair_ID <- make_pair_id(rc_roots$sample_1, rc_roots$sample_2)
rc_roots$RCbray_process <- ifelse(
  rc_roots$rc_bray > 0.95,
  "Dispersal limitation/divergent",
  ifelse(
    rc_roots$rc_bray < -0.95,
    "Homogenizing dispersal/convergent",
    "Drift/undominated"
  )
)

bnti_roots <- subset(bnti, Treatment1 == "Roots" & Treatment2 == "Roots")
bnti_roots$Pair_ID <- make_pair_id(bnti_roots$Sample1, bnti_roots$Sample2)

merged <- merge(
  rc_roots[, c("Pair_ID", "sample_1", "sample_2", "observed_bray", "rc_bray", "RCbray_process")],
  bnti_roots[, c("Pair_ID", "Observed_betaMNTD", "Null_mean_betaMNTD", "Null_sd_betaMNTD", "betaNTI", "Assembly_process")],
  by = "Pair_ID",
  all = FALSE
)

if (nrow(merged) == 0) stop("No matching root pairwise comparisons were found.")

merged$Integrated_process <- ifelse(
  merged$betaNTI > 2,
  "Variable selection",
  ifelse(
    merged$betaNTI < -2,
    "Homogeneous selection",
    ifelse(
      merged$rc_bray > 0.95,
      "Dispersal limitation",
      ifelse(merged$rc_bray < -0.95, "Homogenizing dispersal", "Drift/undominated")
    )
  )
)

process_levels <- c(
  "Homogeneous selection",
  "Variable selection",
  "Homogenizing dispersal",
  "Dispersal limitation",
  "Drift/undominated",
  "RCbray homogenizing/convergent",
  "RCbray dispersal limitation/divergent"
)

scatter_colors <- c(
  "Homogeneous selection" = "#2B6CB0",
  "Variable selection" = "#8C2D04",
  "Homogenizing dispersal" = "#549B40",
  "Dispersal limitation" = "#963C35",
  "Drift/undominated" = "#8A8A8A"
)

bar_colors <- c(
  "Homogeneous selection" = "#2B6CB0",
  "Variable selection" = "#8C2D04",
  "Homogenizing dispersal" = "#549B40",
  "Dispersal limitation" = "#963C35",
  "Drift/undominated" = "#8A8A8A",
  "RCbray homogenizing/convergent" = "#549B40",
  "RCbray dispersal limitation/divergent" = "#963C35"
)

rc_process_display <- ifelse(
  merged$RCbray_process == "Homogenizing dispersal/convergent",
  "RCbray homogenizing/convergent",
  ifelse(
    merged$RCbray_process == "Dispersal limitation/divergent",
    "RCbray dispersal limitation/divergent",
    merged$RCbray_process
  )
)

bar_df <- rbind(
  data.frame(Method = "RCbray only", Process = rc_process_display, stringsAsFactors = FALSE),
  data.frame(Method = "βNTI + RCbray", Process = merged$Integrated_process, stringsAsFactors = FALSE)
)
bar_df$Method <- factor(bar_df$Method, levels = c("RCbray only", "βNTI + RCbray"))
bar_df$Process <- factor(bar_df$Process, levels = process_levels)
bar_summary <- as.data.frame(table(bar_df$Method, bar_df$Process), stringsAsFactors = FALSE)
names(bar_summary) <- c("Method", "Process", "n")
bar_summary <- subset(bar_summary, n > 0)
bar_summary$Percent <- ave(bar_summary$n, bar_summary$Method, FUN = function(x) x / sum(x) * 100)

summary_df <- data.frame(
  Metric = c(
    "Number of root pairs",
    "Mean βNTI",
    "Median βNTI",
    "Mean RCbray",
    "Median RCbray",
    "βNTI homogeneous selection (%)",
    "RCbray dispersal limitation/divergent (%)",
    "Integrated homogeneous selection (%)",
    "Integrated dispersal limitation (%)",
    "Integrated drift/undominated (%)"
  ),
  Value = c(
    nrow(merged),
    mean(merged$betaNTI),
    median(merged$betaNTI),
    mean(merged$rc_bray),
    median(merged$rc_bray),
    mean(merged$betaNTI < -2) * 100,
    mean(merged$rc_bray > 0.95) * 100,
    mean(merged$Integrated_process == "Homogeneous selection") * 100,
    mean(merged$Integrated_process == "Dispersal limitation") * 100,
    mean(merged$Integrated_process == "Drift/undominated") * 100
  )
)

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.tag = element_text(face = "bold", size = base_size + 2)
    )
}

p_scatter <- ggplot(merged, aes(x = betaNTI, y = rc_bray, color = Integrated_process)) +
  annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = "#2B6CB0", alpha = 0.06) +
  annotate("rect", xmin = -2, xmax = 2, ymin = -Inf, ymax = -0.95, fill = "#549B40", alpha = 0.08) +
  annotate("rect", xmin = -2, xmax = 2, ymin = 0.95, ymax = Inf, fill = "#963C35", alpha = 0.08) +
  geom_hline(yintercept = c(-0.95, 0.95), linetype = "dashed", color = "grey35", linewidth = 0.35) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey35", linewidth = 0.35) +
  geom_point(size = 2.2, alpha = 0.82) +
  scale_color_manual(values = scatter_colors, name = "Integrated process") +
  scale_x_continuous(breaks = pretty_breaks(6)) +
  scale_y_continuous(breaks = pretty_breaks(6)) +
  coord_cartesian(ylim = c(-1.05, 1.05), xlim = c(min(merged$betaNTI) - 0.2, 2.35)) +
  labs(
    x = "βNTI",
    y = "RCbray",
    title = "Root pairwise βNTI and RCbray"
  ) +
  theme_pub(10) +
  theme(legend.position = "none")

p_bar <- ggplot(bar_summary, aes(x = Method, y = Percent, fill = Process)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.25) +
  geom_text(
    aes(label = ifelse(Percent >= 6, paste0(round(Percent, 1), "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "white"
  ) +
  scale_fill_manual(values = bar_colors, name = "Process") +
  guides(fill = guide_legend(ncol = 1, title.position = "top")) +
  scale_y_continuous(labels = percent_format(scale = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = NULL,
    y = "Root pairwise comparisons (%)",
    title = "Assembly classification in roots"
  ) +
  theme_pub(10) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold")
  )

panel <- p_scatter | p_bar
panel <- panel + plot_layout(widths = c(1.12, 1)) + plot_annotation(tag_levels = "A")

ggsave(file.path(out_dir, "roots_rcbray_bnti_panel.png"), panel, width = 12.5, height = 5.9, dpi = 300)
ggsave(file.path(out_dir, "roots_rcbray_bnti_panel.pdf"), panel, width = 12.5, height = 5.9)
ggsave(file.path(out_dir, "roots_rcbray_bnti_scatter.png"), p_scatter, width = 7, height = 5.8, dpi = 300)
ggsave(file.path(out_dir, "roots_rcbray_bnti_scatter.pdf"), p_scatter, width = 7, height = 5.8)

write.csv(merged, file.path(out_dir, "roots_pairwise_rcbray_bnti.csv"), row.names = FALSE)
write.csv(bar_summary, file.path(out_dir, "roots_rcbray_bnti_process_summary.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "roots_rcbray_bnti_numeric_summary.csv"), row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "Pairwise")
writeData(wb, "Pairwise", merged)
addWorksheet(wb, "Process_summary")
writeData(wb, "Process_summary", bar_summary)
addWorksheet(wb, "Numeric_summary")
writeData(wb, "Numeric_summary", summary_df)
saveWorkbook(wb, file.path(out_dir, "roots_rcbray_bnti_results.xlsx"), overwrite = TRUE)

cat("Wrote outputs to:", out_dir, "\n")
print(summary_df)
print(bar_summary)


==========================================================================================
# SCRIPT: PGPg_finder/plot_pgpg_heatmap.R
==========================================================================================

library(pheatmap)
library(RColorBrewer)

input_file <- "PGPg_finder/tables/summary/normalized_summary_table.txt"
out_dir <- "PGPg_finder/figures/codex_heatmaps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(
  input_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

mat <- as.matrix(dat)
storage.mode(mat) <- "numeric"
mat[is.na(mat)] <- 0

group <- ifelse(
  grepl("^Leaves", colnames(mat)),
  "Leaves",
  ifelse(grepl("^Roots", colnames(mat)), "Roots", "Other")
)

annotation_col <- data.frame(Compartment = factor(group, levels = c("Leaves", "Roots")))
rownames(annotation_col) <- colnames(mat)

annotation_colors <- list(
  Compartment = c(
    Leaves = "#549B40",
    Roots = "#613F8F"
  )
)

row_zscore <- t(apply(mat, 1, function(x) {
  sd_x <- sd(x)
  if (is.na(sd_x) || sd_x == 0) {
    rep(0, length(x))
  } else {
    (x - mean(x)) / sd_x
  }
}))
rownames(row_zscore) <- rownames(mat)
colnames(row_zscore) <- colnames(mat)
row_zscore[row_zscore > 2.5] <- 2.5
row_zscore[row_zscore < -2.5] <- -2.5

heat_colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)

png(
  file.path(out_dir, "pgpg_normalized_heatmap_by_mag.png"),
  width = 4200,
  height = 2400,
  res = 300
)
pheatmap(
  row_zscore,
  color = heat_colors,
  breaks = seq(-2.5, 2.5, length.out = 102),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  fontsize_row = 8,
  border_color = NA,
  main = "Normalized plant growth-promoting functional profiles",
  angle_col = 45,
  legend_breaks = c(-2, 0, 2),
  legend_labels = c("Lower", "Mean", "Higher")
)
dev.off()

pdf(
  file.path(out_dir, "pgpg_normalized_heatmap_by_mag.pdf"),
  width = 14,
  height = 8
)
pheatmap(
  row_zscore,
  color = heat_colors,
  breaks = seq(-2.5, 2.5, length.out = 102),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  fontsize_row = 8,
  border_color = NA,
  main = "Normalized plant growth-promoting functional profiles",
  legend_breaks = c(-2, 0, 2),
  legend_labels = c("Lower", "Mean", "Higher")
)
dev.off()

cluster_within_group <- function(matrix_data, group_name) {
  selected <- which(group == group_name)
  if (length(selected) <= 1) {
    return(selected)
  }
  local_matrix <- matrix_data[, selected, drop = FALSE]
  local_tree <- hclust(dist(t(local_matrix)), method = "complete")
  selected[local_tree$order]
}

leaves_order <- cluster_within_group(row_zscore, "Leaves")
roots_order <- cluster_within_group(row_zscore, "Roots")
other_order <- cluster_within_group(row_zscore, "Other")
separated_order <- c(leaves_order, roots_order, other_order)
separated_order <- separated_order[!is.na(separated_order)]

separated_mat <- row_zscore[, separated_order, drop = FALSE]
separated_annotation <- annotation_col[colnames(separated_mat), , drop = FALSE]
group_sizes <- table(factor(group[separated_order], levels = c("Leaves", "Roots", "Other")))
gaps <- cumsum(as.numeric(group_sizes[group_sizes > 0]))
gaps <- gaps[gaps < ncol(separated_mat)]

png(
  file.path(out_dir, "pgpg_normalized_heatmap_by_mag_separated.png"),
  width = 4200,
  height = 2400,
  res = 300
)
pheatmap(
  separated_mat,
  color = heat_colors,
  breaks = seq(-2.5, 2.5, length.out = 102),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = gaps,
  annotation_col = separated_annotation,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  fontsize_row = 8,
  border_color = NA,
  main = "Normalized plant growth-promoting functional profiles by compartment",
  legend_breaks = c(-2, 0, 2),
  legend_labels = c("Lower", "Mean", "Higher")
)
dev.off()

pdf(
  file.path(out_dir, "pgpg_normalized_heatmap_by_mag_separated.pdf"),
  width = 14,
  height = 8
)
pheatmap(
  separated_mat,
  color = heat_colors,
  breaks = seq(-2.5, 2.5, length.out = 102),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = gaps,
  annotation_col = separated_annotation,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  fontsize_row = 8,
  border_color = NA,
  main = "Normalized plant growth-promoting functional profiles by compartment",
  legend_breaks = c(-2, 0, 2),
  legend_labels = c("Lower", "Mean", "Higher")
)
dev.off()

compartments <- c("Leaves", "Roots")
mean_mat <- sapply(compartments, function(compartment) {
  rowMeans(mat[, group == compartment, drop = FALSE])
})

mean_scaled <- t(apply(mean_mat, 1, function(x) {
  if (all(x == 0)) {
    rep(0, length(x))
  } else {
    x / max(x)
  }
}))
rownames(mean_scaled) <- rownames(mean_mat)
colnames(mean_scaled) <- colnames(mean_mat)

png(
  file.path(out_dir, "pgpg_normalized_mean_by_compartment.png"),
  width = 2400,
  height = 2400,
  res = 300
)
pheatmap(
  mean_scaled,
  color = colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B"))(101),
  breaks = seq(0, 1, length.out = 102),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 8,
  fontsize_col = 10,
  border_color = "white",
  main = "Mean normalized PGP functional profiles",
  legend_breaks = c(0, 0.5, 1),
  legend_labels = c("Lower", "Intermediate", "Higher")
)
dev.off()

pdf(
  file.path(out_dir, "pgpg_normalized_mean_by_compartment.pdf"),
  width = 8,
  height = 8
)
pheatmap(
  mean_scaled,
  color = colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B"))(101),
  breaks = seq(0, 1, length.out = 102),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 8,
  fontsize_col = 10,
  border_color = "white",
  main = "Mean normalized PGP functional profiles",
  legend_breaks = c(0, 0.5, 1),
  legend_labels = c("Lower", "Intermediate", "Higher")
)
dev.off()

write.csv(
  data.frame(Function = rownames(mean_mat), mean_mat, row.names = NULL),
  file.path(out_dir, "pgpg_normalized_means_by_compartment.csv"),
  row.names = FALSE
)

message("Heatmaps written to: ", normalizePath(out_dir))


==========================================================================================
# SCRIPT: PGPg_finder/plot_pgpg_significant_functions.R
==========================================================================================

library(ggplot2)

input_file <- "PGPg_finder/tables/summary/normalized_summary_table.txt"
out_dir <- "PGPg_finder/figures/codex_heatmaps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(
  input_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

mat <- as.matrix(dat)
storage.mode(mat) <- "numeric"
mat[is.na(mat)] <- 0

leaves_cols <- grepl("^Leaves", colnames(mat))
roots_cols <- grepl("^Roots", colnames(mat))

if (!any(leaves_cols) || !any(roots_cols)) {
  stop("Leaves or Roots MAGs were not identified from column names.")
}

positive_values <- mat[mat > 0]
pseudocount <- if (length(positive_values) > 0) min(positive_values) / 2 else 1e-8

results <- lapply(seq_len(nrow(mat)), function(i) {
  leaves <- as.numeric(mat[i, leaves_cols])
  roots <- as.numeric(mat[i, roots_cols])
  test <- suppressWarnings(wilcox.test(leaves, roots, exact = FALSE))

  mean_leaves <- mean(leaves)
  mean_roots <- mean(roots)
  log2_fc <- log2((mean_roots + pseudocount) / (mean_leaves + pseudocount))

  data.frame(
    Function = rownames(mat)[i],
    n_Leaves = length(leaves),
    n_Roots = length(roots),
    mean_Leaves = mean_leaves,
    mean_Roots = mean_roots,
    median_Leaves = median(leaves),
    median_Roots = median(roots),
    log2FC_Roots_vs_Leaves = log2_fc,
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
})

results <- do.call(rbind, results)
results$FDR <- p.adjust(results$p_value, method = "BH")
results$Enriched_in <- ifelse(
  results$log2FC_Roots_vs_Leaves > 0,
  "Roots",
  ifelse(results$log2FC_Roots_vs_Leaves < 0, "Leaves", "No difference")
)
results$Significant <- results$FDR < 0.05
results <- results[order(results$FDR, -abs(results$log2FC_Roots_vs_Leaves)), ]

write.csv(
  results,
  file.path(out_dir, "pgpg_leaves_vs_roots_statistics.csv"),
  row.names = FALSE
)

sig <- results[results$Significant, , drop = FALSE]

if (nrow(sig) == 0) {
  stop("No functions were significant after FDR correction.")
}

sig$Function_label <- gsub("_", " ", sig$Function, fixed = TRUE)
sig$Function_label <- factor(
  sig$Function_label,
  levels = sig$Function_label[order(sig$log2FC_Roots_vs_Leaves)]
)

colors <- c(
  Leaves = "#549B40",
  Roots = "#613F8F"
)

p <- ggplot(
  sig,
  aes(
    x = log2FC_Roots_vs_Leaves,
    y = Function_label,
    color = Enriched_in,
    size = -log10(FDR)
  )
) +
  geom_vline(xintercept = 0, color = "gray45", linewidth = 0.5) +
  geom_segment(
    aes(x = 0, xend = log2FC_Roots_vs_Leaves, yend = Function_label),
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  geom_point(alpha = 0.95) +
  scale_color_manual(values = colors) +
  scale_size_continuous(range = c(3, 8), name = expression(-log[10]("FDR"))) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.title.y = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Significantly different PGP functions between Leaves and Roots",
    subtitle = "Wilcoxon rank-sum tests with Benjamini-Hochberg correction (FDR < 0.05)",
    x = expression(log[2](" fold change (Roots / Leaves)")),
    color = "Enriched in"
  )

ggsave(
  file.path(out_dir, "pgpg_significant_functions_leaves_vs_roots.png"),
  p,
  width = 10,
  height = max(5, 0.35 * nrow(sig) + 2.2),
  dpi = 300
)

ggsave(
  file.path(out_dir, "pgpg_significant_functions_leaves_vs_roots.pdf"),
  p,
  width = 10,
  height = max(5, 0.35 * nrow(sig) + 2.2)
)

message(
  "Significant functions: ", nrow(sig),
  " of ", nrow(results),
  ". Results written to: ", normalizePath(out_dir)
)


==========================================================================================
# SCRIPT: PGPg_finder/plot_pgpg_functional_tree.R
==========================================================================================

library(ape)

input_file <- "PGPg_finder/tables/summary/normalized_summary_table.txt"
out_dir <- "PGPg_finder/figures/codex_heatmaps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(
  input_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

mat <- as.matrix(dat)
storage.mode(mat) <- "numeric"
mat[is.na(mat)] <- 0

# Standardize each function so that high-abundance functions do not dominate
# the clustering solely because of their numerical scale.
functional_z <- t(apply(mat, 1, function(x) {
  sd_x <- sd(x)
  if (is.na(sd_x) || sd_x == 0) {
    rep(0, length(x))
  } else {
    (x - mean(x)) / sd_x
  }
}))
rownames(functional_z) <- rownames(mat)
colnames(functional_z) <- colnames(mat)

mag_profiles <- t(functional_z)
functional_distance <- dist(mag_profiles, method = "euclidean")
functional_cluster <- hclust(functional_distance, method = "ward.D2")
functional_tree <- as.phylo(functional_cluster)

tip_group <- ifelse(
  grepl("^Leaves", functional_tree$tip.label),
  "Leaves",
  ifelse(grepl("^Roots", functional_tree$tip.label), "Roots", NA_character_)
)

group_colors <- c(
  Leaves = "#549B40",
  Roots = "#613F8F"
)
tip_colors <- unname(group_colors[tip_group])
tip_colors[is.na(tip_colors)] <- "#777777"

write.tree(
  functional_tree,
  file = file.path(out_dir, "pgpg_functional_similarity_tree.nwk")
)

png(
  file.path(out_dir, "pgpg_functional_similarity_tree_circular.png"),
  width = 3000,
  height = 3000,
  res = 300
)
par(mar = c(1, 1, 4, 1))
plot(
  functional_tree,
  type = "fan",
  show.tip.label = FALSE,
  edge.color = "#555555",
  edge.width = 0.7,
  no.margin = FALSE
)
tiplabels(
  pch = 21,
  bg = tip_colors,
  col = tip_colors,
  cex = 0.75
)
title(
  main = "Functional similarity among MAGs based on PGP profiles",
  sub = "Ward clustering of Euclidean distances calculated from standardized functions",
  cex.main = 1.3,
  cex.sub = 0.9
)
legend(
  "topleft",
  legend = names(group_colors),
  pch = 21,
  pt.bg = group_colors,
  col = group_colors,
  pt.cex = 1.2,
  bty = "n",
  title = "Compartment"
)
dev.off()

png(
  file.path(out_dir, "pgpg_functional_similarity_tree_circular_labeled.png"),
  width = 6000,
  height = 6000,
  res = 300
)
par(mar = c(1, 1, 5, 1))
plot(
  functional_tree,
  type = "fan",
  show.tip.label = TRUE,
  tip.color = tip_colors,
  font = 1,
  cex = 0.38,
  label.offset = 0.25,
  edge.color = "#555555",
  edge.width = 0.6,
  no.margin = FALSE
)
tiplabels(
  pch = 21,
  bg = tip_colors,
  col = tip_colors,
  cex = 0.45,
  adj = 0.5
)
title(
  main = "Functional similarity among MAGs based on PGP profiles",
  sub = "Tip labels show MAG identifiers; colors indicate sampling compartment",
  cex.main = 1.4,
  cex.sub = 1
)
legend(
  "topleft",
  legend = names(group_colors),
  pch = 21,
  pt.bg = group_colors,
  col = group_colors,
  pt.cex = 1.3,
  cex = 1.1,
  bty = "n",
  title = "Compartment"
)
dev.off()

pdf(
  file.path(out_dir, "pgpg_functional_similarity_tree_circular.pdf"),
  width = 12,
  height = 12
)
par(mar = c(1, 1, 4, 1))
plot(
  functional_tree,
  type = "fan",
  show.tip.label = TRUE,
  tip.color = tip_colors,
  cex = 0.22,
  label.offset = 0.15,
  edge.color = "#555555",
  edge.width = 0.55,
  no.margin = FALSE
)
tiplabels(
  pch = 21,
  bg = tip_colors,
  col = tip_colors,
  cex = 0.45,
  adj = 0.5
)
title(
  main = "Functional similarity among MAGs based on PGP profiles",
  sub = "Ward clustering of Euclidean distances calculated from standardized functions",
  cex.main = 1.2,
  cex.sub = 0.85
)
legend(
  "topleft",
  legend = names(group_colors),
  pch = 21,
  pt.bg = group_colors,
  col = group_colors,
  pt.cex = 1.1,
  bty = "n",
  title = "Compartment"
)
dev.off()

pdf(
  file.path(out_dir, "pgpg_functional_similarity_tree_rectangular_labeled.pdf"),
  width = 16,
  height = 42
)
par(mar = c(4, 1, 4, 14))
plot(
  functional_tree,
  type = "phylogram",
  direction = "rightwards",
  show.tip.label = TRUE,
  tip.color = tip_colors,
  cex = 0.45,
  label.offset = 0.1,
  edge.color = "#555555",
  edge.width = 0.6,
  no.margin = FALSE
)
title(
  main = "Functional similarity among MAGs based on PGP profiles",
  sub = "Tip labels show MAG identifiers; colors indicate sampling compartment"
)
legend(
  "bottomleft",
  legend = names(group_colors),
  pch = 21,
  pt.bg = group_colors,
  col = group_colors,
  bty = "n",
  title = "Compartment"
)
dev.off()

png(
  file.path(out_dir, "pgpg_functional_similarity_dendrogram.png"),
  width = 4200,
  height = 1800,
  res = 300
)
par(mar = c(3, 4, 4, 1))
plot(
  functional_cluster,
  labels = FALSE,
  hang = -1,
  main = "Functional clustering of MAGs based on PGP profiles",
  xlab = "MAGs",
  ylab = "Functional dissimilarity",
  sub = ""
)
ordered_groups <- ifelse(
  grepl("^Leaves", functional_cluster$labels[functional_cluster$order]),
  "Leaves",
  "Roots"
)
usr <- par("usr")
points(
  seq_along(ordered_groups),
  rep(usr[3], length(ordered_groups)),
  pch = 15,
  col = group_colors[ordered_groups],
  cex = 0.8,
  xpd = NA
)
legend(
  "topright",
  legend = names(group_colors),
  pch = 15,
  col = group_colors,
  bty = "n",
  title = "Compartment"
)
dev.off()

cluster_table <- data.frame(
  MAG = functional_cluster$labels,
  Compartment = ifelse(
    grepl("^Leaves", functional_cluster$labels),
    "Leaves",
    "Roots"
  ),
  Leaf_order = match(functional_cluster$labels, functional_cluster$labels[functional_cluster$order]),
  stringsAsFactors = FALSE
)
write.csv(
  cluster_table,
  file.path(out_dir, "pgpg_functional_tree_tip_order.csv"),
  row.names = FALSE
)

message("Functional tree written to: ", normalizePath(out_dir))


==========================================================================================
# SCRIPT: PGPg_finder/plot_pgpg_heatmap_with_taxonomy.R
==========================================================================================

library(pheatmap)

input_file <- "PGPg_finder/tables/summary/normalized_summary_table.txt"
taxonomy_file <- "PGPg_finder/mags_table.csv"
out_dir <- "PGPg_finder/figures/codex_heatmaps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(
  input_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

mat <- as.matrix(dat)
storage.mode(mat) <- "numeric"
mat[is.na(mat)] <- 0

taxonomy <- read.csv(
  taxonomy_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

taxonomy <- taxonomy[match(colnames(mat), taxonomy$genome), ]
if (any(is.na(taxonomy$genome))) {
  stop("Some heatmap MAG identifiers were not found in mags_table.csv.")
}

row_zscore <- t(apply(mat, 1, function(x) {
  sd_x <- sd(x)
  if (is.na(sd_x) || sd_x == 0) {
    rep(0, length(x))
  } else {
    (x - mean(x)) / sd_x
  }
}))
rownames(row_zscore) <- rownames(mat)
colnames(row_zscore) <- colnames(mat)
row_zscore[row_zscore > 2.5] <- 2.5
row_zscore[row_zscore < -2.5] <- -2.5

cluster_within_group <- function(matrix_data, selected) {
  if (length(selected) <= 1) {
    return(selected)
  }
  local_tree <- hclust(
    dist(t(matrix_data[, selected, drop = FALSE])),
    method = "complete"
  )
  selected[local_tree$order]
}

leaves_indices <- which(taxonomy$Compartment == "Leaves")
roots_indices <- which(taxonomy$Compartment == "Roots")
leaves_order <- cluster_within_group(row_zscore, leaves_indices)
roots_order <- cluster_within_group(row_zscore, roots_indices)
column_order <- c(leaves_order, roots_order)

ordered_mat <- row_zscore[, column_order, drop = FALSE]
ordered_taxonomy <- taxonomy[column_order, ]

annotation_col <- data.frame(
  Compartment = factor(ordered_taxonomy$Compartment, levels = c("Leaves", "Roots")),
  Domain = factor(ordered_taxonomy$Domain),
  Phylum = factor(ordered_taxonomy$Phylum),
  row.names = ordered_taxonomy$genome,
  check.names = FALSE
)

domain_levels <- levels(annotation_col$Domain)
phylum_levels <- levels(annotation_col$Phylum)

domain_base <- c(
  Bacteria = "#3B82A0",
  Archaea = "#D08B3E"
)
domain_colors <- setNames(
  domain_base[domain_levels],
  domain_levels
)

phylum_colors <- setNames(
  hcl.colors(length(phylum_levels), palette = "Dynamic"),
  phylum_levels
)

annotation_colors <- list(
  Compartment = c(
    Leaves = "#549B40",
    Roots = "#613F8F"
  ),
  Domain = domain_colors,
  Phylum = phylum_colors
)

heatmap_colors <- colorRampPalette(
  c("#2166AC", "#F7F7F7", "#B2182B")
)(101)

gap_position <- length(leaves_order)

draw_taxonomy_heatmap <- function() {
  pheatmap(
    ordered_mat,
    color = heatmap_colors,
    breaks = seq(-2.5, 2.5, length.out = 102),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    gaps_col = gap_position,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    annotation_legend = TRUE,
    show_colnames = FALSE,
    fontsize_row = 8,
    border_color = NA,
    main = "Normalized plant growth-promoting functional profiles with MAG taxonomy",
    legend_breaks = c(-2, 0, 2),
    legend_labels = c("Lower", "Mean", "Higher")
  )
}

png(
  file.path(out_dir, "pgpg_normalized_heatmap_by_mag_with_taxonomy.png"),
  width = 5400,
  height = 3000,
  res = 300
)
draw_taxonomy_heatmap()
dev.off()

pdf(
  file.path(out_dir, "pgpg_normalized_heatmap_by_mag_with_taxonomy.pdf"),
  width = 18,
  height = 10
)
draw_taxonomy_heatmap()
dev.off()

taxonomy_key <- ordered_taxonomy[, c(
  "genome",
  "Compartment",
  "Domain",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species",
  "Completeness",
  "Contamination",
  "mag_size"
)]
taxonomy_key$Heatmap_order <- seq_len(nrow(taxonomy_key))
taxonomy_key <- taxonomy_key[, c(
  "Heatmap_order",
  setdiff(names(taxonomy_key), "Heatmap_order")
)]

write.csv(
  taxonomy_key,
  file.path(out_dir, "pgpg_heatmap_MAG_taxonomy_order.csv"),
  row.names = FALSE,
  na = "NA"
)

write.csv(
  data.frame(
    Phylum = names(phylum_colors),
    Color = unname(phylum_colors),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "pgpg_heatmap_phylum_color_key.csv"),
  row.names = FALSE
)

message(
  "Taxonomy-annotated heatmap written to: ",
  normalizePath(out_dir),
  ". MAGs: ", ncol(ordered_mat),
  "; phyla: ", length(phylum_levels)
)


==========================================================================================
# SCRIPT: PGPg_finder/plot_pgpg_heatmap_family_genus.R
==========================================================================================

library(pheatmap)

input_file <- "PGPg_finder/tables/summary/normalized_summary_table.txt"
taxonomy_file <- "PGPg_finder/mags_table.csv"
out_dir <- "PGPg_finder/figures/codex_heatmaps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(
  input_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

mat <- as.matrix(dat)
storage.mode(mat) <- "numeric"
mat[is.na(mat)] <- 0

taxonomy <- read.csv(
  taxonomy_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
taxonomy <- taxonomy[match(colnames(mat), taxonomy$genome), ]

if (any(is.na(taxonomy$genome))) {
  stop("Some heatmap MAG identifiers were not found in mags_table.csv.")
}

row_zscore <- t(apply(mat, 1, function(x) {
  sd_x <- sd(x)
  if (is.na(sd_x) || sd_x == 0) {
    rep(0, length(x))
  } else {
    (x - mean(x)) / sd_x
  }
}))
rownames(row_zscore) <- rownames(mat)
colnames(row_zscore) <- colnames(mat)
row_zscore[row_zscore > 2.5] <- 2.5
row_zscore[row_zscore < -2.5] <- -2.5

cluster_within_group <- function(matrix_data, selected) {
  if (length(selected) <= 1) {
    return(selected)
  }
  tree <- hclust(
    dist(t(matrix_data[, selected, drop = FALSE])),
    method = "complete"
  )
  selected[tree$order]
}

leaves_indices <- which(taxonomy$Compartment == "Leaves")
roots_indices <- which(taxonomy$Compartment == "Roots")
column_order <- c(
  cluster_within_group(row_zscore, leaves_indices),
  cluster_within_group(row_zscore, roots_indices)
)

ordered_mat <- row_zscore[, column_order, drop = FALSE]
ordered_taxonomy <- taxonomy[column_order, ]

clean_taxon <- function(x, fallback) {
  x[is.na(x) | trimws(x) == ""] <- fallback
  x
}

family <- clean_taxon(ordered_taxonomy$Family, "Unclassified family")
genus <- clean_taxon(ordered_taxonomy$Genus, "Unclassified genus")
taxonomy_labels <- paste(family, genus, sep = " | ")

annotation_col <- data.frame(
  Compartment = factor(
    ordered_taxonomy$Compartment,
    levels = c("Leaves", "Roots")
  ),
  Domain = factor(ordered_taxonomy$Domain),
  row.names = ordered_taxonomy$genome,
  check.names = FALSE
)

domain_levels <- levels(annotation_col$Domain)
domain_base <- c(
  Bacteria = "#3B82A0",
  Archaea = "#D08B3E"
)

annotation_colors <- list(
  Compartment = c(
    Leaves = "#549B40",
    Roots = "#613F8F"
  ),
  Domain = setNames(domain_base[domain_levels], domain_levels)
)

heatmap_colors <- colorRampPalette(
  c("#2166AC", "#F7F7F7", "#B2182B")
)(101)

gap_position <- length(leaves_indices)

draw_heatmap <- function(fontsize_col) {
  pheatmap(
    ordered_mat,
    color = heatmap_colors,
    breaks = seq(-2.5, 2.5, length.out = 102),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    gaps_col = gap_position,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    annotation_legend = TRUE,
    show_colnames = TRUE,
    labels_col = taxonomy_labels,
    angle_col = 90,
    fontsize_col = fontsize_col,
    fontsize_row = 8,
    border_color = NA,
    main = "Normalized PGP functional profiles annotated by MAG family and genus",
    legend_breaks = c(-2, 0, 2),
    legend_labels = c("Lower", "Mean", "Higher")
  )
}

png(
  file.path(out_dir, "pgpg_normalized_heatmap_by_mag_family_genus.png"),
  width = 12000,
  height = 4800,
  res = 300
)
draw_heatmap(fontsize_col = 5.2)
dev.off()

pdf(
  file.path(out_dir, "pgpg_normalized_heatmap_by_mag_family_genus.pdf"),
  width = 40,
  height = 16,
  useDingbats = FALSE
)
draw_heatmap(fontsize_col = 5)
dev.off()

taxonomy_key <- data.frame(
  Heatmap_order = seq_len(nrow(ordered_taxonomy)),
  MAG = ordered_taxonomy$genome,
  Compartment = ordered_taxonomy$Compartment,
  Domain = ordered_taxonomy$Domain,
  Phylum = ordered_taxonomy$Phylum,
  Class = ordered_taxonomy$Class,
  Order = ordered_taxonomy$Order,
  Family = family,
  Genus = genus,
  Heatmap_label = taxonomy_labels,
  stringsAsFactors = FALSE
)

write.csv(
  taxonomy_key,
  file.path(out_dir, "pgpg_heatmap_family_genus_labels.csv"),
  row.names = FALSE,
  na = "NA"
)

message(
  "Family/genus heatmap written to: ",
  normalizePath(out_dir),
  ". MAGs: ", ncol(ordered_mat)
)


==========================================================================================
# SCRIPT: PGPg_finder/plot_pgpg_heatmap_selected_phyla.R
==========================================================================================

library(pheatmap)

input_file <- "PGPg_finder/tables/summary/normalized_summary_table.txt"
taxonomy_file <- "PGPg_finder/mags_table.csv"
out_dir <- "PGPg_finder/figures/codex_heatmaps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

selected_phyla <- c(
  "Chloroflexota",
  "Actinomycetota",
  "Cyanobacteriota",
  "Bacillota",
  "Spirochaetota",
  "Pseudomonadota",
  "Desulfobacterota",
  "Myxococcota",
  "Thermoproteota"
)

dat <- read.delim(
  input_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
mat <- as.matrix(dat)
storage.mode(mat) <- "numeric"
mat[is.na(mat)] <- 0

taxonomy <- read.csv(
  taxonomy_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
taxonomy <- taxonomy[match(colnames(mat), taxonomy$genome), ]

keep <- taxonomy$Phylum %in% selected_phyla
mat <- mat[, keep, drop = FALSE]
taxonomy <- taxonomy[keep, , drop = FALSE]

if (ncol(mat) == 0) {
  stop("No MAGs matched the selected phyla.")
}

row_zscore <- t(apply(mat, 1, function(x) {
  sd_x <- sd(x)
  if (is.na(sd_x) || sd_x == 0) {
    rep(0, length(x))
  } else {
    (x - mean(x)) / sd_x
  }
}))
rownames(row_zscore) <- rownames(mat)
colnames(row_zscore) <- colnames(mat)
row_zscore[row_zscore > 2.5] <- 2.5
row_zscore[row_zscore < -2.5] <- -2.5

cluster_within_group <- function(matrix_data, selected) {
  if (length(selected) <= 1) {
    return(selected)
  }
  tree <- hclust(
    dist(t(matrix_data[, selected, drop = FALSE])),
    method = "complete"
  )
  selected[tree$order]
}

leaves_indices <- which(taxonomy$Compartment == "Leaves")
roots_indices <- which(taxonomy$Compartment == "Roots")
column_order <- c(
  cluster_within_group(row_zscore, leaves_indices),
  cluster_within_group(row_zscore, roots_indices)
)

ordered_mat <- row_zscore[, column_order, drop = FALSE]
ordered_taxonomy <- taxonomy[column_order, , drop = FALSE]

clean_taxon <- function(x, fallback) {
  x[is.na(x) | trimws(x) == ""] <- fallback
  x
}

family <- clean_taxon(ordered_taxonomy$Family, "Unclassified family")
genus <- clean_taxon(ordered_taxonomy$Genus, "Unclassified genus")
taxonomy_labels <- paste(family, genus, sep = " | ")

annotation_col <- data.frame(
  Compartment = factor(
    ordered_taxonomy$Compartment,
    levels = c("Leaves", "Roots")
  ),
  Domain = factor(ordered_taxonomy$Domain),
  Phylum = factor(ordered_taxonomy$Phylum, levels = selected_phyla),
  row.names = ordered_taxonomy$genome,
  check.names = FALSE
)

domain_levels <- levels(annotation_col$Domain)
domain_base <- c(
  Bacteria = "#3B82A0",
  Archaea = "#D08B3E"
)

phylum_colors <- setNames(
  hcl.colors(length(selected_phyla), palette = "Dynamic"),
  selected_phyla
)

annotation_colors <- list(
  Compartment = c(
    Leaves = "#549B40",
    Roots = "#613F8F"
  ),
  Domain = setNames(domain_base[domain_levels], domain_levels),
  Phylum = phylum_colors
)

heatmap_colors <- colorRampPalette(
  c("#2166AC", "#F7F7F7", "#B2182B")
)(101)

gap_position <- length(leaves_indices)

draw_heatmap <- function(fontsize_col) {
  pheatmap(
    ordered_mat,
    color = heatmap_colors,
    breaks = seq(-2.5, 2.5, length.out = 102),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    gaps_col = gap_position,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    annotation_legend = TRUE,
    show_colnames = TRUE,
    labels_col = taxonomy_labels,
    angle_col = 90,
    fontsize_col = fontsize_col,
    fontsize_row = 8,
    border_color = NA,
    main = "PGP functional profiles of MAGs from selected phyla",
    legend_breaks = c(-2, 0, 2),
    legend_labels = c("Lower", "Mean", "Higher")
  )
}

png(
  file.path(out_dir, "pgpg_heatmap_selected_phyla_family_genus.png"),
  width = 10500,
  height = 4700,
  res = 300
)
draw_heatmap(fontsize_col = 5.5)
dev.off()

pdf(
  file.path(out_dir, "pgpg_heatmap_selected_phyla_family_genus.pdf"),
  width = 35,
  height = 15.5,
  useDingbats = FALSE
)
draw_heatmap(fontsize_col = 5.2)
dev.off()

taxonomy_key <- data.frame(
  Heatmap_order = seq_len(nrow(ordered_taxonomy)),
  MAG = ordered_taxonomy$genome,
  Compartment = ordered_taxonomy$Compartment,
  Domain = ordered_taxonomy$Domain,
  Phylum = ordered_taxonomy$Phylum,
  Class = ordered_taxonomy$Class,
  Order = ordered_taxonomy$Order,
  Family = family,
  Genus = genus,
  Heatmap_label = taxonomy_labels,
  stringsAsFactors = FALSE
)

write.csv(
  taxonomy_key,
  file.path(out_dir, "pgpg_heatmap_selected_phyla_taxonomy.csv"),
  row.names = FALSE,
  na = "NA"
)

write.csv(
  data.frame(
    Phylum = names(phylum_colors),
    Color = unname(phylum_colors),
    MAG_count = as.integer(table(factor(
      ordered_taxonomy$Phylum,
      levels = selected_phyla
    ))),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "pgpg_heatmap_selected_phyla_color_key.csv"),
  row.names = FALSE
)

message(
  "Selected-phyla heatmap written to: ",
  normalizePath(out_dir),
  ". MAGs: ", ncol(ordered_mat)
)


==========================================================================================
# SCRIPT: PGPg_finder/plot_pgpg_heatmap_phylum_groups.R
==========================================================================================

library(pheatmap)

input_file <- "PGPg_finder/tables/summary/normalized_summary_table.txt"
taxonomy_file <- "PGPg_finder/mags_table.csv"
out_dir <- "PGPg_finder/figures/codex_heatmaps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

selected_phyla <- c(
  "Chloroflexota",
  "Actinomycetota",
  "Cyanobacteria",
  "Bacillota",
  "Spirochaetota",
  "Pseudomonadota",
  "Desulfobacterota",
  "Myxococcota",
  "Thermoproteota"
)

dat <- read.delim(
  input_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
mat <- as.matrix(dat)
storage.mode(mat) <- "numeric"
mat[is.na(mat)] <- 0

taxonomy <- read.csv(
  taxonomy_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
taxonomy <- taxonomy[match(colnames(mat), taxonomy$genome), ]

if (any(is.na(taxonomy$genome))) {
  stop("Some heatmap MAG identifiers were not found in mags_table.csv.")
}

taxonomy$Phylum_group <- taxonomy$Phylum
taxonomy$Phylum_group[taxonomy$Phylum_group == "Cyanobacteriota"] <- "Cyanobacteria"
taxonomy$Phylum_group[
  is.na(taxonomy$Phylum_group) |
    !taxonomy$Phylum_group %in% selected_phyla
] <- "Others"

row_zscore <- t(apply(mat, 1, function(x) {
  sd_x <- sd(x)
  if (is.na(sd_x) || sd_x == 0) {
    rep(0, length(x))
  } else {
    (x - mean(x)) / sd_x
  }
}))
rownames(row_zscore) <- rownames(mat)
colnames(row_zscore) <- colnames(mat)
row_zscore[row_zscore > 2.5] <- 2.5
row_zscore[row_zscore < -2.5] <- -2.5

cluster_within_group <- function(matrix_data, selected) {
  if (length(selected) <= 1) {
    return(selected)
  }
  tree <- hclust(
    dist(t(matrix_data[, selected, drop = FALSE])),
    method = "complete"
  )
  selected[tree$order]
}

leaves_indices <- which(taxonomy$Compartment == "Leaves")
roots_indices <- which(taxonomy$Compartment == "Roots")
column_order <- c(
  cluster_within_group(row_zscore, leaves_indices),
  cluster_within_group(row_zscore, roots_indices)
)

ordered_mat <- row_zscore[, column_order, drop = FALSE]
ordered_taxonomy <- taxonomy[column_order, , drop = FALSE]

phylum_levels <- c(selected_phyla, "Others")
annotation_col <- data.frame(
  Compartment = factor(
    ordered_taxonomy$Compartment,
    levels = c("Leaves", "Roots")
  ),
  Phylum = factor(
    ordered_taxonomy$Phylum_group,
    levels = phylum_levels
  ),
  row.names = ordered_taxonomy$genome,
  check.names = FALSE
)

phylum_colors <- c(
  Chloroflexota = "#D98C7C",
  Actinomycetota = "#C9A66B",
  Cyanobacteria = "#8DB255",
  Bacillota = "#55B88A",
  Spirochaetota = "#30B8B0",
  Pseudomonadota = "#48A7C5",
  Desulfobacterota = "#7A8FD1",
  Myxococcota = "#A681CF",
  Thermoproteota = "#D07EB5",
  Others = "#B8B8B8"
)

annotation_colors <- list(
  Compartment = c(
    Leaves = "#549B40",
    Roots = "#613F8F"
  ),
  Phylum = phylum_colors
)

heatmap_colors <- colorRampPalette(
  c("#2166AC", "#F7F7F7", "#B2182B")
)(101)

gap_position <- length(leaves_indices)

draw_heatmap <- function() {
  pheatmap(
    ordered_mat,
    color = heatmap_colors,
    breaks = seq(-2.5, 2.5, length.out = 102),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    gaps_col = gap_position,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    annotation_legend = TRUE,
    show_colnames = FALSE,
    fontsize_row = 8,
    border_color = NA,
    main = "Normalized PGP functional profiles of all MAGs by phylum group",
    legend_breaks = c(-2, 0, 2),
    legend_labels = c("Lower", "Mean", "Higher")
  )
}

png(
  file.path(out_dir, "pgpg_heatmap_all_MAGs_selected_phyla_and_others.png"),
  width = 4800,
  height = 2700,
  res = 300
)
draw_heatmap()
dev.off()

pdf(
  file.path(out_dir, "pgpg_heatmap_all_MAGs_selected_phyla_and_others.pdf"),
  width = 16,
  height = 9,
  useDingbats = FALSE
)
draw_heatmap()
dev.off()

taxonomy_key <- data.frame(
  Heatmap_order = seq_len(nrow(ordered_taxonomy)),
  MAG = ordered_taxonomy$genome,
  Compartment = ordered_taxonomy$Compartment,
  Original_phylum = ordered_taxonomy$Phylum,
  Displayed_phylum = ordered_taxonomy$Phylum_group,
  stringsAsFactors = FALSE
)

write.csv(
  taxonomy_key,
  file.path(out_dir, "pgpg_heatmap_all_MAGs_phylum_groups.csv"),
  row.names = FALSE,
  na = "NA"
)

write.csv(
  data.frame(
    Phylum_group = phylum_levels,
    Color = unname(phylum_colors[phylum_levels]),
    MAG_count = as.integer(table(factor(
      ordered_taxonomy$Phylum_group,
      levels = phylum_levels
    ))),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "pgpg_heatmap_all_MAGs_phylum_group_key.csv"),
  row.names = FALSE
)

message(
  "All-MAG phylum-group heatmap written to: ",
  normalizePath(out_dir),
  ". MAGs: ", ncol(ordered_mat),
  "; Others: ", sum(ordered_taxonomy$Phylum_group == "Others")
)


==========================================================================================
# SCRIPT: PGPg_finder/plot_pgpg_heatmap_filtered_root_taxa.R
==========================================================================================

library(grid)
library(gtable)
library(pheatmap)

input_file <- "PGPg_finder/tables/summary/normalized_summary_table.txt"
taxonomy_file <- "PGPg_finder/mags_table.csv"
out_dir <- Sys.getenv("PGPG_OUT_DIR", "PGPg_finder/figures/codex_heatmaps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

removed_mags_env <- Sys.getenv("PGPG_REMOVE_MAGS", "")
removed_mags <- if (nzchar(removed_mags_env)) {
  trimws(strsplit(removed_mags_env, ",", fixed = TRUE)[[1]])
} else {
  character(0)
}

selected_phyla <- c(
  "Chloroflexota",
  "Actinomycetota",
  "Cyanobacteria",
  "Bacillota",
  "Spirochaetota",
  "Pseudomonadota",
  "Desulfobacterota",
  "Myxococcota",
  "Thermoproteota"
)

functions_to_remove <- c(
  "NITRIFICATION",
  "MOTILITY CHEMOTAXIS",
  "INSECTICIDAL COMPOUNDS",
  "NEMATICIDAL COMPOUNDS",
  "HEAVY METAL DETOXIFICATION",
  "ORGANIC VOLATILES",
  "P-SOLUBILISATION GLUCONIC ACID-PQQ",
  "SURFACE ATTACHMENT",
  "SECRETION SYSTEMS"
)

dat <- read.delim(
  input_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(dat) <- trimws(rownames(dat))
dat <- dat[!rownames(dat) %in% functions_to_remove, , drop = FALSE]

mat <- as.matrix(dat)
storage.mode(mat) <- "numeric"
mat[is.na(mat)] <- 0
if (length(removed_mags) > 0) {
  mat <- mat[, !colnames(mat) %in% removed_mags, drop = FALSE]
}

taxonomy <- read.csv(
  taxonomy_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
taxonomy <- taxonomy[match(colnames(mat), taxonomy$genome), ]

if (any(is.na(taxonomy$genome))) {
  stop("Some heatmap MAG identifiers were not found in mags_table.csv.")
}

taxonomy$Phylum_group <- taxonomy$Phylum
taxonomy$Phylum_group[taxonomy$Phylum_group == "Cyanobacteriota"] <- "Cyanobacteria"
taxonomy$Phylum_group[
  is.na(taxonomy$Phylum_group) |
    !taxonomy$Phylum_group %in% selected_phyla
] <- "Others"

taxonomy$Root_taxon <- "Other MAGs"
taxonomy$Root_taxon[
  taxonomy$Compartment == "Roots" &
    taxonomy$Family == "Sedimenticolaceae"
] <- "Sedimenticolaceae"
taxonomy$Root_taxon[
  taxonomy$Compartment == "Roots" &
    taxonomy$Order == "Desulfobacterales"
] <- "Desulfobacterales"

row_zscore <- t(apply(mat, 1, function(x) {
  sd_x <- sd(x)
  if (is.na(sd_x) || sd_x == 0) {
    rep(0, length(x))
  } else {
    (x - mean(x)) / sd_x
  }
}))
rownames(row_zscore) <- rownames(mat)
colnames(row_zscore) <- colnames(mat)
row_zscore[row_zscore > 2.5] <- 2.5
row_zscore[row_zscore < -2.5] <- -2.5

cluster_within_group <- function(matrix_data, selected) {
  if (length(selected) <= 1) {
    return(selected)
  }
  tree <- hclust(
    dist(t(matrix_data[, selected, drop = FALSE])),
    method = "complete"
  )
  selected[tree$order]
}

leaves_indices <- which(taxonomy$Compartment == "Leaves")
roots_indices <- which(taxonomy$Compartment == "Roots")
column_order <- c(
  cluster_within_group(row_zscore, leaves_indices),
  cluster_within_group(row_zscore, roots_indices)
)

ordered_mat <- row_zscore[, column_order, drop = FALSE]
ordered_taxonomy <- taxonomy[column_order, , drop = FALSE]

phylum_levels <- c(selected_phyla, "Others")
root_taxon_levels <- c("Sedimenticolaceae", "Desulfobacterales", "Other MAGs")

annotation_col <- data.frame(
  Compartment = factor(
    ordered_taxonomy$Compartment,
    levels = c("Leaves", "Roots")
  ),
  Phylum = factor(
    ordered_taxonomy$Phylum_group,
    levels = phylum_levels
  ),
  row.names = ordered_taxonomy$genome,
  check.names = FALSE
)

phylum_colors <- c(
  Chloroflexota = "#D98C7C",
  Actinomycetota = "#C9A66B",
  Cyanobacteria = "#8DB255",
  Bacillota = "#55B88A",
  Spirochaetota = "#30B8B0",
  Pseudomonadota = "#48A7C5",
  Desulfobacterota = "#7A8FD1",
  Myxococcota = "#A681CF",
  Thermoproteota = "#D07EB5",
  Others = "#B8B8B8"
)

annotation_colors <- list(
  Compartment = c(
    Leaves = "#549B40",
    Roots = "#613F8F"
  ),
  Phylum = phylum_colors
)

root_taxon_colors <- c(
  Sedimenticolaceae = "#E76F00",
  Desulfobacterales = "#005F73",
  `Other MAGs` = "#D8D8D8"
)

heatmap_colors <- colorRampPalette(
  c("#2166AC", "#F7F7F7", "#B2182B")
)(101)

gap_position <- length(leaves_indices)

add_png_padding <- function(input_path, output_path = input_path, pad_px = 260) {
  image <- png::readPNG(input_path)
  dims <- dim(image)
  if (length(dims) == 2) {
    padded <- matrix(1, nrow = dims[1] + 2 * pad_px, ncol = dims[2] + 2 * pad_px)
    padded[
      (pad_px + 1):(pad_px + dims[1]),
      (pad_px + 1):(pad_px + dims[2])
    ] <- image
  } else {
    padded <- array(
      1,
      dim = c(dims[1] + 2 * pad_px, dims[2] + 2 * pad_px, dims[3])
    )
    padded[
      (pad_px + 1):(pad_px + dims[1]),
      (pad_px + 1):(pad_px + dims[2]),
      seq_len(dims[3])
    ] <- image
  }
  png::writePNG(padded, target = output_path)
}

make_root_taxon_strip <- function(root_taxon_vector) {
  n_cols <- length(root_taxon_vector)
  rects <- lapply(seq_len(n_cols), function(i) {
    rectGrob(
      x = (i - 0.5) / n_cols,
      y = 0.5,
      width = 1 / n_cols,
      height = 1,
      gp = gpar(
        fill = unname(root_taxon_colors[root_taxon_vector[i]]),
        col = NA
      )
    )
  })
  do.call(grobTree, rects)
}

make_root_taxon_legend <- function() {
  labels <- names(root_taxon_colors)
  grobs <- list(
    textGrob(
      "Root taxon",
      x = unit(0, "npc"),
      y = unit(1, "npc"),
      just = c("left", "top"),
      gp = gpar(fontface = "bold", fontsize = 11)
    )
  )
  for (i in seq_along(labels)) {
    y_pos <- unit(1, "npc") - unit(18 + (i - 1) * 16, "pt")
    grobs[[length(grobs) + 1]] <- rectGrob(
      x = unit(0, "npc"),
      y = y_pos,
      width = unit(10, "pt"),
      height = unit(10, "pt"),
      just = c("left", "center"),
      gp = gpar(fill = unname(root_taxon_colors[labels[i]]), col = NA)
    )
    grobs[[length(grobs) + 1]] <- textGrob(
      labels[i],
      x = unit(14, "pt"),
      y = y_pos,
      just = c("left", "center"),
      gp = gpar(fontsize = 10)
    )
  }
  grobTree(children = do.call(gList, grobs))
}

add_bottom_root_taxon <- function(plot_gtable, root_taxon_vector) {
  matrix_layout <- plot_gtable$layout[plot_gtable$layout$name == "matrix", ]
  row_names_layout <- plot_gtable$layout[plot_gtable$layout$name == "row_names", ]
  legend_layout <- plot_gtable$layout[plot_gtable$layout$name == "annotation_legend", ]

  plot_gtable <- gtable_add_rows(
    plot_gtable,
    heights = unit(16, "pt"),
    pos = matrix_layout$b
  )
  bottom_row <- matrix_layout$b + 1

  plot_gtable <- gtable_add_grob(
    plot_gtable,
    grobs = make_root_taxon_strip(root_taxon_vector),
    t = bottom_row,
    l = matrix_layout$l,
    b = bottom_row,
    r = matrix_layout$r,
    name = "bottom_root_taxon"
  )
  plot_gtable <- gtable_add_grob(
    plot_gtable,
    grobs = textGrob(
      "Root taxon",
      x = unit(0, "npc"),
      y = unit(0.5, "npc"),
      just = c("left", "center"),
      gp = gpar(fontface = "bold", fontsize = 10)
    ),
    t = bottom_row,
    l = row_names_layout$l,
    b = bottom_row,
    r = row_names_layout$r,
    name = "bottom_root_taxon_name"
  )

  if (nrow(legend_layout) == 1) {
    plot_gtable <- gtable_add_cols(
      plot_gtable,
      widths = unit(240, "pt"),
      pos = legend_layout$r
    )
    root_legend_col <- legend_layout$r + 1
    plot_gtable <- gtable_add_grob(
      plot_gtable,
      grobs = make_root_taxon_legend(),
      t = legend_layout$t,
      l = root_legend_col,
      b = matrix_layout$t,
      r = root_legend_col,
      name = "root_taxon_legend"
    )
  }
  plot_gtable <- gtable_add_rows(
    plot_gtable,
    heights = unit(10, "pt"),
    pos = bottom_row
  )
  plot_gtable <- gtable_add_rows(
    plot_gtable,
    heights = unit(18, "pt"),
    pos = 0
  )
  plot_gtable
}

draw_heatmap <- function() {
  heatmap <- pheatmap(
    ordered_mat,
    color = heatmap_colors,
    breaks = seq(-2.5, 2.5, length.out = 102),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    gaps_col = gap_position,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    annotation_legend = TRUE,
    show_colnames = FALSE,
    fontsize_row = 9,
    border_color = NA,
    main = "Filtered PGP functional profiles of all MAGs with root taxon highlights",
    legend_breaks = c(-2, 0, 2),
    legend_labels = c("Lower", "Mean", "Higher"),
    silent = TRUE
  )
  heatmap$gtable <- add_bottom_root_taxon(
    heatmap$gtable,
    as.character(ordered_taxonomy$Root_taxon)
  )
  grid.newpage()
  pushViewport(viewport(
    x = 0.5,
    y = 0.5,
    width = unit(0.90, "npc"),
    height = unit(0.88, "npc")
  ))
  grid.draw(heatmap$gtable)
  popViewport()
}

png(
  png_path <- file.path(out_dir, "pgpg_heatmap_all_MAGs_filtered_root_taxa_bottom_highlight.png"),
  width = 7200,
  height = 3600,
  res = 300
)
draw_heatmap()
dev.off()

add_png_padding(png_path, png_path)
add_png_padding(
  png_path,
  file.path(out_dir, "pgpg_heatmap_all_MAGs_filtered_root_taxa_bottom_highlight_uncut.png"),
  pad_px = 0
)

pdf(
  file.path(out_dir, "pgpg_heatmap_all_MAGs_filtered_root_taxa_bottom_highlight.pdf"),
  width = 24,
  height = 12,
  useDingbats = FALSE
)
draw_heatmap()
dev.off()

write.csv(
  data.frame(
    Removed_function = functions_to_remove,
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "pgpg_heatmap_filtered_removed_functions.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    Heatmap_order = seq_len(nrow(ordered_taxonomy)),
    MAG = ordered_taxonomy$genome,
    Compartment = ordered_taxonomy$Compartment,
    Original_phylum = ordered_taxonomy$Phylum,
    Displayed_phylum = ordered_taxonomy$Phylum_group,
    Order = ordered_taxonomy$Order,
    Family = ordered_taxonomy$Family,
    Root_taxon_highlight = ordered_taxonomy$Root_taxon,
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "pgpg_heatmap_filtered_root_taxa_order.csv"),
  row.names = FALSE,
  na = "NA"
)

write.csv(
  data.frame(
    Removed_MAG = removed_mags,
    Found_in_input = removed_mags %in% colnames(dat),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "pgpg_heatmap_removed_MAGs.csv"),
  row.names = FALSE
)

message(
  "Filtered all-MAG heatmap written to: ",
  normalizePath(out_dir),
  ". Functions retained: ", nrow(ordered_mat),
  "; MAGs: ", ncol(ordered_mat),
  "; highlighted Sedimenticolaceae Roots MAGs: ",
  sum(ordered_taxonomy$Root_taxon == "Sedimenticolaceae"),
  "; highlighted Desulfobacterales Roots MAGs: ",
  sum(ordered_taxonomy$Root_taxon == "Desulfobacterales")
)


==========================================================================================
# SCRIPT: PGPg_finder/plot_pgpg_nitrogen_levels.R
==========================================================================================

library(ggplot2)
library(patchwork)

lv3_file <- "PGPg_finder/tables/normalized/normalized_gene_counts_Lv3.txt"
lv4_file <- "PGPg_finder/tables/normalized/normalized_gene_counts_Lv4.txt"
out_dir <- Sys.getenv("PGPG_OUT_DIR", "PGPg_finder/figures/codex_heatmaps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

removed_mags_env <- Sys.getenv("PGPG_REMOVE_MAGS", "")
removed_mags <- if (nzchar(removed_mags_env)) {
  trimws(strsplit(removed_mags_env, ",", fixed = TRUE)[[1]])
} else {
  character(0)
}

group_colors <- c(
  Leaves = "#549B40",
  Roots = "#613F8F"
)

read_function_table <- function(path) {
  dat <- read.delim(
    path,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  mat <- as.matrix(dat)
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  mat
}

get_group <- function(mag_names) {
  ifelse(
    grepl("^Leaves", mag_names),
    "Leaves",
    ifelse(grepl("^Roots", mag_names), "Roots", NA_character_)
  )
}

lv3 <- read_function_table(lv3_file)
lv4 <- read_function_table(lv4_file)
if (length(removed_mags) > 0) {
  lv3 <- lv3[, !colnames(lv3) %in% removed_mags, drop = FALSE]
  lv4 <- lv4[, !colnames(lv4) %in% removed_mags, drop = FALSE]
}

lv3_row <- "NITROGEN_ACQUISITION"
if (!lv3_row %in% rownames(lv3)) {
  stop("NITROGEN_ACQUISITION was not found in the Lv3 table.")
}

lv3_group <- get_group(colnames(lv3))
lv3_df <- data.frame(
  MAG = colnames(lv3),
  Compartment = lv3_group,
  Normalized_abundance = as.numeric(lv3[lv3_row, ]),
  stringsAsFactors = FALSE
)
lv3_df <- lv3_df[!is.na(lv3_df$Compartment), ]
lv3_df$Compartment <- factor(lv3_df$Compartment, levels = c("Leaves", "Roots"))

lv3_test <- suppressWarnings(
  wilcox.test(
    Normalized_abundance ~ Compartment,
    data = lv3_df,
    exact = FALSE
  )
)

nitrogen_rows <- grepl("^N-AQUISITION-", rownames(lv4))
nitrogen_lv4 <- lv4[nitrogen_rows, , drop = FALSE]
lv4_group <- get_group(colnames(nitrogen_lv4))
keep_cols <- !is.na(lv4_group)
nitrogen_lv4 <- nitrogen_lv4[, keep_cols, drop = FALSE]
lv4_group <- lv4_group[keep_cols]

positive_values <- nitrogen_lv4[nitrogen_lv4 > 0]
pseudocount <- if (length(positive_values) > 0) min(positive_values) / 2 else 1e-8

lv4_stats <- lapply(seq_len(nrow(nitrogen_lv4)), function(i) {
  leaves <- as.numeric(nitrogen_lv4[i, lv4_group == "Leaves"])
  roots <- as.numeric(nitrogen_lv4[i, lv4_group == "Roots"])
  test <- suppressWarnings(wilcox.test(leaves, roots, exact = FALSE))
  mean_leaves <- mean(leaves)
  mean_roots <- mean(roots)

  data.frame(
    Function = rownames(nitrogen_lv4)[i],
    n_Leaves = length(leaves),
    n_Roots = length(roots),
    mean_Leaves = mean_leaves,
    mean_Roots = mean_roots,
    median_Leaves = median(leaves),
    median_Roots = median(roots),
    log2FC_Roots_vs_Leaves = log2(
      (mean_roots + pseudocount) /
        (mean_leaves + pseudocount)
    ),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
})

lv4_stats <- do.call(rbind, lv4_stats)
lv4_stats$FDR <- p.adjust(lv4_stats$p_value, method = "BH")
lv4_stats$Enriched_in <- ifelse(
  lv4_stats$log2FC_Roots_vs_Leaves >= 0,
  "Roots",
  "Leaves"
)
lv4_stats$Significance <- ifelse(lv4_stats$FDR < 0.05, "FDR < 0.05", "Not significant")
lv4_stats$Function_label <- sub("^N-AQUISITION-", "", lv4_stats$Function)
lv4_stats$Function_label <- gsub("[_|]", " ", lv4_stats$Function_label)
lv4_stats$Function_label <- gsub("ATMOSHPHERIC", "ATMOSPHERIC", lv4_stats$Function_label)
lv4_stats$Function_label <- gsub("CYNATE", "CYANATE", lv4_stats$Function_label)
lv4_stats$Function_label <- gsub("FORMAIDE", "FORMAMIDE", lv4_stats$Function_label)
lv4_stats$Function_label <- factor(
  lv4_stats$Function_label,
  levels = lv4_stats$Function_label[order(lv4_stats$log2FC_Roots_vs_Leaves)]
)

write.csv(
  lv4_stats,
  file.path(out_dir, "pgpg_nitrogen_lv4_leaves_vs_roots_statistics.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    Level = "Lv3",
    Function = lv3_row,
    n_Leaves = sum(lv3_df$Compartment == "Leaves"),
    n_Roots = sum(lv3_df$Compartment == "Roots"),
    mean_Leaves = mean(lv3_df$Normalized_abundance[lv3_df$Compartment == "Leaves"]),
    mean_Roots = mean(lv3_df$Normalized_abundance[lv3_df$Compartment == "Roots"]),
    p_value = lv3_test$p.value
  ),
  file.path(out_dir, "pgpg_nitrogen_lv3_leaves_vs_roots_statistics.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    Removed_MAG = removed_mags,
    Found_in_Lv3 = removed_mags %in% colnames(read_function_table(lv3_file)),
    Found_in_Lv4 = removed_mags %in% colnames(read_function_table(lv4_file)),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "pgpg_nitrogen_removed_MAGs.csv"),
  row.names = FALSE
)

p_lv3 <- ggplot(
  lv3_df,
  aes(x = Compartment, y = Normalized_abundance, fill = Compartment, color = Compartment)
) +
  geom_violin(trim = FALSE, alpha = 0.22, linewidth = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.65, linewidth = 0.5) +
  geom_jitter(width = 0.1, size = 1.1, alpha = 0.55) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Overall nitrogen acquisition (Lv3)",
    subtitle = paste0("Wilcoxon P = ", format.pval(lv3_test$p.value, digits = 2)),
    x = NULL,
    y = "Normalized functional abundance"
  )

p_lv4 <- ggplot(
  lv4_stats,
  aes(
    x = log2FC_Roots_vs_Leaves,
    y = Function_label,
    color = Enriched_in,
    alpha = Significance,
    size = -log10(FDR)
  )
) +
  geom_vline(xintercept = 0, color = "gray45", linewidth = 0.5) +
  geom_segment(
    aes(x = 0, xend = log2FC_Roots_vs_Leaves, yend = Function_label),
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  geom_point() +
  scale_color_manual(values = group_colors) +
  scale_alpha_manual(values = c("FDR < 0.05" = 1, "Not significant" = 0.38)) +
  scale_size_continuous(range = c(2.5, 7), name = expression(-log[10]("FDR"))) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Nitrogen-related subfunctions (Lv4)",
    subtitle = "Wilcoxon tests with Benjamini-Hochberg correction",
    x = expression(log[2](" fold change (Roots / Leaves)")),
    color = "Enriched in",
    alpha = "Statistical support"
  )

panel <- p_lv3 + p_lv4 +
  plot_layout(widths = c(0.75, 1.75)) +
  plot_annotation(
    title = "Nitrogen-related plant growth-promoting functions in Leaves and Roots MAGs",
    tag_levels = "A"
  ) &
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave(
  file.path(out_dir, "pgpg_nitrogen_functions_leaves_vs_roots.png"),
  panel,
  width = 14,
  height = 7.5,
  dpi = 300
)

ggsave(
  file.path(out_dir, "pgpg_nitrogen_functions_leaves_vs_roots.pdf"),
  panel,
  width = 14,
  height = 7.5
)

message(
  "Nitrogen figure written to: ", normalizePath(out_dir),
  ". Significant Lv4 functions: ", sum(lv4_stats$FDR < 0.05),
  " of ", nrow(lv4_stats)
)


==========================================================================================
# SCRIPT: PGPg_finder/plot_sedimenticolaceae_vs_desulfobacterales.R
==========================================================================================

library(ggplot2)
library(patchwork)
library(pheatmap)

input_file <- "PGPg_finder/tables/summary/normalized_summary_table.txt"
taxonomy_file <- "PGPg_finder/mags_table.csv"
out_dir <- Sys.getenv("PGPG_OUT_DIR", "PGPg_finder/figures/codex_taxon_comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

removed_mags_env <- Sys.getenv("PGPG_REMOVE_MAGS", "")
removed_mags <- if (nzchar(removed_mags_env)) {
  trimws(strsplit(removed_mags_env, ",", fixed = TRUE)[[1]])
} else {
  character(0)
}

taxon_colors <- c(
  Sedimenticolaceae = "#2A9D8F",
  Desulfobacterales = "#7B5AA6"
)
compartment_colors <- c(Leaves = "#549B40", Roots = "#613F8F")

dat <- read.delim(
  input_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
mat <- as.matrix(dat)
storage.mode(mat) <- "numeric"
mat[is.na(mat)] <- 0
if (length(removed_mags) > 0) {
  mat <- mat[, !colnames(mat) %in% removed_mags, drop = FALSE]
}

taxonomy <- read.csv(
  taxonomy_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
raw_taxonomy <- taxonomy
taxonomy <- taxonomy[match(colnames(mat), taxonomy$genome), ]
if (any(is.na(taxonomy$genome))) {
  stop("Some functional-table MAG identifiers were not found in mags_table.csv.")
}

taxonomy$Comparison_group <- NA_character_
taxonomy$Comparison_group[taxonomy$Family == "Sedimenticolaceae"] <-
  "Sedimenticolaceae"
taxonomy$Comparison_group[
  is.na(taxonomy$Comparison_group) &
    taxonomy$Order == "Desulfobacterales"
] <- "Desulfobacterales"

selected <- which(
  !is.na(taxonomy$Comparison_group) &
    taxonomy$Compartment == "Roots"
)
selected_mat <- mat[, selected, drop = FALSE]
selected_taxonomy <- taxonomy[selected, , drop = FALSE]

cluster_group <- function(matrix_data, indices) {
  if (length(indices) <= 1) {
    return(indices)
  }
  indices[hclust(dist(t(matrix_data[, indices, drop = FALSE])))$order]
}

z_mat <- t(apply(selected_mat, 1, function(x) {
  sx <- sd(x)
  if (is.na(sx) || sx == 0) rep(0, length(x)) else (x - mean(x)) / sx
}))
rownames(z_mat) <- rownames(selected_mat)
colnames(z_mat) <- colnames(selected_mat)
z_mat[z_mat > 2.5] <- 2.5
z_mat[z_mat < -2.5] <- -2.5

group_levels <- c("Sedimenticolaceae", "Desulfobacterales")
column_order <- unlist(lapply(group_levels, function(group_name) {
  idx <- which(selected_taxonomy$Comparison_group == group_name)
  cluster_group(z_mat, idx)
}), use.names = FALSE)

heat_mat <- z_mat[, column_order, drop = FALSE]
heat_taxonomy <- selected_taxonomy[column_order, , drop = FALSE]
annotation_col <- data.frame(
  Taxon = factor(heat_taxonomy$Comparison_group, levels = group_levels),
  row.names = heat_taxonomy$genome
)

heatmap_result <- pheatmap(
  heat_mat,
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
  breaks = seq(-2.5, 2.5, length.out = 102),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = sum(heat_taxonomy$Comparison_group == "Sedimenticolaceae"),
  annotation_col = annotation_col,
  annotation_colors = list(
    Taxon = taxon_colors
  ),
  show_colnames = FALSE,
  fontsize_row = 9,
  border_color = NA,
  main = "Root MAG functional profiles: Sedimenticolaceae vs. Desulfobacterales",
  legend_breaks = c(-2, 0, 2),
  legend_labels = c("Lower", "Mean", "Higher"),
  silent = TRUE
)

png(
  file.path(out_dir, "sedimenticolaceae_vs_desulfobacterales_roots_only_heatmap.png"),
  width = 3900,
  height = 2700,
  res = 300
)
grid::grid.newpage()
grid::grid.draw(heatmap_result$gtable)
dev.off()

pdf(
  file.path(out_dir, "sedimenticolaceae_vs_desulfobacterales_roots_only_heatmap.pdf"),
  width = 13,
  height = 9,
  useDingbats = FALSE
)
grid::grid.newpage()
grid::grid.draw(heatmap_result$gtable)
dev.off()

root_mat <- selected_mat
root_taxonomy <- selected_taxonomy
root_groups <- root_taxonomy$Comparison_group

comparison_rows <- lapply(seq_len(nrow(root_mat)), function(i) {
  x <- root_mat[i, root_groups == "Sedimenticolaceae"]
  y <- root_mat[i, root_groups == "Desulfobacterales"]
  p_value <- tryCatch(
    wilcox.test(x, y, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  pooled_sd <- sd(c(x, y))
  standardized_difference <- if (
    is.na(pooled_sd) || pooled_sd == 0
  ) {
    0
  } else {
    (mean(x) - mean(y)) / pooled_sd
  }
  data.frame(
    Function = rownames(root_mat)[i],
    Sedimenticolaceae_mean = mean(x),
    Desulfobacterales_mean = mean(y),
    Standardized_difference = standardized_difference,
    P_value = p_value,
    stringsAsFactors = FALSE
  )
})
stats <- do.call(rbind, comparison_rows)
stats$Adjusted_P_value <- p.adjust(stats$P_value, method = "BH")
stats$Significance <- ifelse(
  is.na(stats$Adjusted_P_value),
  "Not testable",
  ifelse(stats$Adjusted_P_value < 0.05, "FDR < 0.05", "Not significant")
)
stats <- stats[order(stats$Standardized_difference, decreasing = TRUE), ]

write.csv(
  stats,
  file.path(out_dir, "sedimenticolaceae_vs_desulfobacterales_roots_only_statistics.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    MAG = selected_taxonomy$genome,
    Comparison_group = selected_taxonomy$Comparison_group,
    Compartment = selected_taxonomy$Compartment,
    Phylum = selected_taxonomy$Phylum,
    Class = selected_taxonomy$Class,
    Order = selected_taxonomy$Order,
    Family = selected_taxonomy$Family,
    Genus = selected_taxonomy$Genus,
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "sedimenticolaceae_vs_desulfobacterales_roots_only_MAGs.csv"),
  row.names = FALSE,
  na = "NA"
)

write.csv(
  data.frame(
    Removed_MAG = removed_mags,
    Found_in_input = removed_mags %in% colnames(dat),
    Affects_roots_taxon_comparison =
      removed_mags %in% raw_taxonomy$genome[
        raw_taxonomy$Compartment == "Roots" &
          (
            raw_taxonomy$Family == "Sedimenticolaceae" |
              raw_taxonomy$Order == "Desulfobacterales"
          )
      ],
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "sedimenticolaceae_vs_desulfobacterales_removed_MAGs.csv"),
  row.names = FALSE
)

difference_plot <- ggplot(
  stats,
  aes(
    x = Standardized_difference,
    y = reorder(Function, Standardized_difference),
    color = Significance
  )
) +
  geom_vline(xintercept = 0, color = "#777777", linewidth = 0.45) +
  geom_segment(
    aes(x = 0, xend = Standardized_difference, yend = reorder(Function, Standardized_difference)),
    linewidth = 0.7
  ) +
  geom_point(size = 2.8) +
  scale_color_manual(values = c(
    "FDR < 0.05" = "#D1495B",
    "Not significant" = "#555555",
    "Not testable" = "#888888"
  )) +
  labs(
    title = "Functional differences in root-associated MAGs",
    subtitle = "Positive values indicate enrichment in Sedimenticolaceae",
    x = "Standardized mean difference",
    y = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    axis.text.y = element_text(size = 8)
  )

ggsave(
  file.path(out_dir, "sedimenticolaceae_vs_desulfobacterales_roots_only_functional_differences.png"),
  difference_plot,
  width = 9,
  height = 7,
  dpi = 300
)
ggsave(
  file.path(out_dir, "sedimenticolaceae_vs_desulfobacterales_roots_only_functional_differences.pdf"),
  difference_plot,
  width = 9,
  height = 7,
  device = cairo_pdf
)

pca_input <- t(root_mat)
variable_functions <- apply(pca_input, 2, sd) > 0
pca <- prcomp(pca_input[, variable_functions, drop = FALSE], scale. = TRUE)
pca_variance <- 100 * summary(pca)$importance[2, 1:2]
pca_data <- data.frame(
  MAG = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Taxon = factor(root_groups, levels = group_levels),
  stringsAsFactors = FALSE
)

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = Taxon)) +
  stat_ellipse(aes(fill = Taxon), geom = "polygon", alpha = 0.12, color = NA) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  labs(
    title = "Functional ordination of root-associated MAGs",
    x = sprintf("PC1 (%.1f%%)", pca_variance[1]),
    y = sprintf("PC2 (%.1f%%)", pca_variance[2]),
    color = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(
  file.path(out_dir, "sedimenticolaceae_vs_desulfobacterales_roots_only_PCA.png"),
  pca_plot,
  width = 7,
  height = 5.5,
  dpi = 300
)
ggsave(
  file.path(out_dir, "sedimenticolaceae_vs_desulfobacterales_roots_only_PCA.pdf"),
  pca_plot,
  width = 7,
  height = 5.5,
  device = cairo_pdf
)

heatmap_panel <- wrap_elements(full = heatmap_result$gtable)
combined_panel <- heatmap_panel / (difference_plot | pca_plot) +
  plot_layout(heights = c(1.25, 1)) +
  plot_annotation(
    title = "Functional comparison of Sedimenticolaceae and Desulfobacterales MAGs",
    subtitle = "All panels include root-associated MAGs only",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11)
    )
  )

ggsave(
  file.path(out_dir, "sedimenticolaceae_vs_desulfobacterales_roots_only_combined_panel.png"),
  combined_panel,
  width = 15,
  height = 13,
  dpi = 300
)
ggsave(
  file.path(out_dir, "sedimenticolaceae_vs_desulfobacterales_roots_only_combined_panel.pdf"),
  combined_panel,
  width = 15,
  height = 13,
  device = cairo_pdf
)

message(
  "Figures written to ", normalizePath(out_dir),
  ". Selected MAGs: ", ncol(selected_mat),
  "; Roots comparison: ", ncol(root_mat),
  "; FDR-significant functions: ", sum(stats$Adjusted_P_value < 0.05, na.rm = TRUE)
)


==========================================================================================
# SCRIPT: PGPg_finder/plot_leaves_gaba_production.R
==========================================================================================

library(ggplot2)
library(ggrepel)
library(patchwork)
library(scales)

input_file <- "PGPg_finder/tables/summary/normalized_summary_table.txt"
taxonomy_file <- "PGPg_finder/mags_table.csv"
out_dir <- "PGPg_finder/figures/codex_gaba_leaves"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

phylum_colors <- c(
  Bacteroidota_A = "#D98C7C",
  Cyanobacteria = "#8DB255",
  Desulfobacterota = "#7A8FD1",
  Halobacteriota = "#C9A66B",
  Myxococcota_A = "#A681CF",
  Pseudomonadota = "#48A7C5",
  Other = "#B8B8B8"
)

dat <- read.delim(
  input_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(dat) <- trimws(rownames(dat))

taxonomy <- read.csv(
  taxonomy_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
taxonomy <- taxonomy[match(colnames(dat), taxonomy$genome), ]
if (any(is.na(taxonomy$genome))) {
  stop("Some functional-table MAG identifiers were not found in mags_table.csv.")
}

if (!"PHYTOHORMONE-GABA" %in% rownames(dat)) {
  stop("PHYTOHORMONE-GABA was not found in the normalized summary table.")
}

gaba <- data.frame(
  MAG = colnames(dat),
  GABA = as.numeric(dat["PHYTOHORMONE-GABA", ]),
  Compartment = taxonomy$Compartment,
  Phylum = taxonomy$Phylum,
  Class = taxonomy$Class,
  Order = taxonomy$Order,
  Family = taxonomy$Family,
  Genus = taxonomy$Genus,
  Completeness = taxonomy$Completeness,
  Contamination = taxonomy$Contamination,
  MAG_size = taxonomy$mag_size,
  stringsAsFactors = FALSE
)
gaba <- gaba[gaba$Compartment == "Leaves", , drop = FALSE]
gaba$Phylum_display <- gaba$Phylum
gaba$Phylum_display[gaba$Phylum_display == "Cyanobacteriota"] <- "Cyanobacteria"
gaba$Phylum_display[
  is.na(gaba$Phylum_display) |
    !gaba$Phylum_display %in% names(phylum_colors)
] <- "Other"
gaba$Family_display <- ifelse(is.na(gaba$Family), "Unclassified family", gaba$Family)
gaba$Genus_display <- ifelse(is.na(gaba$Genus), "Unclassified genus", gaba$Genus)
gaba$Taxon_label <- paste(gaba$Family_display, gaba$Genus_display, sep = " / ")
gaba$MAG_rank <- rank(-gaba$GABA, ties.method = "first")
gaba <- gaba[order(gaba$GABA), ]
gaba$MAG_order <- factor(gaba$MAG, levels = gaba$MAG)

write.csv(
  gaba[order(-gaba$GABA), ],
  file.path(out_dir, "leaves_gaba_MAG_values.csv"),
  row.names = FALSE,
  na = "NA"
)

family_summary <- aggregate(
  GABA ~ Phylum_display + Family_display,
  data = gaba,
  FUN = function(x) c(mean = mean(x), max = max(x), n = length(x))
)
family_summary <- do.call(data.frame, family_summary)
names(family_summary) <- c("Phylum", "Family", "Mean_GABA", "Max_GABA", "MAG_count")
family_summary <- family_summary[order(family_summary$Mean_GABA), ]
family_summary$Family_order <- factor(family_summary$Family, levels = family_summary$Family)

write.csv(
  family_summary[order(-family_summary$Mean_GABA), ],
  file.path(out_dir, "leaves_gaba_family_summary.csv"),
  row.names = FALSE,
  na = "NA"
)

rank_plot <- ggplot(gaba, aes(x = MAG_order, y = GABA, fill = Phylum_display)) +
  geom_col(width = 0.78, color = "grey25", linewidth = 0.12) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = phylum_colors, drop = FALSE) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.001),
    expand = expansion(mult = c(0, 0.16))
  ) +
  labs(
    title = "GABA production potential in leaf-associated MAGs",
    subtitle = "Normalized PHYTOHORMONE-GABA values from PGPg_finder",
    x = NULL,
    y = "Normalized GABA value",
    fill = "Phylum"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 7.5),
    legend.position = "right",
    plot.margin = margin(8, 30, 8, 8)
  )

phylum_plot <- ggplot(gaba, aes(x = Phylum_display, y = GABA, fill = Phylum_display)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.35, color = "grey25") +
  geom_jitter(aes(size = Completeness), width = 0.12, alpha = 0.88, shape = 21, color = "grey20") +
  scale_fill_manual(values = phylum_colors, drop = FALSE) +
  scale_size_continuous(range = c(2.2, 5.5)) +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  labs(
    title = "Distribution by phylum",
    x = NULL,
    y = "Normalized GABA value",
    size = "Completeness (%)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "right"
  ) +
  guides(fill = "none")

family_plot <- ggplot(
  family_summary,
  aes(x = Family_order, y = Mean_GABA, fill = Phylum)
) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.12) +
  geom_text(
    aes(label = paste0("n=", MAG_count)),
    hjust = -0.12,
    size = 3
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = phylum_colors, drop = FALSE) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.001),
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "Mean GABA potential by family",
    x = NULL,
    y = "Mean normalized GABA value",
    fill = "Phylum"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    plot.margin = margin(8, 28, 8, 8)
  )

heatmap_data <- gaba[order(gaba$GABA), ]
heatmap_data$y <- "PHYTOHORMONE-GABA"
gaba_heatmap <- ggplot(heatmap_data, aes(x = MAG_order, y = y, fill = GABA)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#B2182B",
    labels = label_number(accuracy = 0.001)
  ) +
  labs(
    title = "GABA production heatmap across leaf MAGs",
    x = "Leaf-associated MAGs ordered by GABA value",
    y = NULL,
    fill = "Normalized\nGABA"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_blank()
  )

combined_panel <- (rank_plot | (phylum_plot / family_plot)) /
  gaba_heatmap +
  plot_layout(heights = c(1, 0.24), widths = c(1.25, 1)) +
  plot_annotation(
    title = "Exploration of GABA production potential in leaf-associated MAGs",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

ggsave(
  file.path(out_dir, "leaves_gaba_ranked_MAGs.png"),
  rank_plot,
  width = 9.5,
  height = 7.2,
  dpi = 300
)
ggsave(
  file.path(out_dir, "leaves_gaba_ranked_MAGs.pdf"),
  rank_plot,
  width = 9.5,
  height = 7.2,
  device = cairo_pdf
)

ggsave(
  file.path(out_dir, "leaves_gaba_taxonomic_summary.png"),
  phylum_plot / family_plot,
  width = 8.5,
  height = 9,
  dpi = 300
)
ggsave(
  file.path(out_dir, "leaves_gaba_taxonomic_summary.pdf"),
  phylum_plot / family_plot,
  width = 8.5,
  height = 9,
  device = cairo_pdf
)

ggsave(
  file.path(out_dir, "leaves_gaba_heatmap.png"),
  gaba_heatmap,
  width = 11,
  height = 2.4,
  dpi = 300
)
ggsave(
  file.path(out_dir, "leaves_gaba_heatmap.pdf"),
  gaba_heatmap,
  width = 11,
  height = 2.4,
  device = cairo_pdf
)

ggsave(
  file.path(out_dir, "leaves_gaba_combined_panel.png"),
  combined_panel,
  width = 15,
  height = 10,
  dpi = 300
)
ggsave(
  file.path(out_dir, "leaves_gaba_combined_panel.pdf"),
  combined_panel,
  width = 15,
  height = 10,
  device = cairo_pdf
)

message(
  "Leaves GABA figures written to: ", normalizePath(out_dir),
  ". Leaf MAGs: ", nrow(gaba),
  "; MAGs with GABA > 0: ", sum(gaba$GABA > 0),
  "; max GABA: ", max(gaba$GABA)
)


==========================================================================================
# SCRIPT: PGPg_finder/plot_leaves_gaba_deep_levels.R
==========================================================================================

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

base_dir <- "/Users/lwmendes/Documents/Codex/Mangrove/PGPg_finder"
gene_file <- file.path(base_dir, "tables", "gene_counts_sum_with_combined_pathways.txt")
mag_file <- file.path(base_dir, "mags_table.csv")
out_dir <- Sys.getenv("PGPG_OUT_DIR", file.path(base_dir, "figures", "codex_gaba_leaves_deep"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

removed_mags_env <- Sys.getenv("PGPG_REMOVE_MAGS", "")
removed_mags <- if (nzchar(removed_mags_env)) {
  trimws(strsplit(removed_mags_env, ",", fixed = TRUE)[[1]])
} else {
  character(0)
}

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold", color = "black"),
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0),
      plot.tag = element_text(face = "bold", size = base_size + 2)
    )
}

read_pgpt_table <- function(path) {
  header <- strsplit(readLines(path, n = 2)[2], "\t", fixed = TRUE)[[1]]
  header[1] <- "PGPT"
  dat <- read.delim(
    path,
    skip = 2,
    header = FALSE,
    sep = "\t",
    check.names = FALSE,
    comment.char = "",
    stringsAsFactors = FALSE
  )
  names(dat) <- header
  dat
}

clean_na <- function(x) {
  x[is.na(x) | x == "" | x == "NA"] <- "Unclassified"
  x
}

gene_table <- read_pgpt_table(gene_file)
mags <- read.csv(mag_file, stringsAsFactors = FALSE, check.names = FALSE)

gaba_rows <- grepl(
  "PHYTOHORMONE-GAMMA-AMINOBUTYRIC_ACID\\|GABA_PRODUCTION",
  gene_table$taxonomy
)
gaba_genes <- gene_table[gaba_rows, , drop = FALSE]

gaba_genes$Deep_category <- ifelse(
  grepl("PHYTOHORMONE-GABA_BIOSYNTHESIS", gaba_genes$taxonomy),
  "GABA biosynthesis",
  ifelse(
    grepl("PHYTOHORMONE-GABA_CONVERSION", gaba_genes$taxonomy),
    "GABA conversion",
    ifelse(
      grepl("PHYTOHORMONE-GABA_TRANSPORT", gaba_genes$taxonomy),
      "GABA transport",
      "Other GABA metabolism"
    )
  )
)
gaba_genes$Gene_label <- sub("^.*# PGPT;", "", gaba_genes$taxonomy)
gaba_genes$Gene_label <- sub("^PGPT[0-9]+-", "", gaba_genes$Gene_label)
gaba_genes$Gene_label <- gsub("\\|", "/", gaba_genes$Gene_label)
gaba_genes$Gene <- sub("-K[0-9]+$", "", gaba_genes$Gene_label)
gaba_genes$KO <- ifelse(
  grepl("-K[0-9]+$", gaba_genes$Gene_label),
  sub("^.*-(K[0-9]+)$", "\\1", gaba_genes$Gene_label),
  NA_character_
)
gaba_genes$Display_gene <- paste0(gaba_genes$Gene, " (", gaba_genes$KO, ")")

leaf_mags <- mags$genome[mags$Compartment == "Leaves" & mags$genome %in% names(gaba_genes)]
leaf_mags <- setdiff(leaf_mags, removed_mags)
leaf_tax <- mags[mags$genome %in% leaf_mags, , drop = FALSE]
leaf_tax$Phylum <- clean_na(leaf_tax$Phylum)
leaf_tax$Family <- clean_na(leaf_tax$Family)
leaf_tax$Genus <- clean_na(leaf_tax$Genus)
leaf_tax$Taxon_label <- paste0(leaf_tax$Family, " / ", leaf_tax$Genus)

long <- do.call(
  rbind,
  lapply(seq_len(nrow(gaba_genes)), function(i) {
    data.frame(
      PGPT = gaba_genes$PGPT[i],
      Gene = gaba_genes$Gene[i],
      KO = gaba_genes$KO[i],
      Display_gene = gaba_genes$Display_gene[i],
      Deep_category = gaba_genes$Deep_category[i],
      MAG = leaf_mags,
      Count = as.numeric(gaba_genes[i, leaf_mags]),
      stringsAsFactors = FALSE
    )
  })
)
long <- merge(long, leaf_tax, by.x = "MAG", by.y = "genome", all.x = TRUE)
long$Short_category <- ifelse(
  long$Deep_category == "GABA biosynthesis",
  "Biosynthesis",
  ifelse(
    long$Deep_category == "GABA conversion",
    "Conversion",
    ifelse(long$Deep_category == "GABA transport", "Transport", "Other")
  )
)
long$Heatmap_gene <- paste0(long$Short_category, ": ", long$Display_gene)

mag_totals <- aggregate(Count ~ MAG + Phylum + Family + Genus + Taxon_label, long, sum)
mag_totals <- mag_totals[order(-mag_totals$Count, mag_totals$MAG), ]
mag_order <- mag_totals$MAG

gene_summary <- aggregate(
  Count ~ Display_gene + Deep_category,
  long,
  function(x) c(total = sum(x), prevalence = sum(x > 0))
)
gene_summary <- do.call(data.frame, gene_summary)
names(gene_summary) <- c("Display_gene", "Deep_category", "Total_count", "Prevalence")
gene_summary <- gene_summary[order(gene_summary$Deep_category, -gene_summary$Prevalence, -gene_summary$Total_count), ]
gene_order <- rev(gene_summary$Display_gene)
gene_summary$Short_category <- ifelse(
  gene_summary$Deep_category == "GABA biosynthesis",
  "Biosynthesis",
  ifelse(
    gene_summary$Deep_category == "GABA conversion",
    "Conversion",
    ifelse(gene_summary$Deep_category == "GABA transport", "Transport", "Other")
  )
)
gene_summary$Heatmap_gene <- paste0(gene_summary$Short_category, ": ", gene_summary$Display_gene)
heatmap_gene_order <- rev(gene_summary$Heatmap_gene)

category_summary <- aggregate(Count ~ MAG + Deep_category, long, sum)
category_summary$MAG <- factor(category_summary$MAG, levels = rev(mag_order))

category_colors <- c(
  "GABA biosynthesis" = "#287271",
  "GABA conversion" = "#E76F51",
  "GABA transport" = "#6D597A",
  "Other GABA metabolism" = "#8D99AE"
)

long$MAG <- factor(long$MAG, levels = mag_order)
long$Display_gene <- factor(long$Display_gene, levels = gene_order)
long$Heatmap_gene <- factor(long$Heatmap_gene, levels = heatmap_gene_order)
gene_summary$Display_gene <- factor(gene_summary$Display_gene, levels = gene_order)

heatmap_plot <- ggplot(long, aes(MAG, Heatmap_gene, fill = Count)) +
  geom_tile(color = "grey96", linewidth = 0.25) +
  scale_fill_gradient(
    low = "#F6F7F7",
    high = "#B2182B",
    breaks = pretty_breaks(4),
    name = "Gene count"
  ) +
  labs(x = NULL, y = NULL, title = "Gene-level GABA metabolism in leaf MAGs") +
  theme_pub(9) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 6.5),
    axis.text.y = element_text(size = 7),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

stacked_plot <- ggplot(category_summary, aes(MAG, Count, fill = Deep_category)) +
  geom_col(width = 0.75, color = "grey98", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = category_colors, name = "Deep category") +
  labs(x = NULL, y = "Total gene count", title = "Contribution of deep categories") +
  theme_pub(9) +
  theme(
    axis.text.y = element_text(size = 6.5),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

gene_summary$Deep_category <- factor(gene_summary$Deep_category, levels = names(category_colors))
prevalence_plot <- ggplot(gene_summary, aes(Prevalence, Display_gene, fill = Deep_category)) +
  geom_col(width = 0.72) +
  geom_text(
    aes(label = paste0("n=", Total_count)),
    hjust = -0.08,
    size = 2.3,
    color = "black"
  ) +
  scale_fill_manual(values = category_colors, guide = "none") +
  scale_x_continuous(
    limits = c(0, max(gene_summary$Prevalence) * 1.22),
    breaks = pretty_breaks(5)
  ) +
  labs(
    x = "Number of leaf MAGs with gene count > 0",
    y = NULL,
    title = "Gene prevalence across leaf MAGs"
  ) +
  theme_pub(9) +
  theme(axis.text.y = element_text(size = 7))

tax_bar <- ggplot(mag_totals, aes(x = factor(MAG, levels = mag_order), y = 1, fill = Phylum)) +
  geom_tile(color = "grey98", linewidth = 0.2) +
  scale_fill_brewer(palette = "Set2", name = "Phylum") +
  labs(x = NULL, y = NULL, title = "MAG phylum order") +
  theme_pub(8) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 8)
  )

combined_panel <- heatmap_plot / tax_bar / (stacked_plot | prevalence_plot) +
  plot_layout(heights = c(4.2, 0.45, 2.7)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(out_dir, "leaves_gaba_deep_gene_heatmap.png"), heatmap_plot, width = 12, height = 5.8, dpi = 300)
ggsave(file.path(out_dir, "leaves_gaba_deep_gene_heatmap.pdf"), heatmap_plot, width = 12, height = 5.8)
ggsave(file.path(out_dir, "leaves_gaba_deep_category_stacked_bar.png"), stacked_plot, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "leaves_gaba_deep_category_stacked_bar.pdf"), stacked_plot, width = 8, height = 6)
ggsave(file.path(out_dir, "leaves_gaba_deep_gene_prevalence.png"), prevalence_plot, width = 7.5, height = 5, dpi = 300)
ggsave(file.path(out_dir, "leaves_gaba_deep_gene_prevalence.pdf"), prevalence_plot, width = 7.5, height = 5)
ggsave(file.path(out_dir, "leaves_gaba_deep_combined_panel.png"), combined_panel, width = 13, height = 10.5, dpi = 300)
ggsave(file.path(out_dir, "leaves_gaba_deep_combined_panel.pdf"), combined_panel, width = 13, height = 10.5)

write.csv(long, file.path(out_dir, "leaves_gaba_deep_gene_counts_long.csv"), row.names = FALSE)
write.csv(mag_totals, file.path(out_dir, "leaves_gaba_deep_MAG_totals.csv"), row.names = FALSE)
write.csv(gene_summary, file.path(out_dir, "leaves_gaba_deep_gene_summary.csv"), row.names = FALSE)
write.csv(
  data.frame(
    Removed_MAG = removed_mags,
    Found_in_input = removed_mags %in% names(gaba_genes),
    Removed_from_leaf_GABA_panel = removed_mags %in% mags$genome[mags$Compartment == "Leaves"],
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "leaves_gaba_deep_removed_MAGs.csv"),
  row.names = FALSE
)

cat("Wrote GABA deep-level figures and tables to:", out_dir, "\n")
cat("Direct GABA genes:", nrow(gaba_genes), "\n")
cat("Leaf MAGs:", length(leaf_mags), "\n")
cat("Top leaf MAGs by direct GABA gene count:\n")
print(head(mag_totals[, c("MAG", "Phylum", "Family", "Genus", "Count")], 8), row.names = FALSE)
cat("Gene summary:\n")
print(gene_summary[order(-gene_summary$Total_count), ], row.names = FALSE)

