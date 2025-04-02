
# --------------------------
# 1. Load Required Libraries
# --------------------------
library(tidyverse)
library(glue)
library(interactions)
library(patchwork)
library(ggplot2)
library(ggpubr)   # for stat_cor(), if needed
library(future)
library(future.apply)
library(parallelly)
library(knitr)
library(GSEABase)  # for getGmt()

# --------------------------
# 2. Load Config & Set Parallel Processing
# --------------------------
config <- config::get(file = "config.yml")
source(config$xena_script)  # Load custom functions

workers <- as.integer(config$workers)
options(
  parallelly.availableCores.system = workers,
  parallelly.hard.CPU.limit = workers
)
plan(multisession, workers = workers)

# --------------------------
# 3. Set Directories & Knit Root
# --------------------------
output_dir    <- config$output_dir
tcga_data_dir <- config$tcga_data_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
# Set the knit root to the TCGA data directory
knitr::opts_knit$set(root.dir = tcga_data_dir)

# --------------------------
# 4. Load GSDM Expression Data
# --------------------------
load(config$gsdm_file)  # loads object: GSDMs.PANCAN.df

# --------------------------
# 5. Modified Function: Read Gene Expression Data
# --------------------------
# The function now accepts a 'base_dir' so that if the filename is not found 
# in its current form, it will look relative to tcga_data_dir.
read_file_id_on_col1 <- function(id_on_col1, filenames, base_dir) {
  dfLst <- lapply(filenames, function(filename) {
    full_filename <- filename
    if (!file.exists(full_filename)) {
      # If the file is not found, try prepending the base directory.
      full_filename <- file.path(base_dir, filename)
    }
    if (!file.exists(full_filename)) {
      warning("File not found: ", full_filename)
      return(NULL)
    }
    message("Processing file: ", full_filename)
    col1 <- data.table::fread(full_filename, select = 1)[[1]]
    if (!id_on_col1 %in% col1) {
      warning("Gene ", id_on_col1, " not found in ", full_filename)
      return(NULL)
    }
    row_idx <- which(id_on_col1 == col1)
    if (length(row_idx) == 0) {
      warning("No matching rows found for ", id_on_col1, " in ", full_filename)
      return(NULL)
    }
    df <- read_row(full_filename, row_idx)
    df$File <- full_filename
    df
  })
  dfLst <- Filter(Negate(is.null), dfLst)
  if (length(dfLst) == 0) {
    warning("No data found for gene ", id_on_col1)
    return(NULL)
  }
  do.call(rbind, dfLst)
}

# --------------------------
# 6. Gene Set / Published Signature Processing
# --------------------------


# Option A: Use provided gene list
gene_set <- config$gene_list
message("Using gene list from config: ", paste(gene_set, collapse = ", "))
load(config$gex_file)  # loads object: TCGA.cohorts.GEX_be.tb

# Filter out cohorts COADREAD, LUNG, and GBMLGG
TCGA.cohorts.GEX_be.flt.tb <- TCGA.cohorts.GEX_be.tb %>%
  filter(!Cohort_symbol %in% c("COADREAD", "LUNG", "GBMLGG"))

gene_expr_ldfLst <- lapply(as.character(gene_set), function(gene_symbol) {
  message("Processing gene: ", gene_symbol)
  read_file_id_on_col1(id_on_col1 = gene_symbol, 
                       filenames = TCGA.cohorts.GEX_be.flt.tb$destfiles,
                       base_dir = tcga_data_dir)
})
gene_expr_ldfLst <- Filter(Negate(is.null), gene_expr_ldfLst)
if (length(gene_expr_ldfLst) == 0) {
  stop("No valid gene expression data found for any gene in the gene list.")
}
gene_expr_df_first <- do.call(cbind, lapply(gene_expr_ldfLst, function(x) {
  dplyr::select(x, -File, -Name)
}))
gene_expr_df_first$Sample_name <- gene_expr_ldfLst[[1]]$Name
save(gene_expr_df_first, file = "gene_expr_df.Rda")



message("Using published signature data from: ", config$published_score_df)
load(config$published_score_df)  # loads panimmune_signature.df (or similar)
sample_name_col <- if (!is.null(config$sample_name_column) && config$sample_name_column != "") {
  config$sample_name_column
} else {
  "Sample"
}
gene_expr_df_panimmune <- panimmune_signature.df %>%
  dplyr::select(!!sym(sample_name_col), !!sym(config$published_score_column)) %>%
  dplyr::rename(Sample_name = !!sym(sample_name_col))
gene_set_panimmune <- config$published_score_column
message("Using published score column: ", gene_set_panimmune)



# --------------------------
# 7. Merge Sample Type & GSDM Data
# --------------------------

# Merge expression with panimmune


gene_expr_df <- left_join(gene_expr_df_first, gene_expr_df_panimmune, by = "Sample_name")

load(config$sample_type_file)

master_expr.df <- left_join(gene_expr_df, TCGA.PANCAN.sample_type.df, by = "Sample_name") %>%
  filter(Sample_type == "Primary Tumor")
save(master_expr.df, file = "master_expr.df.Rda")


# --------------------------
# 8. Load & Process TMB Data, Merge with Expression Data
# --------------------------
load(config$tmb_data_file)  # loads object: tcga_all
TMB_APM_TIGS.df <- tcga_all %>%
  mutate(
    nAPM = (APM - min(APM, na.rm = TRUE)) / (max(APM, na.rm = TRUE) - min(APM, na.rm = TRUE)),
    TMB_perMB = TMB_NonsynVariants / 38,
    TIGS = log(TMB_perMB + 1) * nAPM
  )
save(TMB_APM_TIGS.df, file = "TMB_APM_TIGS.df.Rda")


final_master.df <- master_expr.df %>%
  mutate(Tumor_Sample_Barcode = str_sub(Sample_name, 1, 12)) %>%
  left_join(TMB_APM_TIGS.df, by = "Tumor_Sample_Barcode") %>%
  mutate(log10_TMB_perMB = log10(TMB_perMB)) %>%
  filter(!is.na(TMB_perMB) & Sample_type == "Primary Tumor")
save(final_master.df, file = "final_master.noNA_primary.df.Rda")

load("final_master.noNA_primary.df.Rda")

# --------------------------
# Make Cytolytic scores
# --------------------------

genes_cytolytic = c("PRF1","GZMA","GZMB" )
std_gene_cols <- paste0(gene_set, "_std")

if (!all(std_gene_cols %in% colnames(final_master.df))) {
  message("Standardized columns not found; computing standardized expression values.")
  std_list <- future_lapply(gene_set, function(gene) {
    as.numeric(scale(final_master.df[[gene]], center = TRUE, scale = TRUE))
  })
  std_df <- as.data.frame(std_list)
  colnames(std_df) <- std_gene_cols
  final_master.df <- bind_cols(final_master.df, std_df)
}

#CAS.PA
final_master.df["CAS.PA"] <- rowMeans(final_master.df[, c("PRF1_std","GZMA_std")], na.rm = TRUE)


#CAS.PB
final_master.df["CAS.PB"] <- rowMeans(final_master.df[, c("PRF1_std","GZMB_std")], na.rm = TRUE)


#################### GENERATE FIGURES NOW ##########################

####### FIGURE 4b ############
# --------------------------
# 9. Set a Base Plot Theme (White background, no grids)
# --------------------------
base_theme <- theme_minimal(base_family = "sans") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.title.x     = element_text(size = 16, family = "sans"),
    axis.title.y     = element_text(size = 16, family = "sans"),
    plot.title       = element_text(size = 18, family = "sans"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# --------------------------
# 10. Create DFNA5 Scatter Plot 
# --------------------------
p_scatter_dfna5 <- ggplot(final_master.df, aes(x = DFNA5, y = log10_TMB_perMB)) +
  geom_point(shape = 1, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  stat_cor(method = "pearson", label.x.npc = 0.05, label.y.npc = 0.98, family = "sans", label.fill = NA) +
  
  labs(
    x = "DFNA5 Gene Expression",
    y = expression(Log[10] ~ "TMB per MB"),
    title = "DFNA5 vs Log10 TMB per MB"
  ) +
  base_theme

# Save DFNA5 scatter plot
ggsave(filename = file.path(output_dir, "Scatter_DFNA5_vs_log10TMB.png"),
       plot = p_scatter_dfna5, height = 6, width = 8, dpi = 300)
ggsave(filename = file.path(output_dir, "Scatter_DFNA5_vs_log10TMB.pdf"),
       plot = p_scatter_dfna5, height = 6, width = 8, dpi = 300)
ggsave(filename = file.path(output_dir, "Scatter_DFNA5_vs_log10TMB.svg"),
       plot = p_scatter_dfna5, height = 6, width = 8, dpi = 300)


# Compute and print regression slope for DFNA5 scatter plot
model_dfna5 <- lm(log10_TMB_perMB ~ DFNA5, data = final_master.df)
summary_dfna5 <- summary(model_dfna5)
slope_dfna5 <- summary_dfna5$coefficients["DFNA5", "Estimate"]
se_dfna5 <- summary_dfna5$coefficients["DFNA5", "Std. Error"]
print(summary_dfna5)
message("DFNA5 slope (all Primary Tumors): ", round(slope_dfna5, 4), 
        " (Std. Error: ", round(se_dfna5, 4), ")")


# --------------------------
# 10. Create CAS.PA Scatter Plot 
# --------------------------

p_scatter_caspa <- ggplot(final_master.df, aes(x = CAS.PA, y = log10_TMB_perMB)) +
  geom_point(shape = 1, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  stat_cor(method = "pearson", label.x.npc = 0.05, label.y.npc = 0.98, family = "sans", label.fill = NA) +
  
  labs(
    x = "Cytolytic genes (PA) score",
    y = expression(Log[10] ~ "TMB per MB"),
    title = "PA vs Log10 TMB per MB"
  ) +
  base_theme

# Save DFNA5 scatter plot
ggsave(filename = file.path(output_dir, "Scatter_PA_vs_log10TMB.png"),
       plot = p_scatter_caspa, height = 6, width = 8, dpi = 300)
ggsave(filename = file.path(output_dir, "Scatter_PA_vs_log10TMB.pdf"),
       plot = p_scatter_caspa, height = 6, width = 8, dpi = 300)
ggsave(filename = file.path(output_dir, "Scatter_PA_vs_log10TMB.svg"),
       plot = p_scatter_caspa, height = 6, width = 8, dpi = 300)


# Compute and print regression slope for CAS.PA scatter plot
model_pa <- lm(log10_TMB_perMB ~ CAS.PA, data = final_master.df)
summary_pa <- summary(model_pa)
slope_pa <- summary_pa$coefficients["CAS.PA", "Estimate"]
se_pa <- summary_pa$coefficients["CAS.PA", "Std. Error"]
print(summary_pa)
message("PA slope (all Primary Tumors): ", round(slope_pa, 4), 
        " (Std. Error: ", round(se_pa, 4), ")")


# --------------------------
# 13. Subgroup Analysis by DFNA5 Levels
# --------------------------
mean_DFNA5 <- mean(final_master.df$DFNA5, na.rm = TRUE)
sd_DFNA5 <- sd(final_master.df$DFNA5, na.rm = TRUE)
final_master.df <- final_master.df %>%
  mutate(DFNA5_group = case_when(
    DFNA5 <= (mean_DFNA5 - sd_DFNA5) ~ "Low",
    DFNA5 >= (mean_DFNA5 + sd_DFNA5) ~ "High",
    TRUE ~ "Medium"
  )) %>%
  filter(!is.na(DFNA5_group))
final_master.df$DFNA5_group <- factor(final_master.df$DFNA5_group, levels = c("Low", "Medium", "High"))

# save the final dataframe for publication 
save(final_master.df, file = "final_master.df_publication.Rda")

load("final_master.df_publication.Rda")

# Analysis 1: Grouping by DFNA5 (Low vs. High)
#   - Scatter plots of CAS.PA vs. log10 TMB per MB by DFNA5 group.
#   - Interaction model: log10_TMB_perMB ~ CAS.PA * DFNA5_group.
#
# Analysis 2: Grouping by CAS.PA (Low vs. High)
#   - Scatter plots of DFNA5 vs. log10 TMB per MB by CAS.PA group.
#   - Interaction model: log10_TMB_perMB ~ DFNA5 * CASPA_group.
#

# for plotting x-axis

min_val_pa <- min(final_master.df$CAS.PA, na.rm = TRUE)
max_val_pa <- max(final_master.df$CAS.PA, na.rm = TRUE)

min_val_gsdme <- min(final_master.df$DFNA5, na.rm = TRUE)
max_val_gsdme <- max(final_master.df$DFNA5, na.rm = TRUE)

min_val_y <- min(final_master.df$log10_TMB_perMB, na.rm = TRUE)
max_val_y <- max(final_master.df$log10_TMB_perMB, na.rm = TRUE)

# =============================================================================
# Analysis 1: Grouping by DFNA5 (Low vs. High)
# =============================================================================

# Group samples by DFNA5 using mean – 1SD thresholds and keep only Low and High
mean_DFNA5 <- mean(final_master.df$DFNA5, na.rm = TRUE)
sd_DFNA5 <- sd(final_master.df$DFNA5, na.rm = TRUE)
df_DFNA5 <- final_master.df %>%
  mutate(DFNA5_group = case_when(
    DFNA5 <= (mean_DFNA5 - sd_DFNA5) ~ "Low",
    DFNA5 >= (mean_DFNA5 + sd_DFNA5) ~ "High",
    TRUE ~ "Medium"
  )) %>%
  filter(DFNA5_group %in% c("Low", "High"))
df_DFNA5$DFNA5_group <- factor(df_DFNA5$DFNA5_group, levels = c("Low", "High"))

# Create scatter plots: CAS.PA (x) vs. log10 TMB per MB (y) for each DFNA5 group
plots_list_DFNA5 <- list()
for (g in levels(df_DFNA5$DFNA5_group)) {
  group_data <- df_DFNA5 %>% filter(DFNA5_group == g)
  p <- ggplot(group_data, aes(x = CAS.PA, y = log10_TMB_perMB)) +
    geom_point(shape = 1, color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 1) +
    stat_cor(method = "pearson", label.x.npc = 0.05, label.y.npc = 0.98, family = "sans") +
    labs(
      x = "CAS.PA",
      y = expression(Log[10] ~ "TMB per MB"),
      title = paste("DFNA5", g, "group")
    ) +
    base_theme+
    scale_x_continuous(limits = c(min_val_pa, max_val_pa)) +
    scale_y_continuous(limits = c(min_val_y, max_val_y))
  plots_list_DFNA5[[g]] <- p
  
  ggsave(filename = file.path(output_dir, paste0("Figure_CASPA_vs_log10TMB_DFNA5_", g, ".png")),
         plot = p, height = 6, width = 8, dpi = 300)
  
  ggsave(filename = file.path(output_dir, paste0("Figure_CASPA_vs_log10TMB_DFNA5_", g, ".pdf")),
         plot = p, height = 6, width = 8, dpi = 300)
  ggsave(filename = file.path(output_dir, paste0("Figure_CASPA_vs_log10TMB_DFNA5_", g, ".svg")),
         plot = p, height = 6, width = 8, dpi = 300)
}

# Fit an interaction model to test slope differences between Low and High DFNA5 groups:
model_DFNA5 <- lm(log10_TMB_perMB ~ CAS.PA * DFNA5_group, data = df_DFNA5)
summary_DFNA5 <- summary(model_DFNA5)
print(summary_DFNA5)

# Extract the p-value for the interaction term comparing Low (reference) vs. High DFNA5 groups.
p_value_interaction <- summary_DFNA5$coefficients["CAS.PA:DFNA5_groupHigh", "Pr(>|t|)"]

# Print a message summarizing the result.
message("Comparison between Low and High DFNA5 groups for the CAS.PA effect yields a p-value of ", 
        format(p_value_interaction, digits = 3))


# Interpretation: The interaction term (CAS.PA:DFNA5_groupHigh) indicates whether the slope for CAS.PA differs significantly
# between Low (reference) and High DFNA5 groups.

# =============================================================================
# Analysis 2: Grouping by CAS.PA (Low vs. High)
# =============================================================================




# Group samples by CAS.PA using mean – 1SD thresholds and keep only Low and High
mean_CASPA <- mean(final_master.df$CAS.PA, na.rm = TRUE)
sd_CASPA <- sd(final_master.df$CAS.PA, na.rm = TRUE)
df_CASPA <- final_master.df %>%
  mutate(CASPA_group = case_when(
    CAS.PA <= (mean_CASPA - sd_CASPA) ~ "Low",
    CAS.PA >= (mean_CASPA + sd_CASPA) ~ "High",
    TRUE ~ "Medium"
  )) %>%
  filter(CASPA_group %in% c("Low", "High"))
df_CASPA$CASPA_group <- factor(df_CASPA$CASPA_group, levels = c("Low", "High"))

# Create scatter plots: DFNA5 (x) vs. log10 TMB per MB (y) for each CAS.PA group
plots_list_CASPA <- list()
for (g in levels(df_CASPA$CASPA_group)) {
  group_data <- df_CASPA %>% filter(CASPA_group == g)
  p <- ggplot(group_data, aes(x = DFNA5, y = log10_TMB_perMB)) +
    geom_point(shape = 1, color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 1) +
    stat_cor(method = "pearson", label.x.npc = 0.05, label.y.npc = 0.98, family = "sans") +
    labs(
      x = "GSDME Expression (DFNA5)",
      y = expression(Log[10] ~ "TMB per MB"),
      title = paste("CAS.PA", g, "group")
    ) +
    base_theme+
    scale_x_continuous(limits = c(min_val_gsdme, max_val_gsdme)) +
    scale_y_continuous(limits = c(min_val_y, max_val_y))
  plots_list_CASPA[[g]] <- p
  
  ggsave(filename = file.path(output_dir, paste0("Figure_DFNA5_vs_log10TMB_CASPA_", g, ".png")),
         plot = p, height = 6, width = 8, dpi = 300)
  ggsave(filename = file.path(output_dir, paste0("Figure_DFNA5_vs_log10TMB_CASPA_", g, ".pdf")),
         plot = p, height = 6, width = 8, dpi = 300)
  ggsave(filename = file.path(output_dir, paste0("Figure_DFNA5_vs_log10TMB_CASPA_", g, ".svg")),
         plot = p, height = 6, width = 8, dpi = 300)
}

# Fit an interaction model to test slope differences between Low and High CAS.PA groups:
model_CASPA <- lm(log10_TMB_perMB ~ DFNA5 * CASPA_group, data = df_CASPA)
summary_CASPA <- summary(model_CASPA)
print(summary_CASPA)


# Extract the p-value for the interaction term comparing Low (reference) vs. High DFNA5 groups.
p_value_interaction <- summary_CASPA$coefficients["DFNA5:CASPA_groupHigh", "Pr(>|t|)"]

# Print a message summarizing the result.
message("Comparison between Low and High CAS.PA  groups for the DFNA5 effect yields a p-value of ", 
        format(p_value_interaction, digits = 3))


# Interpretation: The interaction term (DFNA5:CASPA_groupHigh) shows if the slope for DFNA5 differs significantly
# between the Low (reference) and High CAS.PA groups.


# --------------------------
# Interaction Model and Plots for CAS.PA
# (Model: log10_TMB_perMB ~ CAS.PA * DFNA5)
# --------------------------
model_formula_pa <- as.formula(glue("log10_TMB_perMB ~ CAS.PA * DFNA5"))
model_interaction_pa <- aov(model_formula_pa, data = final_master.df)
print(summary(model_interaction_pa))

# Determine DFNA5 quantiles for vertical lines in the J-N plot
DFNA5_q1_q3 <- quantile(final_master.df$DFNA5, na.rm = TRUE)[c(2, 4)]

# Interaction Plot with points for CAS.PA
p_int_pa <- interact_plot(
  model = model_interaction_pa,
  pred = "CAS.PA",   # Passing as a string
  modx = "DFNA5",    # Passing as a string
  plot.points = TRUE
) +
  labs(
    x = "CAS.PA",
    y = bquote(Log[10] ~ "TMB per MB"),
    title = "Interaction Plot: TMB vs CAS.PA \n by DFNA5 (Primary Tumors)"
  ) +
  base_theme

# Johnson–Neyman Plot for CAS.PA
p_jn_pa <- johnson_neyman(
  model = model_interaction_pa,
  pred = "CAS.PA",
  modx = "DFNA5",
  alpha = 0.05
)
p_jn_pa$plot <- p_jn_pa$plot +
  geom_vline(xintercept = DFNA5_q1_q3, linetype = "dashed", color = "grey50", size = 0.5) +
  ggtitle("Johnson–Neyman Plot")

# Save the Johnson–Neyman plot for CAS.PA alone
ggsave(filename = file.path(output_dir, "Figure_CASPA_JohnsonNeyman_Primary.png"),
       plot = p_jn_pa$plot, height = 5, width = 7, dpi = 300)
ggsave(filename = file.path(output_dir, "Figure_CASPA_JohnsonNeyman_Primary.pdf"),
       plot = p_jn_pa$plot, height = 5, width = 7, dpi = 300)
ggsave(filename = file.path(output_dir, "Figure_CASPA_JohnsonNeyman_Primary.svg"),
       plot = p_jn_pa$plot, height = 5, width = 7, dpi = 300)

# Combine interaction and J-N plots for CAS.PA
combined_plot_pa <- wrap_plots(p_int_pa, p_jn_pa$plot, nrow = 1) +
  plot_annotation(
    title = "TMB vs CAS.PA by DFNA5 (Primary Tumors)",
    subtitle = "Interaction & Johnson–Neyman Plots",
    caption = "Model: log10(TMB per MB) ~ CAS.PA * DFNA5"
  )

# Save the combined plot for CAS.PA
ggsave(filename = file.path(output_dir, "Figure_CASPA_Interaction_Primary.png"),
       plot = combined_plot_pa, height = 5, width = 11, dpi = 300)
ggsave(filename = file.path(output_dir, "Figure_CASPA_Interaction_Primary.pdf"),
       plot = combined_plot_pa, height = 5, width = 11, dpi = 300)
ggsave(filename = file.path(output_dir, "Figure_CASPA_Interaction_Primary.svg"),
       plot = combined_plot_pa, height = 5, width = 11, dpi = 300)

# --------------------------
# Interaction Model and Plots for CAS.PB
# (Model: log10_TMB_perMB ~ CAS.PB * DFNA5)
# --------------------------
model_formula_pb <- as.formula(glue("log10_TMB_perMB ~ CAS.PB * DFNA5"))
model_interaction_pb <- aov(model_formula_pb, data = final_master.df)
print(summary(model_interaction_pb))

# Interaction Plot with points for CAS.PB
p_int_pb <- interact_plot(
  model = model_interaction_pb,
  pred = "CAS.PB",
  modx = "DFNA5",
  plot.points = TRUE
) +
  labs(
    x = "CAS.PB",
    y = bquote(Log[10] ~ "TMB per MB"),
    title = "Interaction Plot: TMB vs CAS.PB \n by DFNA5 (Primary Tumors)"
  ) +
  base_theme

# Johnson–Neyman Plot for CAS.PB
p_jn_pb <- johnson_neyman(
  model = model_interaction_pb,
  pred = "CAS.PB",
  modx = "DFNA5",
  alpha = 0.05
)
p_jn_pb$plot <- p_jn_pb$plot +
  geom_vline(xintercept = DFNA5_q1_q3, linetype = "dashed", color = "grey50", size = 0.5) +
  ggtitle("Johnson–Neyman Plot")

# Save the Johnson–Neyman plot for CAS.PB alone
ggsave(filename = file.path(output_dir, "Figure_CASPB_JohnsonNeyman_Primary.png"),
       plot = p_jn_pb$plot, height = 5, width = 7, dpi = 300)
ggsave(filename = file.path(output_dir, "Figure_CASPB_JohnsonNeyman_Primary.pdf"),
       plot = p_jn_pb$plot, height = 5, width = 7, dpi = 300)
ggsave(filename = file.path(output_dir, "Figure_CASPB_JohnsonNeyman_Primary.svg"),
       plot = p_jn_pb$plot, height = 5, width = 7, dpi = 300)

# Combine interaction and J-N plots for CAS.PB
combined_plot_pb <- wrap_plots(p_int_pb, p_jn_pb$plot, nrow = 1) +
  plot_annotation(
    title = "TMB vs CAS.PB by DFNA5 (Primary Tumors)",
    subtitle = "Interaction & Johnson–Neyman Plots",
    caption = "Model: log10(TMB per MB) ~ CAS.PB * DFNA5"
  )

# Save the combined plot for CAS.PB
ggsave(filename = file.path(output_dir, "Figure_CASPB_Interaction_Primary.png"),
       plot = combined_plot_pb, height = 5, width = 11, dpi = 300)
ggsave(filename = file.path(output_dir, "Figure_CASPB_Interaction_Primary.pdf"),
       plot = combined_plot_pb, height = 5, width = 11, dpi = 300)
ggsave(filename = file.path(output_dir, "Figure_CASPB_Interaction_Primary.svg"),
       plot = combined_plot_pb, height = 5, width = 11, dpi = 300)

# --------------------------
# Interaction Model and Plots for T.cells.CD8
# (Model: log10_TMB_perMB ~ T.cells.CD8 * DFNA5)
# --------------------------
# Note: Use backticks if the column name has special characters.
model_formula_cd8 <- as.formula(glue("log10_TMB_perMB ~ `T.cells.CD8` * DFNA5"))
model_interaction_cd8 <- aov(model_formula_cd8, data = final_master.df)
print(summary(model_interaction_cd8))

# Interaction Plot with points for T.cells.CD8
p_int_cd8 <- interact_plot(
  model = model_interaction_cd8,
  pred = "T.cells.CD8",
  modx = "DFNA5",
  plot.points = TRUE
) +
  labs(
    x = "T.cells.CD8",
    y = bquote(Log[10] ~ "TMB per MB"),
    title = "Interaction Plot: TMB vs T.cells.CD8 \n by DFNA5 (Primary Tumors)"
  ) +
  base_theme

# Johnson–Neyman Plot for T.cells.CD8
p_jn_cd8 <- johnson_neyman(
  model = model_interaction_cd8,
  pred = "T.cells.CD8",
  modx = "DFNA5",
  alpha = 0.05
)
p_jn_cd8$plot <- p_jn_cd8$plot +
  geom_vline(xintercept = DFNA5_q1_q3, linetype = "dashed", color = "grey50", size = 0.5) +
  ggtitle("Johnson–Neyman Plot")

# Save the Johnson–Neyman plot for T.cells.CD8 alone
ggsave(filename = file.path(output_dir, "Figure_TcellsCD8_JohnsonNeyman_Primary.png"),
       plot = p_jn_cd8$plot, height = 5, width = 7, dpi = 300)
ggsave(filename = file.path(output_dir, "Figure_TcellsCD8_JohnsonNeyman_Primary.pdf"),
       plot = p_jn_cd8$plot, height = 5, width = 7, dpi = 300)
ggsave(filename = file.path(output_dir, "Figure_TcellsCD8_JohnsonNeyman_Primary.svg"),
       plot = p_jn_cd8$plot, height = 5, width = 7, dpi = 300)

# Combine interaction and J-N plots for T.cells.CD8
combined_plot_cd8 <- wrap_plots(p_int_cd8, p_jn_cd8$plot, nrow = 1) +
  plot_annotation(
    title = "TMB vs T.cells.CD8 by DFNA5 (Primary Tumors)",
    subtitle = "Interaction & Johnson–Neyman Plots",
    caption = "Model: log10(TMB per MB) ~ T.cells.CD8 * DFNA5"
  )

# Save the combined plot for T.cells.CD8
ggsave(filename = file.path(output_dir, "Figure_TcellsCD8_Interaction_Primary.png"),
       plot = combined_plot_cd8, height = 5, width = 11, dpi = 300)
ggsave(filename = file.path(output_dir, "Figure_TcellsCD8_Interaction_Primary.pdf"),
       plot = combined_plot_cd8, height = 5, width = 11, dpi = 300)
ggsave(filename = file.path(output_dir, "Figure_TcellsCD8_Interaction_Primary.svg"),
       plot = combined_plot_cd8, height = 5, width = 11, dpi = 300)


################# SOME OTHER PLOTS GASDERMIN FAMILIES #################3

library(dplyr)
library(tidyr)
library(ggplot2)

# Ensure the columns exist (adjust if your column names differ)
genes <- c("GSDMA", "GSDMB", "GSDMC", "GSDMD")

# Pivot the data from wide to long format so we have a single "Gene" column and "Expression" column.
df_long <- final_master.df %>%
  dplyr::select(all_of(genes), log10_TMB_perMB) %>%
  pivot_longer(
    cols = all_of(genes),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  mutate(Gene = factor(Gene, levels = genes))

# ------------------------------------------------------------------------------
# Create the facet plot with removed grid lines and hollow circles
p_gsdm <- ggplot(df_long, aes(x = Expression, y = log10_TMB_perMB)) +
  geom_point(alpha = 0.8, shape = 1, color = "black") +  # hollow circles with black boundary
  geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 1) +
  stat_cor(
    method = "pearson",
    label.x.npc = 0.05,      # position of correlation label on x-axis
    label.y.npc = 0.95,      # position of correlation label on y-axis
    size = 4
  ) +
  facet_wrap(~ Gene, nrow = 2, ncol = 2, scales = "free_x") +
  labs(
    x = "Gene Expression",
    y = expression(log[10] ~ "TMB"),
    title = "GSDM Family Gene Expression vs. log10 TMB"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    panel.grid.major = element_blank(),  # remove major grid lines
    panel.grid.minor = element_blank(),  # remove minor grid lines
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )

# ------------------------------------------------------------------------------
# Save the plot
output_dir <- "GSDM_family"  # change to your desired folder
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(
  filename = file.path(output_dir, "GSDM_Family_vs_TMB.png"),
  plot = p_gsdm, width = 10, height = 8, dpi = 300
)
ggsave(
  filename = file.path(output_dir, "GSDM_Family_vs_TMB.pdf"),
  plot = p_gsdm, width = 10, height = 8, dpi = 300
)
ggsave(
  filename = file.path(output_dir, "GSDM_Family_vs_TMB.svg"),
  plot = p_gsdm, width = 10, height = 8, dpi = 300
)

message("Facet plot for GSDMA/B/C/D vs. log10 TMB saved successfully!")
