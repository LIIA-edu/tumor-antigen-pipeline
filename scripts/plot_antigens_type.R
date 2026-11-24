#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(optparse)
  library(stringr)
  library(forcats)
  library(tools)
})

options(error = function() {
  traceback(20)
  quit(status = 1)
})

# =========================
# ARGUMENTS
# =========================

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Input CSV file"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "specificity_groups_by_type_heatmap.svg",
    help = "Output heatmap SVG"
  ),
  make_option(
    c("--merged_output"),
    type = "character",
    default = "merged_summed_genes.csv",
    help = "Output CSV with merged rows"
  ),
  make_option(
    c("--venn_output"),
    type = "character",
    default = "cta_malignant_vs_others_venn.svg",
    help = "Output CTA Venn diagram SVG"
  ),
  make_option(
    c("--barplot_output"),
    type = "character",
    default = "cta_enriched_categories_barplot.svg",
    help = "Output CTA barplot SVG"
  ),
  make_option(
    c("--malignant_pattern"),
    type = "character",
    default = "malignant",
    help = "Regex used to identify malignant cell labels"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# =========================
# VALIDATE INPUT
# =========================

if (is.null(opt$input) || !file.exists(opt$input)) {
  stop("Input CSV does not exist. Provide a valid --input path.")
}

# =========================
# LOAD DATA
# =========================

df <- read_csv(
  opt$input,
  show_col_types = FALSE,
  na = c("", "NA", "na", "Na", "N/A", "n/a")
)

required_columns <- c(
  "celltypist_labels",
  "type",
  "specificity_profile",
  "gene_name"
)

missing_columns <- setdiff(required_columns, names(df))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# =========================
# CLEAN LABELS
# =========================

df <- df %>%
  mutate(
    celltypist_labels = as.character(celltypist_labels),
    celltypist_labels = str_trim(celltypist_labels),
    celltypist_labels = replace_na(celltypist_labels, "Unknown"),
    celltypist_labels = if_else(
      celltypist_labels == "",
      "Unknown",
      celltypist_labels
    ),
    type = as.character(type),
    specificity_profile = as.character(specificity_profile),
    gene_name = as.character(gene_name)
  )

# =========================
# FORMAT CELL TYPE NAMES
# =========================

df <- df %>%
  mutate(
    celltypist_labels = case_when(
      celltypist_labels == "epithelial cell of thymus" ~ "TEC",
      celltypist_labels == "plasmacytoid dendritic cell, human" ~ "pDCs",
      celltypist_labels == "hematopoietic multipotent progenitor cell" ~ "MPPs",
      str_detect(str_to_lower(celltypist_labels), "cd4") ~
        str_to_upper(celltypist_labels),
      str_detect(str_to_lower(celltypist_labels), "cd8") ~
        str_to_upper(celltypist_labels),
      TRUE ~ tools::toTitleCase(celltypist_labels)
    ),
    
    type = case_when(
      str_to_lower(str_trim(type)) == "cancer-testis" ~ "Cancer-Testis",
      str_to_lower(str_trim(type)) == "diferentiation" ~ "Differentiation",
      TRUE ~ tools::toTitleCase(type)
    )
  )

# =========================
# MERGE HETEROGENEOUS GROUPS
# =========================

df <- df %>%
  mutate(
    specificity_profile = case_when(
      specificity_profile %in% c(
        "Enriched",
        "Heterogeneous Enriched"
      ) ~ "Enriched",
      
      specificity_profile %in% c(
        "Enhanced",
        "Heterogeneous Enhanced"
      ) ~ "Enhanced",
      
      specificity_profile %in% c(
        "Group enriched",
        "Group Enriched",
        "Heterogeneous Group Enriched"
      ) ~ "Group Enriched",
      
      TRUE ~ specificity_profile
    )
  )

# =========================
# FILTER ENRICHED CATEGORIES
# =========================

df_merged <- df %>%
  filter(
    specificity_profile %in% c(
      "Enriched",
      "Enhanced",
      "Group Enriched"
    )
  ) %>%
  filter(
    !str_detect(str_to_lower(type), "hla")
  )

if (nrow(df_merged) == 0) {
  stop(
    "No rows remain after filtering specificity profiles and HLA categories."
  )
}

# =========================
# OUTPUT DIRECTORIES
# =========================

output_files <- c(
  opt$output,
  opt$merged_output,
  opt$venn_output,
  opt$barplot_output
)

for (file in output_files) {
  dir.create(
    dirname(file),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# =========================
# SAVE MERGED RAW DATA
# =========================

write_csv(
  df_merged,
  opt$merged_output
)

# =========================
# HEATMAP SUMMARY
# =========================

heatmap_df <- df_merged %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  group_by(type, celltypist_labels) %>%
  summarise(
    total = n_distinct(gene_name),
    .groups = "drop"
  ) %>%
  complete(
    type,
    celltypist_labels,
    fill = list(total = 0)
  ) %>%
  mutate(
    type = fct_rev(factor(type))
  )

if (nrow(heatmap_df) == 0) {
  stop("No genes available to construct the heatmap.")
}

# =========================
# DEFINE TEXT COLOR
# =========================

dark_threshold <- quantile(
  heatmap_df$total,
  probs = 0.80,
  na.rm = TRUE
)

heatmap_df <- heatmap_df %>%
  mutate(
    text_color = if_else(
      total >= dark_threshold & total > 0,
      "white",
      "black"
    )
  )

# =========================
# PLOT HEATMAP
# =========================

p <- ggplot(
  heatmap_df,
  aes(
    x = celltypist_labels,
    y = type,
    fill = total
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.3
  ) +
  geom_text(
    aes(
      label = ifelse(total > 0, total, ""),
      color = text_color
    ),
    size = 3.5,
    family = "Arial"
  ) +
  scale_color_identity() +
  scale_fill_gradientn(
    colours = c(
      "#F2F2F2",
      "#D9D9D9",
      "#BFBFBF",
      "#A3A3A3",
      "#858585",
      "#666666",
      "#454545"
    )
  ) +
  theme_minimal(base_size = 14, base_family = "Arial") +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 1,
      size = 14,
      colour = "black",
      family = "Arial"
    ),
    axis.text.y = element_text(
      size = 14,
      colour = "black",
      family = "Arial"
    ),
    axis.title.x = element_text(
      colour = "black",
      family = "Arial",
      size = 14
    ),
    axis.title.y = element_text(
      colour = "black",
      family = "Arial",
      size = 14
    ),
    legend.title = element_text(
      size = 14,
      colour = "black",
      family = "Arial"
    ),
    legend.text = element_text(
      size = 14,
      colour = "black",
      family = "Arial"
    ),
    panel.grid = element_blank(),
    plot.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "transparent",
      colour = NA
    )
  ) +
  labs(
    x = NULL,
    y = "Antigen category",
    fill = "Enriched gene categories*"
  )

ggsave(
  filename = opt$output,
  plot = p,
  device = "svg",
  width = 13,
  height = 4,
  units = "in",
  bg = "transparent"
)

# ============================================================
# CTA DATA
# ============================================================

cta_df <- df_merged %>%
  filter(
    str_to_upper(str_trim(type)) == "CANCER-TESTIS",
    !is.na(gene_name),
    gene_name != "",
    !is.na(celltypist_labels),
    celltypist_labels != ""
  ) %>%
  distinct(gene_name, celltypist_labels)

if (nrow(cta_df) == 0) {
  stop("No CTA genes available for Venn diagram and barplot.")
}

# ============================================================
# 1. CTA VENN DIAGRAM
# ============================================================

malignant_genes <- cta_df %>%
  filter(
    str_detect(
      celltypist_labels,
      regex(opt$malignant_pattern, ignore_case = TRUE)
    )
  ) %>%
  pull(gene_name) %>%
  unique()

other_genes <- cta_df %>%
  filter(
    !str_detect(
      celltypist_labels,
      regex(opt$malignant_pattern, ignore_case = TRUE)
    )
  ) %>%
  pull(gene_name) %>%
  unique()

intersection_genes <- intersect(
  malignant_genes,
  other_genes
)

malignant_only <- setdiff(
  malignant_genes,
  other_genes
)

others_only <- setdiff(
  other_genes,
  malignant_genes
)

n_malignant_only <- length(malignant_only)
n_others_only <- length(others_only)
n_intersection <- length(intersection_genes)

# =========================
# VENN COLORS
# =========================

malignant_color <- "#9F3B4A"
others_color <- "#9E9E9E"
intersection_color <- "#9F6D75"

# =========================
# CIRCLE COORDINATES
# =========================

theta <- seq(
  0,
  2 * pi,
  length.out = 500
)

radius <- 0.34
circle_df <- bind_rows(
  data.frame(
    x = 0.39 + radius * cos(theta),
    y = 0.50 + radius * sin(theta),
    group = "Malignant"
  ),
  data.frame(
    x = 0.61 + radius * cos(theta),
    y = 0.50 + radius * sin(theta),
    group = "Non-Malignant"
  )
)

# =========================
# VENN PLOT
# =========================

venn_plot <- ggplot() +
  
  # Malignant circle
  geom_polygon(
    data = circle_df %>% filter(group == "Malignant"),
    aes(x = x, y = y),
    fill = "#9F3B4A",
    alpha = 0.65,
    colour = NA
  ) +
  
  # Others circle
  geom_polygon(
    data = circle_df %>% filter(group == "Non-Malignant"),
    aes(x = x, y = y),
    fill = "#9E9E9E",
    alpha = 0.65,
    colour = NA
  ) +
  
  # Malignant only
  annotate(
    "text",
    x = 0.15,
    y = 0.50,
    label = n_malignant_only,
    size = 4,
    family = "Arial",
    colour = "white"
  ) +
  # Intersection
  annotate(
    "text",
    x = 0.50,
    y = 0.50,
    label = n_intersection,
    size = 4,
    family = "Arial",
    colour = "black"
  ) +
  # Others only
  annotate(
    "text",
    x = 0.85,
    y = 0.50,
    label = n_others_only,
    size = 4,
    family = "Arial",
    colour = "black"
  ) +
  # Labels
  annotate(
    "text",
    x = 0.35,
    y = 0.87,
    label = "Malignant",
    size = 4,
    family = "Arial",
    colour = "black"
  ) +
  annotate(
    "text",
    x = 0.65,
    y = 0.87,
    label = "Non-Malignant",
    size = 4,
    family = "Arial",
    colour = "black"
  ) +
  coord_fixed(
    xlim = c(0.02, 0.98),
    ylim = c(0.10, 0.90),
    clip = "off"
  ) +
  theme_void(
    base_size = 12,
    base_family = "Arial"
  ) +
  theme(
    text = element_text(
      family = "Arial",
      size = 12
    ),
    plot.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "transparent",
      colour = NA
    )
  )

ggsave(
  filename = opt$venn_output,
  plot = venn_plot,
  device = "svg",
  width = 9,
  height = 3.5,
  units = "in",
  bg = "transparent"
)

# ============================================================
# 2. CTA ENRICHED CATEGORIES BARPLOT
# ============================================================

cta_barplot_df <- cta_df %>%
  group_by(celltypist_labels) %>%
  summarise(
    total = n_distinct(gene_name),
    .groups = "drop"
  ) %>%
  arrange(total) %>%
  mutate(
    celltypist_labels = factor(
      celltypist_labels,
      levels = celltypist_labels
    )
  )

# ============================================================
# BARPLOT — FIGURE 1C DESIGN
# ============================================================

barplot <- ggplot(
  cta_barplot_df,
  aes(
    x = reorder(celltypist_labels, total),
    y = total
  )
) +
  
  # Background reference bars
  geom_col(
    aes(y = max(total)),
    fill = "#f2f2f2",
    width = 0.72
  ) +
  
  # Main bars
  geom_col(
    fill = "#b7b7b7",
    width = 0.8
  ) +
  
  # Numerical labels
  geom_text(
    aes(label = scales::comma(total)),
    hjust = -0.09,
    size = 3,
    color = "black",
    family = "Arial"
  ) +
  
  coord_flip() +
  
  scale_y_continuous(
    labels = NULL,
    breaks = NULL,
    expand = expansion(
      mult = c(0, 0.25)
    )
  ) +
  
  labs(
    x = "Cell type",
    y = "Number of enriched CTA genes"
  ) +
  
  theme_classic(base_size = 8) +
  
  theme(
    text = element_text(
      family = "Arial",
      color = "black"
    ),
    
    axis.title.x = element_text(
      family = "Arial",
      color = "black",
      size = 7
    ),    
    axis.title.y = element_text(
      family = "Arial",
      color = "black",
      size = 7
    ),    
    axis.text.x = element_text(
      family = "Arial",
      color = "black",
      size = 7
    ),    
    axis.text.y = element_text(
      family = "Arial",
      color = "black",
      size = 7
    ),    
    axis.ticks.y = element_blank(),
    axis.line = element_line(
      linewidth = 0.5,
      color = "black"
    ),
    legend.position = "none",
    panel.background = element_rect(
      fill = "transparent",
      colour = NA
    ),    
    plot.background = element_rect(
      fill = "transparent",
      colour = NA
    ),    
    legend.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    legend.box.background = element_rect(
      fill = "transparent",
      colour = NA
    )
  )

ggsave(
  filename = opt$barplot_output,
  plot = barplot,
  device = "svg",
  width = 2,
  height = 3.5,
  units = "in",
  bg = "transparent"
)
