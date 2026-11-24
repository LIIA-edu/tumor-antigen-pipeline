
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste(
      "Usage:",
      "Rscript epitope_plots.R",
      "<epitopes.csv>",
      "<malignant_candidates.csv>",
      "<output_dir>"
    )
  )
}

epitopes_file <- args[1]
candidates_file <- args[2]
output_dir <- args[3]

# =========================
# VALIDATE INPUTS
# =========================

if (!file.exists(epitopes_file)) {
  stop("Epitope file not found: ", epitopes_file)
}

if (!file.exists(candidates_file)) {
  stop("Candidate file not found: ", candidates_file)
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

epi <- read_csv(
  epitopes_file,
  show_col_types = FALSE
)

candidates <- read_csv(
  candidates_file,
  show_col_types = FALSE
)

# =========================
# VALIDATE COLUMNS
# =========================

required_epi_columns <- c(
  "gene_name",
  "cancer",
  "epitope_id_clean"
)

missing_epi <- setdiff(
  required_epi_columns,
  colnames(epi)
)

if (length(missing_epi) > 0) {
  stop(
    "Missing columns in epitope file: ",
    paste(missing_epi, collapse = ", ")
  )
}

if (!"gene_name" %in% colnames(candidates)) {
  stop(
    "Candidate file must contain a gene_name column."
  )
}

# =========================
# GENE COUNTS
# =========================

malignant_genes <- candidates %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  distinct(gene_name)

n_malignant <- nrow(malignant_genes)

epitope_genes <- epi %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  distinct(gene_name)

# Count only genes present in the candidate universe
n_epitope <- sum(
  malignant_genes$gene_name %in% epitope_genes$gene_name
)

n_not_found <- n_malignant - n_epitope

if (n_not_found < 0) {
  stop("Invalid gene counts. Check the input files.")
}

# =========================
# COLORS
# =========================

col_mid  <- "#A6A6A6"
col_dark <- "#89988D"
col_bg   <- "#F2F2F2"

# =========================
# COMMON THEME
# =========================

common_theme <- theme_classic(
  base_family = "Arial",
  base_size = 9
) +
  theme(
    text = element_text(
      family = "Arial",
      color = "black",
      size = 9
    ),
    
    axis.text = element_text(
      family = "Arial",
      color = "black",
      size = 9
    ),
    
    axis.title = element_text(
      family = "Arial",
      color = "black",
      size = 12
    ),
    
    plot.title = element_text(
      family = "Arial",
      color = "black",
      hjust = 0.5,
      face = "plain",
      size = 9,
      margin = margin(b = 12)
    ),
    
    plot.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    
    panel.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    
    legend.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    
    legend.key = element_rect(
      fill = "transparent",
      colour = NA
    ),
    
    plot.margin = margin(
      t = 12,
      r = 15,
      b = 12,
      l = 15
    )
  )

save_plot <- function(plot, filename, width, height) {
  
  ggsave(
    filename = file.path(
      output_dir,
      filename
    ),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = "svg",
    bg = "transparent"
  )
}

# =====================================================
# A. GENES WITH KNOWN EPITOPE ASSOCIATIONS
# =====================================================

plot_a_df <- data.frame(
  Category = factor(
    c(
      "Epitope-validated",
      "No epitope match"
    ),
    levels = c(
      "Epitope-validated",
      "No epitope match"
    )
  ),
  Count = c(
    n_epitope,
    n_not_found
  )
)

p1 <- ggplot(
  plot_a_df,
  aes(
    x = Count,
    y = "Genes",
    fill = Category
  )
) +
  
  geom_col(
    width = 0.45,
    position = "stack"
  ) +
  
  geom_text(
    aes(label = Count),
    position = position_stack(vjust = 0.5),
    size = 12 / .pt,
    family = "Arial",
    color = "black"
  ) +
  
  scale_fill_manual(
    values = c(
      "Epitope-validated" = col_dark,
      "No epitope match" = col_mid
    ),
    breaks = c(
      "Epitope-validated",
      "No epitope match"
    )
  ) +
  
  scale_x_continuous(
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  
  labs(
    title = "Genes with known epitope associations",
    x = "Number of genes",
    y = NULL,
    fill = NULL
  ) +
  
  common_theme +
  
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(
      family = "Arial",
      size = 12,
      color = "black"
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_line(
      linewidth = 0.4
    ),
    axis.title.x = element_text(
      margin = margin(t = 8)
    )
  )

save_plot(
  p1,
  "plot_A_genes_with_known_epitope_associations.svg",
  width = 6,
  height = 3
)

# =====================================================
# B. CANCER ASSOCIATIONS OF EPITOPE-VALIDATED GENES
# =====================================================

format_cancer_name <- function(x) {
  
  x %>%
    str_replace_all("_", " ") %>%
    str_squish() %>%
    str_to_sentence()
}

cancer_counts <- epi %>%
  filter(
    !is.na(gene_name),
    gene_name != "",
    !is.na(cancer),
    cancer != ""
  ) %>%
  distinct(
    gene_name,
    cancer
  ) %>%
  count(
    cancer,
    name = "n_genes"
  ) %>%
  mutate(
    cancer_label = format_cancer_name(cancer)
  ) %>%
  arrange(n_genes, cancer_label) %>%
  mutate(
    cancer_label = factor(
      cancer_label,
      levels = cancer_label
    ),
    background = max(n_genes)
  )

p2 <- ggplot(
  cancer_counts,
  aes(
    y = cancer_label
  )
) +
  
  # Uniform background bars
  geom_col(
    aes(x = background),
    fill = col_bg,
    width = 0.55
  ) +
  
  # Actual gene counts
  geom_col(
    aes(x = n_genes),
    fill = col_dark,
    width = 0.55
  ) +
  
  # Count labels
  geom_text(
    aes(
      x = n_genes,
      label = n_genes
    ),
    hjust = -0.45,
    size = 12 / .pt,
    family = "Arial",
    color = "black"
  ) +
  
  scale_x_continuous(
    limits = c(
      0,
      max(cancer_counts$n_genes) * 1.18
    ),
    breaks = scales::pretty_breaks(),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  
  scale_y_discrete(
    expand = expansion(
      add = c(0.02, 0.02)
    )
  ) +
  
  labs(
    title = "Cancer associations of epitope-validated genes",
    x = "Number of genes",
    y = NULL
  ) +
  
  common_theme +
  
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(
      linewidth = 0.4
    ),
    panel.grid = element_blank()
  )

save_plot(
  p2,
  "plot_B_cancer_associations.svg",
  width = 6,
  height = 4.5
)

# =====================================================
# GENES GROUPED BY EPITOPE COUNT
# =====================================================

epitope_counts <- epi %>%
  filter(
    !is.na(gene_name),
    gene_name != "",
    !is.na(epitope_id_clean),
    epitope_id_clean != ""
  ) %>%
  group_by(gene_name) %>%
  summarise(
    n_epitopes = n_distinct(epitope_id_clean),
    .groups = "drop"
  )

genes_by_epitope_count <- epitope_counts %>%
  group_by(n_epitopes) %>%
  summarise(
    genes = paste(
      sort(unique(gene_name)),
      collapse = ", "
    ),
    n_genes = n(),
    .groups = "drop"
  ) %>%
  arrange(n_epitopes)

write.csv(
  genes_by_epitope_count,
  file.path(output_dir, "genes_by_epitope_count.csv"),
  row.names = FALSE
)