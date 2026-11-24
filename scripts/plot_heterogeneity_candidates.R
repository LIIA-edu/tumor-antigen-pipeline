
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  parsed <- list()
  i <- 1
  
  while (i <= length(args)) {
    if (startsWith(args[i], "--")) {
      key <- substring(args[i], 3)
      
      if (i == length(args) || startsWith(args[i + 1], "--")) {
        parsed[[key]] <- TRUE
        i <- i + 1
      } else {
        parsed[[key]] <- args[i + 1]
        i <- i + 2
      }
    } else {
      i <- i + 1
    }
  }
  
  parsed
}

parsed <- parse_args(args)

# ==========================================================
# OUTPUT DIRECTORY
# ==========================================================

OUTPUT_PATH <- parsed$output_dir %||% "."

dir.create(
  OUTPUT_PATH,
  recursive = TRUE,
  showWarnings = FALSE
)

# ==========================================================
# READ INPUT FILES
# ==========================================================

modules_file <- parsed$modules
mean_expression_file <- parsed[["mean-expression"]]
pct_expression_file <- parsed[["pct-expression"]]
epitopes_file <- parsed$epitopes

required_files <- list(
  modules = modules_file,
  mean_expression = mean_expression_file,
  pct_expression = pct_expression_file,
  epitopes = epitopes_file
)

for (file_name in names(required_files)) {
  file_path <- required_files[[file_name]]
  
  if (
    is.null(file_path) ||
    !is.character(file_path) ||
    length(file_path) != 1 ||
    !file.exists(file_path)
  ) {
    stop(
      "Input file missing or not found for '",
      file_name,
      "': ",
      if (is.null(file_path)) "<not provided>" else file_path
    )
  }
}

modules <- read.csv(
  modules_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

expr <- read.csv(
  mean_expression_file,
  row.names = 1,
  check.names = FALSE
)

pct <- read.csv(
  pct_expression_file,
  row.names = 1,
  check.names = FALSE
)

epitopes <- read.csv(
  epitopes_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"gene_name" %in% colnames(epitopes)) {
  stop("Column 'gene_name' not found in the epitope file.")
}

epitope_genes <- epitopes %>%
  transmute(gene_symbol = trimws(as.character(gene_name))) %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  distinct() %>%
  pull(gene_symbol)

if (length(epitope_genes) == 0) {
  stop("No valid gene names found in the epitope file column 'gene_name'.")
}

# ==========================================================
# CLEAN GENE IDENTIFIERS
# ==========================================================

if ("gene" %in% colnames(modules) &&
    !"gene_symbol" %in% colnames(modules)) {
  modules <- modules %>%
    rename(gene_symbol = gene)
}

if (!"gene_symbol" %in% colnames(modules)) {
  stop(
    "Gene column not found. Expected 'gene' or 'gene_symbol'."
  )
}

if (!"module" %in% colnames(modules)) {
  stop("Column 'module' not found in the module table.")
}

modules$gene_symbol <- trimws(
  as.character(modules$gene_symbol)
)

modules$gene_symbol <- gsub(
  " \\(ENSG.*\\)",
  "",
  modules$gene_symbol
)

rownames(expr) <- gsub(
  " \\(ENSG.*\\)",
  "",
  rownames(expr)
)

rownames(pct) <- gsub(
  " \\(ENSG.*\\)",
  "",
  rownames(pct)
)

modules <- modules %>%
  filter(
    !is.na(gene_symbol),
    gene_symbol != "",
    !is.na(module)
  ) %>%
  distinct(gene_symbol, .keep_all = TRUE)

# Retain only genes present in the epitope file (gene_name column)
n_modules_before_epitope_filter <- n_distinct(modules$gene_symbol)

modules <- modules %>%
  filter(gene_symbol %in% epitope_genes)

n_modules_after_epitope_filter <- n_distinct(modules$gene_symbol)

if (nrow(modules) == 0) {
  stop("No genes in the module table are present in the epitope file column 'gene_name'.")
}

# ==========================================================
# MATCH GENES AND TUMOR TYPES
# ==========================================================

genes <- Reduce(
  intersect,
  list(
    modules$gene_symbol,
    rownames(expr),
    rownames(pct)
  )
)

if (length(genes) == 0) {
  stop(
    "No overlapping genes found among the three input files."
  )
}

modules <- modules %>%
  filter(gene_symbol %in% genes) %>%
  arrange(module, gene_symbol)

gene_order <- modules$gene_symbol

expr <- expr[gene_order, , drop = FALSE]
pct <- pct[gene_order, , drop = FALSE]

common_diseases <- intersect(
  colnames(expr),
  colnames(pct)
)

if (length(common_diseases) == 0) {
  stop("No matching tumor-type columns found.")
}

expr <- expr[, common_diseases, drop = FALSE]
pct <- pct[, common_diseases, drop = FALSE]


# ==========================================================
# CONVERT MATRICES TO LONG FORMAT
# ==========================================================

expr_long <- expr %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("gene_symbol") %>%
  pivot_longer(
    cols = -gene_symbol,
    names_to = "disease",
    values_to = "mean_expression"
  )

pct_long <- pct %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("gene_symbol") %>%
  pivot_longer(
    cols = -gene_symbol,
    names_to = "disease",
    values_to = "pct_expression"
  )

modules_plot <- modules %>%
  select(gene_symbol, module) %>%
  distinct(gene_symbol, .keep_all = TRUE)

plot_df <- expr_long %>%
  inner_join(
    pct_long,
    by = c("gene_symbol", "disease")
  ) %>%
  inner_join(
    modules_plot,
    by = "gene_symbol"
  )

# ==========================================================
# VALIDATE DATA
# ==========================================================

required_columns <- c(
  "gene_symbol",
  "disease",
  "mean_expression",
  "pct_expression",
  "module"
)

missing_columns <- setdiff(
  required_columns,
  colnames(plot_df)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing columns after joins: ",
    paste(missing_columns, collapse = ", ")
  )
}

plot_df$mean_expression <- as.numeric(
  plot_df$mean_expression
)

plot_df$pct_expression <- as.numeric(
  plot_df$pct_expression
)

plot_df <- plot_df %>%
  filter(
    is.finite(mean_expression),
    is.finite(pct_expression),
    !is.na(module)
  )

if (nrow(plot_df) == 0) {
  stop("No valid observations remain for plotting.")
}

# ==========================================================
# NORMALIZE PERCENTAGES
# ==========================================================

max_pct <- max(
  plot_df$pct_expression,
  na.rm = TRUE
)

if (max_pct <= 1) {
  plot_df$pct_expression <- plot_df$pct_expression * 100
}

plot_df$pct_expression <- pmin(
  pmax(plot_df$pct_expression, 0),
  100
)

# ==========================================================
# REMOVE ZERO-EXPRESSION OBSERVATIONS
# No dot will be drawn when mean expression is zero.
# ==========================================================

n_before_filter <- nrow(plot_df)

plot_df <- plot_df %>%
  filter(mean_expression > 0)

if (nrow(plot_df) == 0) {
  stop(
    "All observations have zero mean expression. ",
    "No dotplots can be generated."
  )
}

# ==========================================================
# ORDERING
# ==========================================================

module_levels <- unique(modules$module)

plot_df$module <- factor(
  plot_df$module,
  levels = module_levels
)

plot_df$gene_symbol <- as.character(
  plot_df$gene_symbol
)

plot_df$disease <- factor(
  plot_df$disease,
  levels = common_diseases
)

# ==========================================================
# SHARED THEME AND SCALES
# ==========================================================

dotplot_theme <- theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(
      colour = "grey70",
      fill = NA,
      linewidth = 0.2
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

dotplot_scales <- list(
  scale_size_continuous(
    name = "% cells expressing",
    range = c(0, 8),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100)
  ),
  scale_colour_gradientn(
    colours = c(
      "#B69FAF",
      "#9F3B4A",
      "#3A0F1A"
    ),
    name = "Mean expression"
  )
)

# ==========================================================
# 1. DOTPLOT WITH MODULE FACETS
# ==========================================================

plot_modules_df <- plot_df %>%
  mutate(
    gene_symbol = factor(
      gene_symbol,
      levels = rev(gene_order)
    )
  )

p_modules <- ggplot(
  plot_modules_df,
  aes(
    x = disease,
    y = gene_symbol
  )
) +
  geom_point(
    aes(
      size = pct_expression,
      colour = mean_expression
    )
  ) +
  facet_grid(
    module ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  dotplot_scales +
  labs(
    x = "Tumor type",
    y = "Candidate antigen"
  ) +
  dotplot_theme +
  theme(
    strip.text.y = element_text(
      angle = 0,
      face = "bold"
    )
  )

modules_width <- 14
modules_height <- max(10, 0.18 * length(gene_order) + 3)

ggsave(
  filename = file.path(OUTPUT_PATH, "DotPlot_gene_modules_by_disease.svg"),
  plot = p_modules,
  width = 14,
  height = max(10, 0.18 * length(gene_order) + 3),
  device = "svg",
  limitsize = FALSE
)

ggsave(
  filename = file.path(OUTPUT_PATH, "DotPlot_gene_modules_by_disease.pdf"),
  plot = p_modules,
  width = 14,
  height = max(10, 0.18 * length(gene_order) + 3),
  device = cairo_pdf,
  limitsize = FALSE
)


# ==========================================================
# 2. DOTPLOT WITHOUT MODULE FACETS
# ==========================================================

gene_order_unfaceted <- modules %>%
  arrange(module, gene_symbol) %>%
  pull(gene_symbol) %>%
  unique()

plot_all_df <- plot_df %>%
  mutate(
    gene_symbol = factor(
      gene_symbol,
      levels = rev(gene_order_unfaceted)
    )
  )

p_all <- ggplot(
  plot_all_df,
  aes(
    x = disease,
    y = gene_symbol
  )
) +
  geom_point(
    aes(
      size = pct_expression,
      colour = mean_expression
    )
  ) +
  dotplot_scales +
  labs(
    x = "Tumor type",
    y = "Candidate antigen"
  ) +
  dotplot_theme

all_genes_width <- 14
all_genes_height <- max(10, 0.18 * length(gene_order_unfaceted) + 3)

all_genes_pdf <- file.path(OUTPUT_PATH, "DotPlot_candidates_by_disease.pdf")
all_genes_svg <- file.path(OUTPUT_PATH, "DotPlot_candidates_by_disease.svg")

ggsave(
  filename = all_genes_pdf,
  plot = p_all,
  width = all_genes_width,
  height = all_genes_height,
  units = "in",
  limitsize = FALSE
)

ggsave(
  filename = all_genes_svg,
  plot = p_all,
  device = "svg",
  width = all_genes_width,
  height = all_genes_height,
  units = "in",
  limitsize = FALSE
)
