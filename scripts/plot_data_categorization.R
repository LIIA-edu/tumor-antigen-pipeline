#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(svglite)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(stringr)
  library(optparse)
})


option_list <- list(
  make_option(
    "--input-h5ad",
    type = "character",
    help = "Input H5AD file"
  ),
  make_option(
    "--output-dir",
    type = "character",
    default = ".",
    help = "Output directory [default: .]"
  )
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

if (is.null(args$`input-h5ad`)) {
  stop("Required argument: --input-h5ad")
}

h5ad_file <- args$`input-h5ad`
output_dir <- args$`output-dir`

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(h5ad_file)) {
  stop("Input H5AD file does not exist: ", h5ad_file)
}

# =========================
# LOAD H5AD
# =========================

sce <- readH5AD(h5ad_file)

# =========================
# COLUMN CONFIGURATION
# =========================

celltype_col <- "cell_type_grouped"
tumor_col <- "disease"
malignant_col <- "malignant cell"

required_columns <- c(
  celltype_col,
  tumor_col
)

if (!all(required_columns %in% colnames(colData(sce)))) {
  stop(
    "Required columns not found in H5AD colData: ",
    paste(
      setdiff(required_columns, colnames(colData(sce))),
      collapse = ", "
    )
  )
}

if (!"X_umap" %in% reducedDimNames(sce)) {
  stop("X_umap was not found in H5AD reduced dimensions.")
}

# =========================
# EXTRACT COUNTS DIRECTLY
# FROM H5AD
# =========================


obs_df <- data.frame(
  cell_type_grouped = as.character(
    colData(sce)[[celltype_col]]
  ),
  disease = as.character(
    colData(sce)[[tumor_col]]
  ),
  stringsAsFactors = FALSE
)

# Figure 1C: Cells per tumor type

tumor_df <- obs_df %>%
  filter(
    !is.na(disease),
    disease != ""
  ) %>%
  count(
    disease,
    name = "n",
    sort = FALSE
  )

# Figure 1D: Cell-type composition

df_fig1d <- obs_df %>%
  filter(
    !is.na(cell_type_grouped),
    cell_type_grouped != ""
  ) %>%
  count(
    cell_type_grouped,
    name = "n",
    sort = FALSE
  )

save_plot <- function(
    plot,
    filename,
    width,
    height,
    dpi = 300,
    svg = TRUE,
    png = TRUE
) {
  
  if (svg) {
    ggsave(
      filename = file.path(
        output_dir,
        paste0(filename, ".svg")
      ),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      device = svglite::svglite,
      bg = "transparent"
    )
  }
  
  if (png) {
    ggsave(
      filename = file.path(
        output_dir,
        paste0(filename, ".png")
      ),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "transparent"
    )
  }
}

# ============================================================
# FIGURE 1C
# Cells per tumor type
# ============================================================

tumor_df <- tumor_df %>%
  mutate(
    disease = str_to_title(disease)
  )

fig1C <- ggplot(
  tumor_df,
  aes(
    x = reorder(disease, n),
    y = n
  )
) +
  
  geom_col(
    aes(y = max(n)),
    fill = "#f2f2f2",
    width = 0.72
  ) +
  
  geom_col(
    fill = "#b7b7b7",
    width = 0.8
  ) +
  
  geom_text(
    aes(label = scales::comma(n)),
    hjust = -0.09,
    size = 3.0,
    color = "black",
    family = "Arial"
  ) +
  
  coord_flip() +
  
  scale_y_continuous(
    labels = NULL,
    breaks = NULL,
    expand = expansion(mult = c(0, 0.25))
  ) +
  
  labs(
    x = "Tumor type",
    y = "Number of cells"
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
      size = 8
    ),
    axis.title.y = element_text(
      family = "Arial",
      color = "black",
      size = 9
    ),
    axis.text.x = element_text(
      family = "Arial",
      color = "black",
      size = 9
    ),
    axis.text.y = element_text(
      family = "Arial",
      color = "black",
      size = 9
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

save_plot(
  fig1C,
  "Figure1C_CellsPerTumorType",
  width = 6,
  height = 2,
  dpi = 300
)

# ============================================================
# FIGURE 1D
# Cell-type composition
# ============================================================

celltype_df <- df_fig1d %>%
  mutate(
    cell_type_grouped = str_to_title(cell_type_grouped)
  ) %>%
  mutate(
    celltype_plot = case_when(
      cell_type_grouped == "Epithelial Cell Of Thymus" ~ "TEC",
      cell_type_grouped == "Plasmacytoid Dendritic Cell, Human" ~ "pDCs",
      cell_type_grouped == "Hematopoietic Multipotent Progenitor Cell" ~ "MPPs",
      TRUE ~ cell_type_grouped
    )
  )

fig1D <- ggplot(
  celltype_df,
  aes(
    x = reorder(celltype_plot, n),
    y = n
  )
) +
  
  geom_col(
    aes(y = max(n)),
    fill = "#f2f2f2",
    width = 0.8
  ) +
  
  geom_col(
    fill = "#b7b7b7",
    width = 0.8
  ) +
  
  geom_text(
    aes(label = scales::comma(n)),
    hjust = -0.09,
    size = 3,
    color = "black",
    family = "Arial"
  ) +
  
  coord_flip() +
  
  scale_y_continuous(
    labels = NULL,
    breaks = NULL,
    expand = expansion(mult = c(0, 0.25))
  ) +
  
  labs(
    x = "Cell type",
    y = "Number of cells"
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
      size = 3
    ),
    axis.title.y = element_text(
      family = "Arial",
      color = "black",
      size = 3
    ),
    axis.text.x = element_text(
      family = "Arial",
      color = "black",
      size = 3
    ),
    axis.text.y = element_text(
      family = "Arial",
      color = "black",
      size = 3
    ),
    axis.ticks.y = element_blank(),
    axis.line = element_line(
      linewidth = 0.5,
      color = "black"
    ),
    legend.position = "none"
  )

ggsave(
  filename = file.path(
    output_dir,
    "Figure1D_CellTypeComposition.svg"
  ),
  plot = fig1D,
  width = 6,
  height = 5,
  units = "in",
  device = svglite::svglite,
  bg = "transparent"
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure1D_CellTypeComposition.png"
  ),
  plot = fig1D,
  width = 4,
  height = 5,
  dpi = 300,
  bg = "white"
)

# ============================================================
# FIGURE 1E
# UMAP malignant vs non-malignant
# ============================================================

umap <- reducedDim(sce, "X_umap")

umap_df <- data.frame(
  UMAP1 = umap[, 1],
  UMAP2 = umap[, 2],
  Status = ifelse(
    colData(sce)[[celltype_col]] == malignant_col,
    "Malignant",
    "Non-Malignant"
  )
)

# Plot non-malignant first
umap_df <- umap_df %>%
  arrange(Status == "Malignant")

fig1E <- ggplot(
  umap_df,
  aes(
    UMAP1,
    UMAP2,
    color = Status
  )
) +
  
  geom_point(
    size = 0.1,
    alpha = 0.45
  ) +
  
  theme_classic() +
  
  labs(
    color = NULL,
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  
  scale_color_manual(
    values = c(
      "Malignant" = "#9F3B4A",
      "Non-Malignant" = "#9E9E9E"
    ),
    guide = guide_legend(
      override.aes = list(
        size = 2.5,
        alpha = 1
      )
    )
  ) +
  
  theme(
    axis.title = element_text(
      size = 12,
      color = "black",
      family = "Arial"
    ),
    axis.text = element_text(
      size = 12,
      color = "black",
      family = "Arial"
    ),
    legend.text = element_text(
      size = 12,
      family = "Arial"
    ),
    legend.position = "right"
  )

save_plot(
  fig1E,
  "Figure1E_MalignantVsNonMalignant_UMAP",
  width = 5,
  height = 5,
  dpi = 400
)

# ============================================================
# FIGURE 1F
# UMAP by disease
# ============================================================

disease_colors <- c(
  
  # Breast cancer
  "HER2 positive breast carcinoma" = "#F5C4D3",
  "triple-negative breast carcinoma" = "#D96B91",
  "estrogen-receptor positive breast cancer" = "#9B2452",
  
  # Skin cancer
  "basal cell carcinoma" = "#B7D7A8",
  "melanoma" = "#4F8A5B",
  "metastatic melanoma" = "#1F4D35",
  
  # Kidney cancer
  "nonpapillary renal cell carcinoma" = "#D9A62E",
  
  # Liver cancer
  "intrahepatic cholangiocarcinoma" = "#8EC5E8",
  "hepatocellular carcinoma" = "#145A8D",
  
  # Squamous cell carcinoma
  "squamous cell carcinoma" = "#F2D6A2"
)

umap_df <- data.frame(
  UMAP1 = umap[, 1],
  UMAP2 = umap[, 2],
  Disease = colData(sce)[[tumor_col]]
)

# Remove cells without disease information

umap_df <- umap_df %>%
  filter(!is.na(Disease)) %>%
  mutate(
    Disease_label = str_to_title(Disease),
    Disease_label = factor(
      Disease_label,
      levels = str_to_title(names(disease_colors))
    )
  ) %>%
  filter(!is.na(Disease_label))

# Fix legend order

umap_df <- umap_df %>%
  mutate(
    Disease_label = factor(
      Disease_label,
      levels = str_to_title(names(disease_colors))
    )
  )

fig1F <- ggplot(
  umap_df,
  aes(
    x = UMAP1,
    y = UMAP2,
    color = Disease_label
  )
) +
  
  geom_point(
    size = 0.1,
    alpha = 0.45
  ) +
  
  scale_color_manual(
    values = setNames(
      disease_colors,
      str_to_title(names(disease_colors))
    ),
    drop = FALSE,
    guide = guide_legend(
      ncol = 1,
      byrow = TRUE,
      override.aes = list(
        size = 2.5,
        alpha = 1
      ),
      keywidth = 0.3,
      keyheight = 0.3
    )
  ) +
  
  labs(
    color = NULL,
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  
  coord_fixed() +
  
  theme_classic() +
  
  theme(
    axis.title = element_text(
      size = 10,
      color = "black"
    ),
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    legend.text = element_text(
      size = 6,
      color = "black"
    ),
    legend.position = "right",
    legend.key = element_blank(),
    legend.key.height = unit(0.28, "cm"),
    legend.key.width = unit(0.28, "cm"),
    legend.spacing.y = unit(0.02, "cm"),
    plot.margin = margin(
      t = 3,
      r = 3,
      b = 3,
      l = 3
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "Figure1F_Disease_UMAP.png"
  ),
  plot = fig1F,
  width = 5,
  height = 5,
  units = "in",
  dpi = 400,
  bg = "transparent"
)