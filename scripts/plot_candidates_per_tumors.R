#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(ComplexUpset)
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  parsed <- list()
  i <- 1
  
  while (i <= length(args)) {
    key <- args[i]
    
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    
    if (i == length(args)) {
      stop("Missing value for argument: ", key)
    }
    
    parsed[[substring(key, 3)]] <- args[i + 1]
    i <- i + 2
  }
  
  parsed
}

params <- parse_args(args)

if (is.null(params$input)) {
  stop("Required argument: --input")
}

if (is.null(params$output_dir)) {
  stop("Required argument: --output_dir")
}

if (is.null(params$epitopes)) {
  stop("Required argument: --epitopes")
}

input_file <- params$input
epitopes_file <- params$epitopes
output_dir <- params$output_dir

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

if (!file.exists(epitopes_file)) {
  stop("Epitope file does not exist: ", epitopes_file)
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ==========================================================
# COLORS
# ==========================================================

COLOR_EPITOPE <- "#89988D"
COLOR_NO_EPITOPE <- "#A6A6A6"

EPITOPE_COLORS <- c(
  "With epitope" = COLOR_EPITOPE,
  "Without epitope" = COLOR_NO_EPITOPE
)

# ==========================================================
# LOAD CANDIDATE MATRIX
# ==========================================================

df <- read_csv(
  input_file,
  show_col_types = FALSE
)

required_columns <- c("Gene_ID", "Gene_Name")

missing_columns <- setdiff(
  required_columns,
  colnames(df)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ==========================================================
# LOAD EPITOPE DATA
# ==========================================================

epitopes_df <- read_csv(
  epitopes_file,
  show_col_types = FALSE
)

if (!"gene_name" %in% colnames(epitopes_df)) {
  stop("Missing required column in epitope file: gene_name")
}

epitope_genes <- epitopes_df %>%
  transmute(
    Gene_Name = trimws(as.character(gene_name))
  ) %>%
  filter(
    !is.na(Gene_Name),
    Gene_Name != ""
  ) %>%
  distinct() %>%
  pull(Gene_Name)


# ==========================================================
# ANNOTATE CANDIDATE GENES
# ==========================================================

df <- df %>%
  mutate(
    Gene_Name = trimws(as.character(Gene_Name)),
    
    Has_epitope = !is.na(Gene_Name) &
      Gene_Name != "" &
      Gene_Name %in% epitope_genes,
    
    Epitope_status = factor(
      ifelse(
        Has_epitope,
        "With epitope",
        "Without epitope"
      ),
      levels = c(
        "With epitope",
        "Without epitope"
      )
    )
  )

# ==========================================================
# DEFINE TUMOR SETS
# ==========================================================

sets <- setdiff(
  colnames(df),
  c(
    "Gene_ID",
    "Gene_Name",
    "Has_epitope",
    "Epitope_status"
  )
)

if (length(sets) == 0) {
  stop("No tumor columns were found in the input matrix.")
}
sets_title <- tools::toTitleCase(tolower(sets))

if (anyDuplicated(sets_title)) {
  stop(
    "Tumor names become duplicated after capitalization: ",
    paste(sets_title, collapse = ", ")
  )
}
colnames(df)[match(sets, colnames(df))] <- sets_title
sets <- sets_title

for (tumor in sets) {
  df[[tumor]] <- as.numeric(df[[tumor]])
  
  if (any(
    !is.na(df[[tumor]]) &
    !df[[tumor]] %in% c(0, 1)
  )) {
    stop(
      "Column '", tumor,
      "' contains values other than 0 and 1."
    )
  }
  
  df[[tumor]][is.na(df[[tumor]])] <- 0
}

# ==========================================================
# CREATE GENE LIST BY TUMOR COMBINATION
# ==========================================================

df$Combination <- apply(
  as.matrix(df[, sets, drop = FALSE]),
  1,
  function(x) {
    paste(sets[x == 1], collapse = "_")
  }
)

df_combination <- df %>%
  filter(Combination != "")

gene_combination_table <- df_combination %>%
  group_by(Combination) %>%
  summarise(
    N_genes = n(),
    N_tumors = n_distinct(
      unlist(
        strsplit(
          first(Combination),
          "_",
          fixed = TRUE
        )
      )
    ),
    N_genes_with_epitopes = sum(Has_epitope),
    N_genes_without_epitopes = sum(!Has_epitope),
    Genes = paste(
      ifelse(
        is.na(Gene_Name) | Gene_Name == "",
        Gene_ID,
        Gene_Name
      ),
      collapse = "; "
    ),
    Genes_with_epitopes = paste(
      ifelse(
        Has_epitope,
        ifelse(
          is.na(Gene_Name) | Gene_Name == "",
          Gene_ID,
          Gene_Name
        ),
        NA_character_
      ) %>%
        na.omit(),
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(N_genes), desc(N_tumors))

summary_file <- file.path(
  output_dir,
  "genes_by_tumor_combination_summary.csv"
)

write_csv(
  gene_combination_table,
  summary_file
)

# ==========================================================
# PREPARE UPSET DATA
# ==========================================================

plot_file <- file.path(
  output_dir,
  "upset_candidate_genes.svg"
)
df_plot <- df
for (tumor in sets) {
  df_plot[[tumor]] <- df_plot[[tumor]] == 1
}

# ==========================================================
# CREATE UPSET PLOT
# ==========================================================

p <- upset(
  df_plot,
  intersect = sets,
  min_size = 3,
  sort_intersections_by = "cardinality",
  width_ratio = 0.22,
  
  base_annotations = list(
    "Intersection size" =
      intersection_size(
        counts = TRUE,
        mapping = aes(fill = Epitope_status)
      ) +
      scale_fill_manual(
        name = "Number of Genes by Epitope Status",
        values = EPITOPE_COLORS,
        breaks = c(
          "With epitope",
          "Without epitope"
        ),
        labels = c(
          "With epitope" = "Found",
          "Without epitope" = "Not Found"
        ),
        drop = FALSE
      ) +
      labs(
        x = NULL,
        y = "Number of Genes"
      ) +
      theme(
        legend.position = "top",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )
  ),
  matrix = intersection_matrix() +
    labs(x = "Group", y = "Cancers") +
    theme(
      axis.title.y = element_text(
        family = "Arial",
        size = 11,
        colour = "black"
      ),
      strip.text = element_blank(),
      strip.background = element_blank()
    ),
  set_sizes =
    upset_set_size() +
    labs(y = "Set Size") +
    theme(
      legend.position = "none"
    )
) +
  theme_minimal(
    base_size = 14,
    base_family = "Arial"
  ) +
  theme(
    panel.grid.major = element_line(
      colour = "#F5F6F7",
      linewidth = 0.3
    ),
    
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    text = element_text(
      family = "Arial",
      colour = "black",
      size = 11
    ),
    
    axis.text = element_text(
      family = "Arial",
      size = 11,
      colour = "black"
    ),
    
    axis.title = element_text(
      family = "Arial",
      size = 11,
      colour = "black"
    ),
    
    strip.text = element_text(
      family = "Arial",
      size = 11,
      colour = "black"
    ),
    
    plot.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    
    panel.background = element_rect(
      fill = "#F5F6F7",
      colour = NA
    )
  )

ggsave(
  filename = plot_file,
  plot = p,
  device = "svg",
  width = 10.7,
  height = 5.7,
  bg = "transparent"
)