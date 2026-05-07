library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Load data
df <- read_tsv("vOTUs_metadata.tsv",
               show_col_types = FALSE)

# Build vOTU x sample RPKM matrix (take max RPKM per vOTU across duplicate entries)
rpkm_cols <- names(df)[startsWith(names(df), "RPKM_")]

mat_full <- df %>%
  select(virus_id, all_of(rpkm_cols)) %>%
  group_by(virus_id) %>%
  summarise(across(everything(), max), .groups = "drop") %>%
  column_to_rownames("virus_id")

colnames(mat_full) <- gsub("^RPKM_", "", colnames(mat_full))

# Column order: Hexactinellida (Aphro18, Farrea) first, then Demospongiae alphabetically
# Povisternata is placed before Phycopsis
hex_cols  <- intersect(c("Aphro18", "Farrea"), colnames(mat_full))
demo_cols <- sort(setdiff(colnames(mat_full), hex_cols))

ph_i <- which(demo_cols == "Phycopsis")
po_i <- which(demo_cols == "Povisternata")
if (length(ph_i) && length(po_i)) demo_cols[c(ph_i, po_i)] <- demo_cols[c(po_i, ph_i)]

col_order <- c(hex_cols, demo_cols)
mat_full  <- mat_full[, col_order]

# log10(x + 1) transform
log10_na <- function(mat) log10(as.matrix(mat) + 1)

# Metadata lookups
quality_lookup <- df %>%
  select(virus_id, checkv_quality) %>%
  distinct(virus_id, .keep_all = TRUE)

# Use the lowest taxonomic rank from the semicolon-delimited string
tax_last <- function(x) {
  sapply(strsplit(as.character(x), ";"),
         function(v) trimws(tail(v[nchar(trimws(v)) > 0], 1)),
         USE.NAMES = FALSE)
}

taxonomy_lookup <- df %>%
  select(virus_id, taxonomy) %>%
  distinct(virus_id, .keep_all = TRUE) %>%
  mutate(tax_label = ifelse(
    is.na(taxonomy) | trimws(taxonomy) %in% c("", "Unclassified"),
    "Unclassified",
    tax_last(taxonomy)
  ))

host_type_lookup <- df %>%
  select(virus_id, `Host type`) %>%
  distinct(virus_id, .keep_all = TRUE) %>%
  mutate(`Host type` = replace_na(`Host type`, "Unknown"))

# Colour palettes
heatmap_col <- colorRamp2(
  c(0, 0.1, 1, 2, 3),
  c("black", "#ffffe0", "#feb24c", "#f03b20", "#bd0026")
)

quality_colors <- c(
  "High-quality"   = "#2ecc71",
  "Medium-quality" = "#3498db",
  "Low-quality"    = "#95a5a6",
  "Not-determined" = "black"
)

class_colors <- c(
  "Hexactinellida" = "#8e44ad",
  "Demospongiae"   = "#00bcd4"
)

set.seed(42)
tax_levels <- sort(unique(taxonomy_lookup$tax_label))
tax_colors <- setNames(
  colorRampPalette(c("#e41a1c","#377eb8","#4daf4a","#984ea3",
                     "#ff7f00","#a65628","#f781bf","#999999",
                     "#66c2a5","#fc8d62","#8da0cb","#e78ac3"))(length(tax_levels)),
  tax_levels
)

host_type_colors <- c(
  "Prokaryote" = "#e74c3c",
  "Eukaryote"  = "#3498db",
  "Unknown"    = "#f1c40f"
)

# Assign sponge class from sample name
assign_class <- function(x) {
  ifelse(grepl("Farrea|Aphro", x, ignore.case = TRUE), "Hexactinellida", "Demospongiae")
}

# Column labels with italic species names and roman specimen numbers
make_col_labels <- function(sample_names) {
  label_exprs <- c(
    "Aphro18"       = "italic('A. beatrix')~'specimen 18'",
    "Farrea"        = "italic('Farrea')~'sp'",
    "Bolosoma"      = "italic('Bolosoma')~'sp.'",
    "Characella"    = "italic('Characella')~'sp.'",
    "Chondrilla"    = "italic('Chondrilla')~'sp.'",
    "Ggar08"        = "italic('G. garoupa')~'specimen 08'",
    "Ggar09"        = "italic('G. garoupa')~'specimen 09'",
    "Gmega16"       = "italic('G. megastrella')~'specimen 16'",
    "Gmega17"       = "italic('G. megastrella')~'specimen 17'",
    "Paratimea"     = "italic('Paratimea')~'sp.'",
    "Pcaliculata07" = "italic('P. caliculata')~'specimen 07'",
    "Pcaliculata13" = "italic('P. caliculata')~'specimen 13'",
    "Pcaliculata14" = "italic('P. caliculata')~'specimen 14'",
    "Phycopsis"     = "italic('Phycopsis')~'sp.'",
    "Povisternata"  = "italic('P. ovisternata')"
  )
  do.call(c, lapply(sample_names, function(s) {
    if (s %in% names(label_exprs)) parse(text = label_exprs[s]) else bquote(.(s))
  }))
}

# Heatmap annotations
make_top_annot <- function(sample_names) {
  HeatmapAnnotation(
    Class                = assign_class(sample_names),
    col                  = list(Class = class_colors),
    height               = unit(0.6, "cm"),
    annotation_name_side = "right",
    annotation_name_gp   = gpar(fontsize = 9, fontface = "bold"),
    show_legend          = FALSE
  )
}

make_row_annot <- function(ids, q_colors) {
  get_vec <- function(lut, col) {
    lut %>%
      filter(virus_id %in% ids) %>%
      slice(match(ids, virus_id)) %>%
      pull({{ col }})
  }
  quals <- get_vec(quality_lookup,   checkv_quality)
  taxes <- get_vec(taxonomy_lookup,  tax_label)
  hosts <- get_vec(host_type_lookup, `Host type`)
  
  rowAnnotation(
    Quality  = quals,
    Taxonomy = taxes,
    Host     = hosts,
    col = list(
      Quality  = q_colors,
      Taxonomy = tax_colors[intersect(names(tax_colors), unique(taxes))],
      Host     = host_type_colors[intersect(names(host_type_colors), unique(hosts))]
    ),
    width                = unit(1.8, "cm"),
    show_legend          = FALSE,
    annotation_name_side = "top",
    annotation_name_gp   = gpar(fontsize = 7, fontface = "bold")
  )
}

# Legends
make_legends <- function(quality_cols, tax_cols, host_cols) {
  lgd_heatmap <- Legend(
    title         = "log10(RPKM+1)",
    col_fun       = heatmap_col,
    at            = c(0, 1, 2, 3),
    labels        = c("0", "1", "2", "3"),
    direction     = "vertical",
    title_gp      = gpar(fontsize = 9, fontface = "bold"),
    labels_gp     = gpar(fontsize = 8),
    legend_height = unit(3.5, "cm")
  )
  lgd_quality <- Legend(
    title     = "CheckV quality",
    legend_gp = gpar(fill = quality_cols),
    labels    = names(quality_cols),
    title_gp  = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
  # Viral taxonomy: family/class names in italics, Unclassified in roman
  tax_label_exprs <- parse(text = sapply(names(tax_cols), function(x) {
    if (x == "Unclassified") "'Unclassified'" else paste0("italic('", x, "')")
  }))
  lgd_taxonomy <- Legend(
    title     = "Viral taxonomy",
    legend_gp = gpar(fill = tax_cols),
    labels    = tax_label_exprs,
    title_gp  = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
  lgd_host <- Legend(
    title     = "Host type",
    legend_gp = gpar(fill = host_cols),
    labels    = names(host_cols),
    title_gp  = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
  lgd_class <- Legend(
    title     = "Sponge class",
    legend_gp = gpar(fill = class_colors),
    labels    = names(class_colors),
    title_gp  = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
  packLegend(lgd_heatmap, lgd_quality, lgd_taxonomy, lgd_host, lgd_class,
             direction = "vertical", gap = unit(0.4, "cm"))
}

# Save heatmap + legends to PDF using a two-column grid layout
save_heatmap_pdf <- function(ht, filename, pdf_width, pdf_height, legends) {
  graphics.off()
  pdf(filename, width = pdf_width, height = pdf_height)
  
  pushViewport(viewport(layout = grid.layout(
    nrow = 1, ncol = 2,
    widths = unit(c(1, 2.8), c("null", "cm"))
  )))
  
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  draw(ht, newpage = FALSE, padding = unit(c(0.3, 0.1, 0.3, 0.3), "cm"))
  popViewport()
  
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  draw(legends, x = unit(-0.2, "npc"), y = unit(0.95, "npc"), just = c("left", "top"))
  popViewport()
  
  popViewport()
  dev.off()
  cat("Saved:", filename, "\n")
}

# ----- Filtered heatmap (High- and Medium-quality vOTUs only) -----------------

quality_keep <- quality_lookup %>%
  filter(checkv_quality %in% c("High-quality", "Medium-quality")) %>%
  pull(virus_id)

mat_filt     <- mat_full[rownames(mat_full) %in% quality_keep, ]
mat_filt_log <- log10_na(mat_filt)
mat_filt_log <- mat_filt_log[rowSums(!is.na(mat_filt_log)) > 0, ]
cat("Filtered vOTUs:", nrow(mat_filt_log), "\n")

mat_filt_plot <- mat_filt_log
mat_filt_plot[mat_filt_plot == 0] <- NA

filt_ids       <- rownames(mat_filt_log)
filt_tax_cols  <- tax_colors[intersect(names(tax_colors),
                                       taxonomy_lookup  %>% filter(virus_id %in% filt_ids) %>% pull(tax_label))]
filt_host_cols <- host_type_colors[intersect(names(host_type_colors),
                                             host_type_lookup %>% filter(virus_id %in% filt_ids) %>% pull(`Host type`))]

ht_filt <- Heatmap(
  mat_filt_plot,
  name                    = "log10(RPKM+1)",
  col                     = heatmap_col,
  na_col                  = "black",
  column_title            = NULL,
  show_heatmap_legend     = FALSE,
  cluster_rows            = FALSE,
  show_row_names          = TRUE,
  row_names_side          = "left",
  row_names_gp            = gpar(fontsize = 12),
  row_names_max_width     = max_text_width(rownames(mat_filt_log), gp = gpar(fontsize = 12)),
  cluster_columns         = FALSE,
  column_order            = col_order,
  column_labels           = make_col_labels(col_order),
  show_column_names       = TRUE,
  column_names_side       = "bottom",
  column_names_gp         = gpar(fontsize = 10),
  column_names_rot        = 45,
  column_names_max_height = unit(3.5, "cm"),
  top_annotation          = make_top_annot(col_order),
  left_annotation         = make_row_annot(filt_ids,
                                           quality_colors[c("High-quality", "Medium-quality")]),
  width                   = unit(ncol(mat_filt_log) * 1.2, "cm"),
  height                  = unit(nrow(mat_filt_log) * 0.8, "cm"),
  border                  = FALSE
)

save_heatmap_pdf(
  ht_filt,
  "vOTU_RPKM_heatmap_filtered.pdf",
  pdf_width  = ncol(mat_filt_log) * 1.2 / 2.54 + 7,
  pdf_height = max(nrow(mat_filt_log) * 0.8 / 2.54 + 2, 8),
  legends    = make_legends(quality_colors[c("High-quality", "Medium-quality")],
                            filt_tax_cols, filt_host_cols)
)

# ----- Unfiltered heatmap (all quality tiers) ---------------------------------

mat_unfilt_log <- log10_na(mat_full)
mat_unfilt_log <- mat_unfilt_log[rowSums(!is.na(mat_unfilt_log)) > 0, ]
cat("Unfiltered vOTUs:", nrow(mat_unfilt_log), "\n")

mat_unfilt_plot <- mat_unfilt_log
mat_unfilt_plot[mat_unfilt_plot == 0] <- NA

unfilt_ids       <- rownames(mat_unfilt_log)
unfilt_tax_cols  <- tax_colors[intersect(names(tax_colors),
                                         taxonomy_lookup  %>% filter(virus_id %in% unfilt_ids) %>% pull(tax_label))]
unfilt_host_cols <- host_type_colors[intersect(names(host_type_colors),
                                               host_type_lookup %>% filter(virus_id %in% unfilt_ids) %>% pull(`Host type`))]

ht_unfilt <- Heatmap(
  mat_unfilt_plot,
  name                    = "log10(RPKM+1)",
  col                     = heatmap_col,
  na_col                  = "black",
  column_title            = NULL,
  show_heatmap_legend     = FALSE,
  cluster_rows            = FALSE,
  show_row_names          = TRUE,
  row_names_side          = "left",
  row_names_gp            = gpar(fontsize = 7),
  row_names_max_width     = max_text_width(rownames(mat_unfilt_log), gp = gpar(fontsize = 7)),
  cluster_columns         = FALSE,
  column_order            = col_order,
  column_labels           = make_col_labels(col_order),
  show_column_names       = TRUE,
  column_names_side       = "bottom",
  column_names_gp         = gpar(fontsize = 10),
  column_names_rot        = 45,
  column_names_max_height = unit(3.5, "cm"),
  top_annotation          = make_top_annot(col_order),
  left_annotation         = make_row_annot(unfilt_ids, quality_colors),
  width                   = unit(ncol(mat_unfilt_log) * 1.2, "cm"),
  height                  = unit(nrow(mat_unfilt_log) * 0.4, "cm"),
  border                  = FALSE
)

save_heatmap_pdf(
  ht_unfilt,
  "vOTU_RPKM_heatmap_unfiltered.pdf",
  pdf_width  = ncol(mat_unfilt_log) * 1.2 / 2.54 + 7,
  pdf_height = max(nrow(mat_unfilt_log) * 0.4 / 2.54 + 2, 12),
  legends    = make_legends(quality_colors, unfilt_tax_cols, unfilt_host_cols)
)