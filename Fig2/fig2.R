# ============================================================================
# Fig.2  Community composition: (A) phylum-level and (B) class-level stacked
#        bars, (C) Z-score heatmap of LEfSe differential taxa (21 taxa x 6
#        samples), assembled into one two-panel-row figure
# ----------------------------------------------------------------------------
# Project : kelpPR - metagenomic study of kelp aquaculture effects on
#           bacterioplankton and proteorhodopsin-based phototrophy
# Author  : Fei Zhang
# Date    : 2026-10-27
# Version : 1.0
# Input   : data/phylum_relabund.tsv, data/class_relabund.tsv,
#           data/{order,family,genus,species}_relabund.tsv (rows = taxa,
#           columns = 6 samples, relative abundance),
#           data/lefse_rep.tsv (full LEfSe feature table),
#           data/group.txt (sample / treat1 = CK | HD)
# Output  : Fig2/Fig2A_phylum_bar.{tiff,pdf}, Fig2/Fig2B_class_bar.{tiff,pdf},
#           Fig2/Fig2C_zscore_heatmap.{tiff,pdf},
#           Fig2/Fig2_final.{tiff,pdf} (190 mm double-column, 300 dpi + vector)
# Requires: ggplot2, dplyr, tidyr, ggtext, pheatmap, grid, gtable, patchwork
# ----------------------------------------------------------------------------

# ---- Step 1: package check (install missing packages from CRAN mirror) -----
need_pkgs <- c("ggplot2", "dplyr", "tidyr", "ggtext", "pheatmap", "gtable",
               "patchwork")
miss_pkgs <- need_pkgs[!is.element(need_pkgs, rownames(installed.packages()))]
if (length(miss_pkgs) > 0) {
  install.packages(miss_pkgs, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN",
                   type = "binary", quiet = TRUE)
}
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(grid)
library(patchwork)

samples <- c("CK-1", "CK-2", "CK-3", "HD-1", "HD-2", "HD-3")

# ---- Step 2: locate working directories ------------------------------------
# Run from the project root (or from Fig2/: the script walks one level up)
locate_dirs <- function(fig_id) {
  up <- if (dir.exists("data")) "." else if (dir.exists(file.path("..", "data"))) ".."
        else stop("Cannot find the project root. Set the working directory to ",
                  "the project root or to ", fig_id, "/")
  list(data_dir = file.path(up, "data"), out_dir = file.path(up, fig_id))
}
dirs <- locate_dirs("Fig2")
data_dir <- dirs$data_dir
out_dir  <- dirs$out_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Step 3: Panels A/B - phylum / class relative-abundance stacked bars
#   Wide table -> long table, merge minor taxa (top 9 + Others + Unassigned),
#   stack ordered by mean abundance (largest at the bottom)
# ============================================================================

# ---- 3.1 read data ----------------------------------------------------------
# check.names = FALSE keeps sample names like "CK-1" intact (otherwise R
# turns "-" into ".")
phylum <- read.delim(file.path(data_dir, "phylum_relabund.tsv"),
                     header = TRUE, row.names = 1, check.names = FALSE)
class  <- read.delim(file.path(data_dir, "class_relabund.tsv"),
                     header = TRUE, row.names = 1, check.names = FALSE)
grp <- read.delim(file.path(data_dir, "group.txt"), header = TRUE,
                  check.names = FALSE)

# ---- 3.2 data health check: per-sample column sums should be ~ 1 -----------
cat("== per-sample relative abundance sums (should be ~1) ==\n")
print(round(colSums(phylum), 4))
stopifnot(all(abs(colSums(phylum) - 1) < 0.01))

# ---- 3.3 reshape: wide -> long, merge minor taxa ----------------------------
# one function used for both ranks (phylum and class)
make_stack <- function(df, top_n = 9) {
  top_names <- df |> rowMeans() |> sort(decreasing = TRUE) |> names()
  keep_n <- top_names[seq_len(top_n)]
  df |>
    tibble::rownames_to_column("taxon") |>
    pivot_longer(cols = -taxon, names_to = "sample",
                 values_to = "relabund") |>
    left_join(grp, by = "sample") |>
    mutate(taxon = ifelse(is.element(taxon, keep_n), taxon,
                   ifelse(taxon == "Unassigned", "Unassigned", "Others"))) |>
    group_by(taxon, sample, treat1) |>
    summarise(relabund = sum(relabund), .groups = "drop")  # re-sum Others
}

phylum_sum <- make_stack(phylum, top_n = 9)
class_sum  <- make_stack(class,  top_n = 9)

# stacking order (mean abundance, decreasing) and sample order (CK left, HD right)
order_taxa <- function(d) {
  ord <- d |> group_by(taxon) |> summarise(m = mean(relabund), .groups = "drop") |>
    arrange(desc(m)) |> pull(taxon)
  d$taxon  <- factor(d$taxon, levels = rev(ord))
  d$sample <- factor(d$sample, levels = samples)
  d
}
phylum_sum <- order_taxa(phylum_sum)
class_sum  <- order_taxa(class_sum)

# assertion: per-sample totals still sum to 1 after merging
stopifnot(all(abs(tapply(phylum_sum$relabund, phylum_sum$sample, sum) - 1) < 0.01))

# ---- 3.4 plot ----------------------------------------------------------------
# 11-colour qualitative palette (colour-blind friendly); x-axis sample labels
# coloured by group (ggtext rich text: ggplot2 >= 4.0 does not support
# element_text(colour = <vector>), use element_markdown as recommended)
pal <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
         "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD", "#B07AA1")
x_label_colours <- c("CK-1" = "#4C72B0", "CK-2" = "#4C72B0", "CK-3" = "#4C72B0",
                     "HD-1" = "#C44E52", "HD-2" = "#C44E52", "HD-3" = "#C44E52")
coloured_x_labels <- function(x) {
  sprintf("<span style='color:%s'>%s</span>", x_label_colours[x], x)
}

plot_stack <- function(d, ttl) {
  ggplot(d, aes(x = sample, y = relabund, fill = taxon)) +
    geom_col(width = 0.75, colour = "white", linewidth = 0.15) +  # white outline
    scale_fill_manual(values = pal) +
    scale_x_discrete(labels = coloured_x_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                       labels = scales::label_percent(accuracy = 1)) +
    labs(x = NULL, y = "Relative abundance", fill = NULL, title = ttl) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1),
          legend.position = "right",
          legend.key.size = unit(0.35, "cm"),
          legend.text = element_text(size = 7),
          plot.title = element_text(size = 10))
}

p_phylum <- plot_stack(phylum_sum, "Phylum level")
p_class  <- plot_stack(class_sum,  "Class level")

# ---- 3.5 export panels (90 x 70 mm single column, 300 dpi + vector PDF) -----
ggsave(file.path(out_dir, "Fig2A_phylum_bar.tiff"), p_phylum,
       width = 90, height = 70, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "Fig2A_phylum_bar.pdf"), p_phylum,
       width = 90, height = 70, units = "mm")
ggsave(file.path(out_dir, "Fig2B_class_bar.tiff"), p_class,
       width = 90, height = 70, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "Fig2B_class_bar.pdf"), p_class,
       width = 90, height = 70, units = "mm")

cat("\n== Fig.2 panels A/B done ==\n")
cat("phylum categories:", paste(rev(levels(phylum_sum$taxon)), collapse = " | "), "\n")
cat("class categories:", paste(rev(levels(class_sum$taxon)), collapse = " | "), "\n")

# ============================================================================
# Step 4: Panel C - Z-score heatmap of LEfSe differential taxa
#   Data quirks handled below:
#   (a) lefse_rep.tsv is a full-feature table, only 26 rows carry group
#       labels, with duplicated entries (k__Bacteria appears 3x):
#       drop kingdom rank + deduplicate (keep max LDA);
#   (b) p__Proteobacteria is labelled both HD and CK - verified against the
#       measured mean direction, artefacts removed;
#   (c) LEfSe lineage writes species names "sp." as "sp__", which does not
#       match the abundance tables' "sp._" - name conversion needed.
#   Final 21 differential taxa (same set as the Fig.4 cladogram / LDA panels)
# ============================================================================

# ---- 4.1 read the LEfSe table: keep only rows with a group label -------------
lefse_raw <- read.delim(file.path(data_dir, "lefse_rep.tsv"),
                        header = TRUE, check.names = FALSE,
                        colClasses = "character")   # read as character, convert manually
lefse <- lefse_raw[is.element(lefse_raw$Groups, c("HD", "CK")), ]
lefse$LDA <- as.numeric(lefse$LDA_value)
cat("differential feature rows with group label:", nrow(lefse), "\n")

# ---- 4.2 parse lineage: rank prefix and name of the last level ---------------
# lineage looks like k__Bacteria.p__Proteobacteria.c__Alpha... - split on "."
tokens   <- strsplit(lefse$Biomarker_names, ".", fixed = TRUE)
lefse$rank  <- vapply(tokens, function(x) substr(x[length(x)], 1, 1), "")
lefse$label <- vapply(tokens, function(x) substr(x[length(x)], 4, 100), "")

lefse <- lefse[lefse$rank != "k", ]                     # quirk (a): drop kingdom
lefse <- lefse[order(lefse$Biomarker_names, -lefse$LDA), ]
lefse <- lefse[!duplicated(lefse$Biomarker_names), ]    # dedupe, keep max LDA
cat("after dedup + kingdom removal:", nrow(lefse), "differential taxa\n")

# ---- 4.3 fetch 6-sample relative abundance at the matching rank --------------
rank_files <- c(p = "phylum_relabund.tsv", c = "class_relabund.tsv",
                o = "order_relabund.tsv",  f = "family_relabund.tsv",
                g = "genus_relabund.tsv",  s = "species_relabund.tsv")
read_rank <- function(rank) {
  read.delim(file.path(data_dir, rank_files[[rank]]),
             header = TRUE, row.names = 1, check.names = FALSE)
}
# quirk (c): species "sp__HIMB1321" -> abundance table's "sp._HIMB1321"
match_row <- function(tab, name) {
  if (is.element(name, rownames(tab))) return(name)
  alt <- gsub("__", "._", name)
  if (is.element(alt, rownames(tab))) return(alt)
  NA_character_
}
abund_list <- lapply(seq_len(nrow(lefse)), function(i) {
  tab  <- read_rank(lefse$rank[i])
  name <- match_row(tab, lefse$label[i])
  if (is.na(name)) {
    cat("  !! not found in abundance tables:", lefse$label[i],
        "(rank", lefse$rank[i], ")\n")
    return(NULL)
  }
  cbind(data.frame(taxon = name, rank = lefse$rank[i], label = lefse$label[i],
                   group = lefse$Groups[i]),
        as.data.frame(tab[name, samples, drop = FALSE]), row.names = NULL)
})
abund <- do.call(rbind, abund_list)

# ---- 4.4 direction check: measured mean direction vs LEfSe label (quirk b) ---
ck_mean <- rowMeans(abund[, c("CK-1", "CK-2", "CK-3")])
hd_mean <- rowMeans(abund[, c("HD-1", "HD-2", "HD-3")])
actual_dir <- ifelse(hd_mean > ck_mean, "HD", "CK")
mismatch <- actual_dir != abund$group
if (any(mismatch)) {
  cat("direction conflicts (removed based on measured data):\n")
  for (i in which(mismatch)) {
    cat(sprintf("  %s (LEfSe label %s, measured CK=%.4f vs HD=%.4f)\n",
                abund$taxon[i], abund$group[i], ck_mean[i], hd_mean[i]))
  }
}
abund <- abund[!mismatch, ]
cat("differential taxa included in the heatmap:", nrow(abund), "\n")
stopifnot(nrow(abund) >= 15)

# ---- 4.5 build matrix and Z-score rows ----------------------------------------
mat <- as.matrix(abund[, samples])
storage.mode(mat) <- "numeric"

# row display names: mark Unclassified with its rank; "sp__" -> "sp. ";
# underscores restored to spaces
disp <- ifelse(abund$label == "Unclassified",
               paste0("Unclassified (", c(p = "phylum", c = "class", o = "order",
                                          f = "family", g = "genus", s = "species")[abund$rank], ")"),
               abund$label)
disp <- gsub("sp__", "sp. ", disp)
disp <- gsub("_", " ", disp)
stopifnot(!anyDuplicated(disp))

# Z-score = (x - row mean) / row SD: scale() works column-wise, so transpose twice
z <- t(scale(t(mat)))
stopifnot(all(abs(z) < 10))

# row order: HD-enriched block on top, CK block below, LDA decreasing within blocks
lda_vec <- vapply(seq_len(nrow(abund)), function(i) {
  m <- which(lefse$label == abund$label[i] & lefse$rank == abund$rank[i])[1]
  lefse$LDA[m]
}, 0)
ord <- order(abund$group == "CK", -lda_vec)
z        <- z[ord, ]
disp     <- disp[ord]
rank_ord <- abund$rank[ord]
grp_row  <- abund$group[ord]
rownames(z) <- disp

# ---- 4.6 italicise taxon names (plotmath expression labels) -------------------
# conventions: genus and species epithet italic; "sp.", strain IDs, higher
# ranks and Unclassified roman; "Candidatus" italic but the following candidiate
# name roman. parse(text = ...) requires valid R expressions - plain text must
# be wrapped as string literals '"text"'
make_label <- function(name, rank) {
  if (grepl("^Unclassified", name)) return(sprintf('"%s"', name))
  if (grepl(" bacterium$", name))   return(sprintf('"%s"', name))
  if (grepl("^Candidatus ", name))
    return(sprintf('italic("Candidatus")~"%s"', sub("^Candidatus ", "", name)))
  if (rank == "g") return(sprintf('italic("%s")', name))
  if (rank == "s") {
    if (grepl(" sp\\. ", name)) {
      parts <- strsplit(name, " ", fixed = TRUE)[[1]]
      return(sprintf('italic("%s")~"%s"', parts[1],
                     paste(parts[-1], collapse = " ")))
    }
    return(sprintf('italic("%s")', name))
  }
  sprintf('"%s"', name)
}
lbl_txt <- vapply(seq_along(disp),
                  function(i) make_label(disp[i], rank_ord[i]), "")
labels_row_expr <- parse(text = lbl_txt)
cat("italicised rows (genus/species level:", sum(grepl("italic", lbl_txt)), ")\n")

# ---- 4.7 pheatmap --------------------------------------------------------------
ann_col <- data.frame(Habitat = factor(c(rep("CK", 3), rep("HD", 3)),
                                       levels = c("CK", "HD")),
                      row.names = samples)
ann_colors <- list(Habitat = c(CK = "#4C72B0", HD = "#C44E52"))
ann_row <- data.frame(Enriched = factor(grp_row, levels = c("HD", "CK")),
                      row.names = disp)
ann_colors$Enriched <- c(CK = "#4C72B0", HD = "#C44E52")

# diverging colour scale: blue-white-red (red = high z, echoing the HD colour)
zmax <- max(abs(z))
brks <- seq(-zmax, zmax, length.out = 101)
pal  <- colorRampPalette(c("#4C72B0", "white", "#C44E52"))(100)

# silent = TRUE grabs the gtable; gtable_add_padding adds margins on all sides
# (protects annotation legends from clipping); export via self-managed devices
ph <- pheatmap(z,
               color = pal, breaks = brks,
               cluster_rows = FALSE, cluster_cols = FALSE,
               annotation_col = ann_col, annotation_row = ann_row,
               annotation_colors = ann_colors,
               labels_row = labels_row_expr,
               border_color = NA,
               fontsize = 9, fontsize_row = 8,
               cellwidth = 14, cellheight = 10,
               angle_col = 45,
               main = "Differential taxa (LEfSe, LDA >= 4)",
               silent = TRUE)
g <- gtable::gtable_add_padding(ph$gtable, unit(4, "mm"))

grid.newpage()
grid.draw(g)

# ---- 4.8 export panel -----------------------------------------------------------
tiff(file.path(out_dir, "Fig2C_zscore_heatmap.tiff"),
     width = 1850, height = 1300, res = 300)
grid.newpage(); grid.draw(g); dev.off()

pdf(file.path(out_dir, "Fig2C_zscore_heatmap.pdf"), width = 6.2, height = 4.4)
grid.newpage(); grid.draw(g); dev.off()

cat("\n== Fig.2 panel C done ==\n")
cat("differential taxa: HD side", sum(grp_row == "HD"),
    "/ CK side", sum(grp_row == "CK"), "\n")

# ============================================================================
# Step 5: assemble the final figure
#   Two object types: ggplot (stacked bars) and gtable (pheatmap). patchwork
#   wraps the gtable with wrap_elements(full = ...); plot_annotation(
#   tag_levels = "A") adds the A/B/C panel tags. Row 1: A + B (equal width),
#   row 2: C full width.
# ============================================================================

# shrink stacked-bar legend fonts to fit the narrower sub-panels
p_phylum <- p_phylum + theme(legend.text = element_text(size = 6),
                             legend.title = element_text(size = 7))
p_class  <- p_class  + theme(legend.text = element_text(size = 6),
                             legend.title = element_text(size = 7))

fig2 <- (p_phylum | p_class) /
        wrap_elements(full = g) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold"))

print(fig2)

# ---- 5.1 export (190 mm double-column full width) ------------------------------
ggsave(file.path(out_dir, "Fig2_final.tiff"), fig2,
       width = 190, height = 210, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Fig2_final.pdf"), fig2,
       width = 190, height = 210, units = "mm", device = cairo_pdf)

cat("\n== Fig.2 final figure done ==\nOutput:", out_dir, "/Fig2_final.{tiff,pdf}\n")
