# ============================================================================
# Fig.5  Functional potential of kelp-farming (HD) vs control (CK)
#        metagenomes: (A) KEGG L3 log2FC vs mean abundance,
#        (B) CAZy families implicated in kelp-polysaccharide degradation
# ----------------------------------------------------------------------------
# Project : kelpPR - metagenomic study of kelp aquaculture effects on
#           bacterioplankton and proteorhodopsin-based phototrophy
# Author  : Fei Zhang
# Date    : 2026-10-27
# Version : 1.0
# Input   : Fig5/kegg_l3_diff.tsv (170 KEGG L3 pathways: group means,
#           log2FC, exact Mann-Whitney p),
#           Fig5/kelp_polysaccharide_families.tsv (CAZy families, HD/CK
#           ratios and tests)
# Output  : Fig5/Fig5_KEGG_CAZy.{tiff,pdf} (190 mm, 300 dpi + vector)
# Requires: ggplot2, patchwork
# ----------------------------------------------------------------------------

# ---- Step 1: package check (install missing packages from CRAN mirror) -----
need_pkgs <- c("ggplot2", "patchwork")
miss_pkgs <- need_pkgs[!is.element(need_pkgs, rownames(installed.packages()))]
if (length(miss_pkgs) > 0) {
  install.packages(miss_pkgs, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN",
                   type = "binary", quiet = TRUE)
}
library(ggplot2)
library(patchwork)
theme_set(theme_classic(base_size = 10))

# ---- Step 2: locate working directories ------------------------------------
# Run from the project root (or from Fig5/: the script walks one level up).
# Figure-specific inputs live in Fig5/.
locate_dirs <- function(fig_id) {
  up <- if (dir.exists("data")) "." else if (dir.exists(file.path("..", "data"))) ".."
        else stop("Cannot find the project root. Set the working directory to ",
                  "the project root or to ", fig_id, "/")
  list(data_dir = file.path(up, "data"), out_dir = file.path(up, fig_id))
}
dirs <- locate_dirs("Fig5")
data_dir <- dirs$data_dir
out_dir  <- dirs$out_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Step 3: Panel A - KEGG L3 overview (volcano-plot variant)
#   x = mean relative abundance (log10), y = log2FC, colour = KEGG level-1
#   category. Statistics: n = 3 vs 3 exact Mann-Whitney; complete separation
#   gives p = 0.100, the floor of this design - the annotation therefore
#   states both "94/170 completely separated" and "none significant after
#   BH-FDR", so the panel is not over-read (functional redundancy is the
#   main message).
# ============================================================================

# fileEncoding = "UTF-8-BOM": the vendor-delivered tables carry a BOM;
# without the explicit declaration the first header becomes garbled on
# Windows (\ufeff prefix)
kegg <- read.delim(file.path(out_dir, "kegg_l3_diff.tsv"),
                   fileEncoding = "UTF-8-BOM", check.names = FALSE,
                   stringsAsFactors = FALSE)
kegg$mean_all <- (kegg$mean_CK + kegg$mean_HD) / 2
kegg$log2FC <- suppressWarnings(as.numeric(kegg$log2FC))
kegg <- kegg[is.finite(kegg$log2FC) & kegg$mean_all > 0, ]  # drop 3 one-side-zero pathways
n_path <- nrow(kegg)
n_sep <- sum(kegg$p_exactMW <= 0.100)
cat("KEGG L3 plotted:", n_path, "| complete separation:", n_sep, "\n")
stopifnot(n_path == 170, n_sep == 94)     # matches the corrected manuscript counts

pal_l1 <- c("Metabolism" = "#E64B35",
            "Genetic Information Processing" = "#4DBBD5",
            "Environmental Information Processing" = "#00A087",
            "Cellular Processes" = "#3C5488",
            "Organismal Systems" = "#F39B7F",
            "Human Diseases" = "#8491B4",
            "Brite Hierarchies" = "#91D1C2")
pal_l1 <- pal_l1[intersect(names(pal_l1), unique(kegg[["KEGG_pathway level 1"]]))]

pA <- ggplot(kegg, aes(x = mean_all, y = log2FC,
                       color = `KEGG_pathway level 1`)) +
  geom_hline(yintercept = c(-0.25, 0.25), linetype = "dashed",
             color = "grey70", linewidth = 0.35) +
  geom_point(size = 1.8, alpha = 0.9) +
  scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2),
                labels = c("0.0001", "0.001", "0.01")) +
  scale_color_manual(values = pal_l1, name = "KEGG level 1") +
  labs(x = "Mean relative abundance (proportion of annotated reads)",
       y = expression(log[2]*" FC (HD / CK)")) +
  annotate("text", x = min(kegg$mean_all), y = max(kegg$log2FC) + 0.32,
           hjust = 0, size = 2.7,
           label = paste0(n_sep, "/", n_path,
                          " pathways separated (exact p = 0.100); none at BH-FDR q < 0.10")) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.25))) +
  coord_cartesian(clip = "off") +
  theme(legend.position = "right",
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 9))

# ============================================================================
# Step 4: Panel B - kelp-polysaccharide CAZy families
#   Families ordered by log2FC as a two-directional bar plot; asterisks mark
#   completely separated families. GH16/GH17 (laminarinases) are always
#   included.
# ============================================================================
kelp <- read.delim(file.path(out_dir, "kelp_polysaccharide_families.tsv"),
                   fileEncoding = "UTF-8-BOM", check.names = FALSE,
                   stringsAsFactors = FALSE)
kelp$lfc <- log2(kelp$ratio_HD_CK)
sel <- kelp[kelp$p_exactMW <= 0.100 |
              is.element(kelp$family, c("GH16", "GH17")), ]  # always keep GH16/17
sel <- sel[order(sel$lfc), ]
sel$lab <- paste0(sel$family, ifelse(sel$p_exactMW <= 0.100, " *", ""))
sel$family <- factor(sel$lab, levels = sel$lab)

pB <- ggplot(sel, aes(x = lfc, y = family, fill = direction)) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.35) +
  geom_col(width = 0.7, color = "grey25", linewidth = 0.15) +
  geom_text(aes(label = sprintf("%.2f", ratio_HD_CK)),
            hjust = ifelse(sel$lfc >= 0, -0.25, 1.25), size = 2.4) +
  scale_fill_manual(values = c(HD = "#C44E52", CK = "#4C72B0"),
                    labels = c(HD = "HD-enriched", CK = "CK-enriched"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.28, 0.28))) +
  labs(x = expression(log[2]*" FC (HD / CK)"),
       y = "CAZy family",
       caption = "* complete separation between groups (exact p = 0.100, n = 3 vs 3)") +
  theme(axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 9),
        legend.position = "bottom",
        legend.margin = margin(t = -2),
        plot.caption = element_text(size = 6.5, color = "grey35"))

# ---- Step 5: assemble and export --------------------------------------------
fig5 <- pA + pB + plot_layout(widths = c(1.25, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold"))

print(fig5)

ggsave(file.path(out_dir, "Fig5_KEGG_CAZy.tiff"), fig5,
       width = 190, height = 95, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "Fig5_KEGG_CAZy.pdf"), fig5,
       width = 190, height = 95, units = "mm", device = cairo_pdf)

cat("\n== Fig.5 done ==\nOutput:", out_dir, "/Fig5_KEGG_CAZy.{tiff,pdf}\n")
