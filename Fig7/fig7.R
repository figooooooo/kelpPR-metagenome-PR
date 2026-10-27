# ============================================================================
# Fig.7  Proteorhodopsin genetic potential and composition (three panels):
#        (A) summed PR TPM / PR GPM boxplots;
#        (B) enriched-side contribution to each sample's PR pool;
#        (C) spectral-type composition, TPM-weighted
# ----------------------------------------------------------------------------
# Project : kelpPR - metagenomic study of kelp aquaculture effects on
#           bacterioplankton and proteorhodopsin-based phototrophy
# Author  : Fei Zhang
# Date    : 2026-10-27
# Version : 1.0
# Input   : Fig7/pr_gpm_summary.tsv (sample, group, PR_TPM_sum, PR_GPM, ...),
#           Fig7/pr_annotation_all.tsv (gene_id, ..., spectral, ...; one or
#           more rows per gene),
#           data/pr_tpm.tsv (gene_id, CK-1..HD-3, sum; 139 genes)
# Output  : Fig7/fig7a_tpm_gpm.{tiff,pdf}, Fig7/fig7b_pool_composition.{tiff,pdf},
#           Fig7/fig7c_spectral.{tiff,pdf},
#           Fig7/Fig7_pr_potential.{tiff,pdf} (170 x 150 mm assembly)
# Requires: ggplot2, patchwork
# ----------------------------------------------------------------------------
# Panel A is the main result: total PR abundance CK > HD. Boxplots with the
# n = 3 points drawn on top (grey, no jitter needed), facetted TPM | GPM;
# both metrics give exact MW p = 0.100 (complete separation, the power floor
# of n = 3 vs 3), annotated uniformly. Values are taken directly from
# pr_gpm_summary.tsv (server-fixed) rather than re-summing pr_tpm.tsv, to
# stay consistent with the GPM normalisation (recA/rpoB/rplB).
# Panel B is the mechanism: 139 genes split by mean(HD) vs mean(CK)
# (expected 61 HD / 78 CK, cross-checked in the console); each sample's PR
# pool = HD-enriched contribution + CK-enriched contribution = 100%.
# Panel C: spectral types (TPM-weighted); annotations matched by gene_id
# first-hit; the type set must be green / blue / unknown (manuscript:
# 128 green + 10 blue + 1 unknown); blue-light shares are printed for the
# manuscript's Section 3.6 cross-check. The x axis holds only the 6 samples,
# no mean reference column.
# ----------------------------------------------------------------------------

# ---- Step 1: package check ---------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

# ---- Step 2: locate working directories ---------------------------------------
locate_dirs <- function(fig_id) {
  up <- if (dir.exists("data")) "." else if (dir.exists(file.path("..", "data"))) ".."
        else stop("Cannot find the project root. Set the working directory to ",
                  "the project root or to ", fig_id, "/")
  list(data_dir = file.path(up, "data"), out_dir = file.path(up, fig_id))
}
dirs <- locate_dirs("Fig7")
data_dir <- dirs$data_dir
out_dir  <- dirs$out_dir

samples <- c("CK-1", "CK-2", "CK-3", "HD-1", "HD-2", "HD-3")
COL_HD <- "#C44E52"; COL_CK <- "#4C72B0"
COL_GREEN <- "#3C9A5F"; COL_BLUE <- "#3C5488"; COL_UNK <- "grey60"

# ---- Step 3: read data ----------------------------------------------------------
gpm <- read.delim(file.path(out_dir, "pr_gpm_summary.tsv"),
                  header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
tpm <- read.delim(file.path(data_dir, "pr_tpm.tsv"),
                  header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
ann <- read.delim(file.path(out_dir, "pr_annotation_all.tsv"),
                  header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(gpm) == 6, nrow(tpm) == 139)

# ============================================================================
# Step 4: Panel A - summed PR TPM / PR GPM boxplots
# ============================================================================
dfA <- data.frame(
  sample = rep(gpm$sample, 2),
  group  = rep(gpm$group, 2),
  metric = rep(c("summed PR TPM", "PR GPM"), each = 6),
  value  = c(gpm$PR_TPM_sum, gpm$PR_GPM),
  stringsAsFactors = FALSE)
dfA$metric <- factor(dfA$metric, levels = c("summed PR TPM", "PR GPM"))
dfA$group  <- factor(dfA$group, levels = c("CK", "HD"))

# hard assertion: CK > HD complete separation (both metrics, same direction)
for (v in list(gpm$PR_TPM_sum, gpm$PR_GPM)) {
  stopifnot(all(v[1:3] > v[4:6]))
}

pA <- ggplot(dfA, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.52, outlier.shape = NA, alpha = 0.85,
               colour = "grey25", linewidth = 0.35) +
  # points vertically aligned (no jitter needed at n = 3), dark grey so the
  # data layer reads above the coloured boxes
  geom_point(colour = "grey15", size = 1.6, show.legend = FALSE) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(CK = COL_CK, HD = COL_HD), guide = "none") +
  labs(x = NULL, y = NULL) +
  annotate("text", x = 1.5, y = Inf, label = "exact p = 0.100",
           size = 2.4, colour = "grey30", vjust = 1.6) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 7.5, colour = "black"),
        strip.text = element_text(size = 8.5, face = "bold"),
        strip.background = element_blank())

ggsave(file.path(out_dir, "fig7a_tpm_gpm.tiff"), pA,
       width = 120, height = 62, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig7a_tpm_gpm.pdf"), pA,
       width = 120, height = 62, units = "mm", device = cairo_pdf)

# ============================================================================
# Step 5: Panel B - PR pool contribution by enriched side x sample
# ============================================================================
mat <- as.matrix(tpm[, samples]); storage.mode(mat) <- "numeric"
hd_mean <- rowMeans(mat[, 4:6]); ck_mean <- rowMeans(mat[, 1:3])
enrich <- ifelse(hd_mean > ck_mean, "HD", "CK")
cat(sprintf("enrichment split: HD-enriched %d / CK-enriched %d (manuscript: 61/78)\n",
            sum(enrich == "HD"), sum(enrich == "CK")))

pool <- rbind(HD = colSums(mat[enrich == "HD", ]),
              CK = colSums(mat[enrich == "CK", ]))
prop <- sweep(pool, 2, colSums(pool), "/")
cat(sprintf("mean HD-enriched gene contribution in HD samples: %.1f%% (manuscript ~60%%)\n",
            100 * mean(prop["HD", 4:6])))
cat(sprintf("mean CK-enriched gene contribution in CK samples: %.1f%% (manuscript ~82%%)\n",
            100 * mean(prop["CK", 1:3])))

# pitfall: as.vector(matrix) reads column-major and would interleave the two
# rows - build the vector explicitly row by row
dfB <- data.frame(
  sample = factor(rep(samples, 2), levels = samples),
  side   = factor(rep(c("HD-enriched genes", "CK-enriched genes"), each = 6),
                  levels = c("HD-enriched genes", "CK-enriched genes")),
  prop   = c(prop["HD", ], prop["CK", ]),
  stringsAsFactors = FALSE)
stopifnot(all(abs(tapply(dfB$prop, dfB$sample, sum) - 1) < 1e-9))

pB <- ggplot(dfB, aes(x = sample, y = prop, fill = side)) +
  geom_col(width = 0.62, colour = "grey25", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * prop), group = side),
            position = position_stack(vjust = 0.5), size = 2.2,
            colour = "white") +
  scale_fill_manual(values = c("HD-enriched genes" = COL_HD,
                               "CK-enriched genes" = COL_CK), name = NULL) +
  scale_y_continuous(labels = function(v) sprintf("%.0f%%", 100 * v),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "share of sample PR pool") +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 7.5, colour = "black"),
        axis.title.y = element_text(size = 8),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.margin = margin(t = -2))

ggsave(file.path(out_dir, "fig7b_pool_composition.tiff"), pB,
       width = 80, height = 85, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig7b_pool_composition.pdf"), pB,
       width = 80, height = 85, units = "mm", device = cairo_pdf)

# ============================================================================
# Step 6: Panel C - spectral-type composition (TPM-weighted)
# ============================================================================
sp <- ann$spectral[match(tpm$gene_id, ann$gene_id)]
stopifnot(!any(is.na(sp)))
spec_sum <- sapply(unique(sp), function(s) colSums(mat[sp == s, , drop = FALSE]))
rownames(spec_sum) <- samples
spec_prop <- sweep(spec_sum, 1, rowSums(spec_sum), "/")
cat("spectral type set:", paste(colnames(spec_prop), collapse = ", "),
    sprintf("(gene counts: %s)", paste(names(table(sp)), table(sp), sep = "=", collapse = ", ")), "\n")
cat(sprintf("mean blue-light share: HD samples %.1f%%; CK samples %.1f%%\n",
            100 * mean(spec_prop[4:6, "blue"]), 100 * mean(spec_prop[1:3, "blue"])))

dfC <- data.frame(
  sample = factor(rep(samples, ncol(spec_prop)), levels = samples),
  spectral = factor(rep(colnames(spec_prop), each = 6),
                    levels = c("green", "blue", "unknown")),
  prop = as.vector(spec_prop[, c("green", "blue", "unknown")]),
  stringsAsFactors = FALSE)
# labels on the grey (unknown) segment use dark text (white on grey60 is only
# 2.8:1 contrast, dark reaches 4.6:1); green/blue segments keep white text
dfC$txt_col <- ifelse(dfC$spectral == "unknown", "grey15", "white")

pC <- ggplot(dfC, aes(x = sample, y = prop, fill = spectral)) +
  # position_stack(reverse = TRUE) locks the stacking order: green at the
  # bottom, blue in the middle, unknown on top (ggplot2 semantics:
  # reverse = TRUE puts the first level at the bottom)
  geom_col(width = 0.62, colour = "grey25", linewidth = 0.2,
           position = position_stack(reverse = TRUE)) +
  geom_text(data = dfC[dfC$prop > 0.02, ],
            aes(label = sprintf("%.1f%%", 100 * prop), group = spectral,
                colour = txt_col),
            position = position_stack(vjust = 0.5, reverse = TRUE), size = 2.0,
            show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_manual(values = c(green = COL_GREEN, blue = COL_BLUE,
                               unknown = COL_UNK), name = NULL) +
  scale_y_continuous(labels = function(v) sprintf("%.0f%%", 100 * v),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "share of sample PR pool (TPM-weighted)") +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 7.5, colour = "black"),
        axis.title.y = element_text(size = 8),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.margin = margin(t = -2))

ggsave(file.path(out_dir, "fig7c_spectral.tiff"), pC,
       width = 80, height = 85, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig7c_spectral.pdf"), pC,
       width = 80, height = 85, units = "mm", device = cairo_pdf)

# ============================================================================
# Step 7: assemble the final figure
#   Layout: row 1 panel A (TPM + GPM facet strip), row 2 B | C;
#   tag_levels = "A" adds the A/B/C panel tags automatically.
# ============================================================================
fig7 <- pA / (pB | pC) +
  plot_layout(heights = c(1, 1.35)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold"))

print(fig7)

ggsave(file.path(out_dir, "Fig7_pr_potential.tiff"), fig7,
       width = 170, height = 150, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Fig7_pr_potential.pdf"), fig7,
       width = 170, height = 150, units = "mm", device = cairo_pdf)

cat("\n== Fig.7 done ->", out_dir, "==\n")
