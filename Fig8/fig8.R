# ============================================================================
# Fig.8  MAG catalogue and PR-rhodopsin host linkage (four panels):
#        (A) CheckM completeness / contamination;
#        (B) habitat log2FC by phylum;
#        (C) PR-carrying / total MAGs per phylum;
#        (D) bin.27 coverage vs OTU110 TPM
# ----------------------------------------------------------------------------
# Project : kelpPR - metagenomic study of kelp aquaculture effects on
#           bacterioplankton and proteorhodopsin-based phototrophy
# Author  : Fei Zhang
# Date    : 2026-10-27
# Version : 1.0
# Input   : Fig8/mag_hd_vs_ck.tsv (131 rows: MAG, GTDB_taxonomy,
#             Completeness, Contamination, meanTPM_CK, meanTPM_HD,
#             log2FC_HD_CK, exact_MW_p),
#           Fig8/mag_pr_summary.tsv (63 rows: MAG, n_PR, PR_proteins, ...),
#           Fig8/bin_abundance_table.tab (Genomic bins, S01_R0001-6;
#             R0001-3 = CK, R0004-6 = HD),
#           data/pr_tpm.tsv (OTU110 centroid = geneset000191922)
# Output  : Fig8/fig8a_checkm.{tiff,pdf}, fig8b_habitat_lfc.{tiff,pdf},
#           fig8c_pr_host.{tiff,pdf}, fig8d_bin27_otu110.{tiff,pdf},
#           Fig8/Fig8_mag_rhodopsin.{tiff,pdf} (180 x 130 mm assembly)
# Requires: ggplot2, cowplot
# ----------------------------------------------------------------------------
# Panel A: 131 MAGs, completeness (x) vs contamination (y), MIMAG threshold
#   lines x = 70/90 and y = 5/10; tiers HQ (>=90 & <5) / MQ (>=80 & <10) / LQ;
#   hard assertions HQ 51 / MQ 45 / LQ 35, PR-carrying 63.
#   Pitfall: the log2FC_HD_CK column carries a "+" prefix - strip with
#   sub("^\\+", "") before as.numeric.
# Panel B: x = log2FC, y = phylum (phyla with < 4 MAGs pooled into
#   "Other (7 phyla)", ordered by median lfc), jittered points coloured by
#   enrichment direction (|log2FC| >= 1 -> Kelp farm / Control, else NS),
#   n = labels on the right.
# Panel C: PR-carrying / total horizontal bars per phylum (total bar in the
#   PR colour underneath, total - n_pr overdrawn in grey); y order matches
#   panel B (ord_b) so B and C read as a pair.
# Panel D: bin.27 (Opacimonas) coverage bars + OTU110 (= geneset000191922)
#   TPM diamonds on a dotted line, dual y axis with a fixed scale factor
#   (sc = 4), ylim 0-55, the secondary axis converts back to TPM.
# ----------------------------------------------------------------------------

# ---- Step 1: package check ---------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(cowplot) })
theme_set(theme_classic(base_size = 10))

# ---- Step 2: locate working directories ---------------------------------------
locate_dirs <- function(fig_id) {
  up <- if (dir.exists("data")) "." else if (dir.exists(file.path("..", "data"))) ".."
        else stop("Cannot find the project root. Set the working directory to ",
                  "the project root or to ", fig_id, "/")
  list(data_dir = file.path(up, "data"), out_dir = file.path(up, fig_id))
}
dirs <- locate_dirs("Fig8")
data_dir <- dirs$data_dir
out_dir  <- dirs$out_dir

col_hd <- "#C44E52"   # kelp farm
col_ck <- "#4C72B0"   # control
col_pr <- "#00A087"   # rhodopsin-carrying

# ---- Step 3: read data, tiers, assertions ---------------------------------------
mag <- read.delim(file.path(out_dir, "mag_hd_vs_ck.tsv"),
                  fileEncoding = "UTF-8-BOM", check.names = FALSE,
                  stringsAsFactors = FALSE)
mag$lfc  <- suppressWarnings(as.numeric(sub("^\\+", "", mag$log2FC_HD_CK)))
stopifnot(nrow(mag) == 131, sum(is.finite(mag$lfc)) == 131)

pr <- read.delim(file.path(out_dir, "mag_pr_summary.tsv"),
                 fileEncoding = "UTF-8-BOM", check.names = FALSE,
                 stringsAsFactors = FALSE)
stopifnot(nrow(pr) == 63)

# quality tier (same rule as the MAG catalogue supplementary table)
tier_of <- function(comp, cont) {
  if (comp >= 90 && cont < 5) return("HQ")
  if (comp >= 80 && cont < 10) return("MQ")
  return("LQ")
}
mag$tier <- mapply(tier_of, mag$Completeness, mag$Contamination)
mag$has_pr <- is.element(mag$MAG, pr$MAG)

# phylum, pooling phyla with < 4 MAGs
ph <- sub("^p__", "", sapply(strsplit(mag$GTDB_taxonomy, ";"),
                             function(x) trimws(x[2])))
names(ph) <- mag$MAG
ph[ph == "Thermoplasmatota"] <- "Thermoplasmatota\n(Poseidoniia)"  # two-line label
mag$phylum <- ph
n_ph <- table(mag$phylum)
mag$phylum2 <- ifelse(n_ph[mag$phylum] >= 4, mag$phylum, "Other (7 phyla)")

# enrichment direction (same rule as the supplementary table)
mag$direction <- ifelse(mag$lfc >= 1, "Kelp farm",
                 ifelse(mag$lfc <= -1, "Control", "NS"))

# hard assertions: tier and PR-MAG counts (final catalogue numbers)
tt <- table(mag$tier)
stopifnot(tt["HQ"] == 51, tt["MQ"] == 45, tt["LQ"] == 35, sum(mag$has_pr) == 63)
cat("tier counts: HQ=", tt["HQ"], " MQ=", tt["MQ"], " LQ=", tt["LQ"],
    " | PR-carrying: ", sum(mag$has_pr), "\n", sep = "")
cat("direction counts: ", paste(names(table(mag$direction)), table(mag$direction),
    sep = "=", collapse = ", "), "\n", sep = "")

# ---- Step 4: Panel A - CheckM quality plot --------------------------------------
n_t <- table(mag$tier)
lab_t <- paste0(names(n_t), " (", as.integer(n_t), ")")
names(lab_t) <- names(n_t)
mag$tier <- factor(mag$tier, levels = c("HQ", "MQ", "LQ"))

pA <- ggplot(mag, aes(x = Completeness, y = Contamination, color = tier)) +
  geom_hline(yintercept = c(5, 10), linetype = "dashed",
             color = "grey65", linewidth = 0.35) +
  geom_vline(xintercept = c(70, 90), linetype = "dashed",
             color = "grey65", linewidth = 0.35) +
  geom_point(size = 1.7, alpha = 0.9, stroke = 0.2) +
  scale_color_manual(values = c("HQ" = "#00A087", "MQ" = "#F39B7F",
                               "LQ" = "#8491B4"),
                    labels = lab_t, name = "MIMAG tier",
                    drop = FALSE) +
  coord_cartesian(xlim = c(60, 100), ylim = c(0, 11.5), clip = "off") +
  scale_x_continuous(breaks = c(60, 80, 100)) +
  scale_y_continuous(breaks = c(0, 3, 6, 9)) +
  labs(x = "Completeness (%)", y = "Contamination (%)") +
  theme(legend.position = "right",
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        axis.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 7.3),
        axis.title = element_text(size = 9))

ggsave(file.path(out_dir, "fig8a_checkm.tiff"), pA,
       width = 95, height = 62, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig8a_checkm.pdf"), pA,
       width = 95, height = 62, units = "mm", device = cairo_pdf)

# ---- Step 5: Panel B - habitat preference by phylum -------------------------------
ord_b <- names(sort(tapply(mag$lfc, mag$phylum2, median)))
mag$phb <- factor(mag$phylum2, levels = ord_b)
n_b <- as.integer(table(factor(mag$phylum2, levels = ord_b)))

pB <- ggplot(mag, aes(x = lfc, y = phb, color = direction)) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.35) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             color = "grey80", linewidth = 0.3) +
  geom_jitter(height = 0.28, size = 1.5, alpha = 0.85,
              show.legend = TRUE) +
  # legend title kept short ("Enriched in"); the |log2FC| >= 1 threshold is
  # stated in the manuscript caption
  scale_color_manual(values = c("Kelp farm" = col_hd, "Control" = col_ck,
                                "NS" = "grey55"),
                     name = "Enriched in",
                     drop = FALSE) +
  annotate("text", x = max(mag$lfc) + 1.1, y = seq_along(ord_b),
           label = paste0("n=", n_b), hjust = 0, size = 2.6,
           color = "grey30") +
  coord_cartesian(xlim = c(min(mag$lfc) - 0.4, max(mag$lfc) + 2.9),
                  clip = "off") +
  labs(x = expression(log[2]*" FC (HD / CK)"), y = NULL) +
  theme(legend.position = "right",
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7.5),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 9))

ggsave(file.path(out_dir, "fig8b_habitat_lfc.tiff"), pB,
       width = 95, height = 62, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig8b_habitat_lfc.pdf"), pB,
       width = 95, height = 62, units = "mm", device = cairo_pdf)

# ---- Step 6: Panel C - PR-carrying MAGs by phylum ----------------------------------
pr_ph <- table(mag$phylum2, mag$has_pr)
pr_df <- data.frame(phylum = ord_b,
                   total = as.integer(table(factor(mag$phylum2,
                                                   levels = ord_b))),
                   n_pr  = as.integer(pr_ph[ord_b, "TRUE"]))
pr_df$phc <- factor(pr_df$phylum, levels = ord_b)

pC <- ggplot(pr_df, aes(x = total, y = phc)) +
  geom_col(aes(fill = "PR-carrying"), width = 0.62,
           color = "grey25", linewidth = 0.15) +
  geom_col(aes(x = total - n_pr, fill = "No PR"), width = 0.62,
           color = "grey25", linewidth = 0.15) +
  geom_text(aes(x = total + 1.2, label = paste0(n_pr, "/", total)),
            hjust = 0, size = 2.6, color = "grey30") +
  scale_fill_manual(values = c("PR-carrying" = col_pr,
                               "No PR" = "grey82"), name = NULL) +
  scale_x_continuous(breaks = c(0, 20, 40, 60)) +
  coord_cartesian(xlim = c(0, 82), clip = "off") +
  labs(x = "Number of MAGs", y = NULL,
       caption = "PR-carrying / total MAGs per phylum") +
  theme(legend.position = "right",
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 7),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 9),
        plot.caption = element_text(size = 6.5, color = "grey35"))

ggsave(file.path(out_dir, "fig8c_pr_host.tiff"), pC,
       width = 95, height = 62, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig8c_pr_host.pdf"), pC,
       width = 95, height = 62, units = "mm", device = cairo_pdf)

# ---- Step 7: Panel D - Opacimonas MAG coverage vs OTU110 TPM -----------------------
cov <- read.delim(file.path(out_dir, "bin_abundance_table.tab"),
                  check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("bin.27" %in% cov$`Genomic bins`)
v27 <- as.numeric(cov[cov$`Genomic bins` == "bin.27",
                       paste0("S01_R000", 1:6)])

tpm <- read.delim(file.path(data_dir, "pr_tpm.tsv"),
                  fileEncoding = "UTF-8-BOM", check.names = FALSE,
                  stringsAsFactors = FALSE)
stopifnot("geneset000191922" %in% tpm$gene_id)
t110 <- as.numeric(tpm[tpm$gene_id == "geneset000191922",
                       c("CK-1", "CK-2", "CK-3",
                         "HD-1", "HD-2", "HD-3")])

d <- data.frame(sample = paste0(rep(c("CK", "HD"), each = 3), "-",
                                rep(1:3, 2)),
                group  = rep(c("Control", "Kelp farm"), each = 3),
                cov    = v27,
                tpm    = t110)
sc <- 4  # secondary-axis scale factor (TPM -> coverage range)
cat("bin.27 mean coverage CK =", round(mean(v27[1:3]), 2),
    " HD =", round(mean(v27[4:6]), 2), "\n")
cat("OTU110 mean TPM CK =", round(mean(t110[1:3]), 2),
    " HD =", round(mean(t110[4:6]), 2), "\n")

pD <- ggplot(d, aes(x = sample, fill = group)) +
  geom_col(aes(y = cov), width = 0.62, color = "grey25",
           linewidth = 0.15, show.legend = FALSE) +
  # black diamonds on a dotted line (kept dotted to stay distinct from the
  # dashed threshold lines of panels A/B)
  geom_point(aes(y = tpm * sc), shape = 18, size = 2.6,
             color = "black", show.legend = FALSE) +
  geom_line(aes(y = tpm * sc, group = 1), color = "black",
            linewidth = 0.55, linetype = "dotted") +
  scale_fill_manual(values = c("Control" = col_ck,
                               "Kelp farm" = col_hd)) +
  scale_y_continuous(name = "MAG coverage (bin.27, \u00d7)",
                     breaks = c(0, 20, 40),
                     sec.axis = sec_axis(
                       ~ . / sc, name = "PR-gene TPM (OTU110)",
                       breaks = c(0, 5, 10))) +
  coord_cartesian(ylim = c(0, 55), clip = "off") +
  labs(x = NULL,
       caption = "bars: MAG coverage; diamonds: PR-gene TPM") +
  theme(axis.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 9),
        plot.caption = element_text(size = 6.5, color = "grey35"))

ggsave(file.path(out_dir, "fig8d_bin27_otu110.tiff"), pD,
       width = 90, height = 62, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig8d_bin27_otu110.pdf"), pD,
       width = 90, height = 62, units = "mm", device = cairo_pdf)

# ============================================================================
# Step 8: assemble the final 2 x 2 figure with cowplot::plot_grid(align = "none")
#   Each panel is scaled as a whole into its grid cell - internal layouts are
#   not re-arranged, so a panel tuned at its single-panel size stays as
#   designed. (patchwork's default cross-panel alignment was measured to
#   compress panel A from 62 to 18 mm and D from 68 to 23 mm; freeing the
#   axes there failed too - hence the cowplot route.)
# ============================================================================
fig8 <- plot_grid(pA, pB, pC, pD, ncol = 2, align = "none",
                  labels = c("A", "B", "C", "D"),
                  label_size = 11, label_fontface = "bold")

print(fig8)

ggsave(file.path(out_dir, "Fig8_mag_rhodopsin.tiff"), fig8,
       width = 180, height = 130, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Fig8_mag_rhodopsin.pdf"), fig8,
       width = 180, height = 130, units = "mm", device = cairo_pdf)

cat("== Fig.8 done ->", out_dir, "==\n")
