# ============================================================================
# Fig.9  PR-CAZy coupling:
#        (A) sample x function z-scored heatmap with a habitat strip
#            (PR total TPM / 10 kelp-polysaccharide CAZy families / top 8
#            NR genera);
#        (B) Spearman rho bars, PR total TPM vs the same families
# ----------------------------------------------------------------------------
# Project : kelpPR - metagenomic study of kelp aquaculture effects on
#           bacterioplankton and proteorhodopsin-based phototrophy
# Author  : Fei Zhang
# Date    : 2026-10-27
# Version : 1.0
# Input   : data/pr_tpm.tsv (gene_id, CK-1..HD-3; column sums = PR total TPM),
#           Fig9/cazy_family_relabundance.tsv (Family, activities, CK-1..HD-3;
#             379 families),
#           Fig9/nr_genus_top.tsv (genus, phylum, CK-1..HD-3, ...; 40 rows,
#             duplicate genus names need aggregation),
#           Fig9/pr_cazy_coupling.tsv (family, spearman_prTpm_vs_family,
#             activity; 379 rows)
# Output  : Fig9/Fig9_pr_cazy_coupling.{tiff,pdf} (190 x 150 mm, 300 dpi + vector)
# Requires: ggplot2, patchwork
# ----------------------------------------------------------------------------
# fams: the 10 kelp-polysaccharide-degrading / laminarin-related families
#   (PL7 PL6 PL17 PL1 PL18 GH16 GH17 GH95 CBM16 CBM35) - alginate (PL
#   enzymes), laminarin (GH16/GH17, always kept - they anchor the
#   carbon-fixation plus polysaccharide narrative), fucoidan (GH95) and
#   binding modules (CBM16/CBM35).
# Panel A: three z-scored heatmap blocks side by side (pHab habitat strip |
#   p_pr single row | p_cz 10 rows | p_gn 8 rows), widths = c(0.07, 0.9, 9, 7)
#   so the CAZy block is the widest; z-scores standardised per row, fill
#   clipped to [-2.6, 2.6] so extreme rows do not flatten the colour scale;
#   duplicate genera are aggregated first (nr_genus_top.tsv repeats genus
#   names); when building z_long remember as.vector(matrix) reads
#   column-major - the rep(feature, times) / rep(sample, each) pairing order
#   must not be swapped.
# Panel B: Spearman rho horizontal bars for the same 10 families, positive
#   red / negative blue, values at the bar ends; the legend sits inside the
#   panel's empty band (right of the short GH16/GH17 bars); the caption
#   states n = 6 samples.
# Panel tags: A is attached to pHab (the leftmost sub-plot of panel A) and B
#   to pB, both flush with their left edges. plot_annotation(tag_levels)
#   cannot be used here - it would cascade into the nested sub-plots of
#   panel A and produce five tags A-E.
# ----------------------------------------------------------------------------

# ---- Step 1: package check ---------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
theme_set(theme_classic(base_size = 10))
samples <- c("CK-1", "CK-2", "CK-3", "HD-1", "HD-2", "HD-3")
grp <- c("CK-1" = "Control", "CK-2" = "Control", "CK-3" = "Control",
         "HD-1" = "Kelp farm", "HD-2" = "Kelp farm", "HD-3" = "Kelp farm")

# ---- Step 2: locate working directories ---------------------------------------
locate_dirs <- function(fig_id) {
  up <- if (dir.exists("data")) "." else if (dir.exists(file.path("..", "data"))) ".."
        else stop("Cannot find the project root. Set the working directory to ",
                  "the project root or to ", fig_id, "/")
  list(data_dir = file.path(up, "data"), out_dir = file.path(up, fig_id))
}
dirs <- locate_dirs("Fig9")
data_dir <- dirs$data_dir
out_dir  <- dirs$out_dir

# ---- Step 3: read data ----------------------------------------------------------
pr <- read.delim(file.path(data_dir, "pr_tpm.tsv"),
                 fileEncoding = "UTF-8-BOM", check.names = FALSE,
                 stringsAsFactors = FALSE)
pr_sum <- colSums(pr[, samples])

cazy <- read.delim(file.path(out_dir, "cazy_family_relabundance.tsv"),
                   fileEncoding = "UTF-8-BOM", check.names = FALSE,
                   stringsAsFactors = FALSE)
fams <- c("PL7", "PL6", "PL17", "PL1", "PL18", "GH16", "GH17",
          "GH95", "CBM16", "CBM35")
stopifnot(all(is.element(fams, cazy$Family)))
cazy_sel <- cazy[match(fams, cazy$Family), samples]
rownames(cazy_sel) <- fams

nr <- read.delim(file.path(out_dir, "nr_genus_top.tsv"),
                 fileEncoding = "UTF-8-BOM", check.names = FALSE,
                 stringsAsFactors = FALSE)
nr$genus <- gsub("_", " ", nr$genus)
nr_agg <- aggregate(nr[, samples], by = list(genus = nr$genus), sum)
nr_agg <- nr_agg[order(-rowMeans(nr_agg[, samples])), ]
top_genera <- head(nr_agg, 8)
gen_mat <- as.matrix(top_genera[, samples])
rownames(gen_mat) <- top_genera$genus

# ============================================================================
# Step 4: Panel A - sample x function heatmap
# ============================================================================
mat <- rbind(PR = pr_sum, cazy_sel, gen_mat)
z <- t(scale(t(mat)))
stopifnot(nrow(z) == 1 + length(fams) + nrow(gen_mat))
z_long <- data.frame(
  feature = rep(rownames(z), times = ncol(z)),
  sample = factor(rep(colnames(z), each = nrow(z)), levels = samples),
  z = as.vector(z),
  stringsAsFactors = FALSE)
z_long$block <- ifelse(z_long$feature == "PR", "PR total TPM",
                ifelse(is.element(z_long$feature, fams), "CAZy family",
                       "Genus (NR)"))
z_long$block <- factor(z_long$block,
                       levels = c("PR total TPM", "CAZy family", "Genus (NR)"))
z_long$feature <- factor(z_long$feature, levels = rev(rownames(z)))

hab <- data.frame(sample = factor(samples, levels = samples),
                  Habitat = factor(grp[samples],
                                   levels = c("Control", "Kelp farm")))

mk_tile <- function(df, head_txt, tag_txt = NULL, keep_legend = FALSE) {
  ggplot(df, aes(x = feature, y = sample, fill = z)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-2.6, 2.6),
                         name = expression(italic(z)*"-score")) +
    labs(x = NULL, y = NULL, title = head_txt, tag = tag_txt) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5,
                                     color = "black"),
          axis.text.y = element_text(size = 8, color = "black"),
          plot.title = element_text(size = 8, hjust = 0.5,
                                    margin = margin(b = 3)),
          plot.tag = element_text(size = 11, face = "bold"),
          legend.position = if (keep_legend) "right" else "none",
          legend.key.height = unit(0.7, "cm"),
          legend.key.width = unit(0.3, "cm"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8),
          plot.margin = margin(2, 2, 2, 2))
}

d_pr <- z_long[z_long$block == "PR total TPM", ]
d_cz <- z_long[z_long$block == "CAZy family", ]
d_cz$feature <- factor(as.character(d_cz$feature), levels = rev(fams))
d_gn <- z_long[z_long$block == "Genus (NR)", ]
d_gn$feature <- factor(as.character(d_gn$feature),
                       levels = rev(rownames(gen_mat)))

pHab <- ggplot(hab, aes(x = "Habitat", y = sample, fill = Habitat)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("Control" = "#4C72B0", "Kelp farm" = "#C44E52"),
                    name = "Habitat") +
  labs(x = NULL, y = NULL, title = "Habitat", tag = "A") +
  theme_void() +
  theme(axis.text.y = element_blank(), legend.position = "none",
        plot.title = element_text(size = 8, hjust = 0.5,
                                  margin = margin(b = 3, r = 16)),
        plot.tag = element_text(size = 11, face = "bold"),
        plot.margin = margin(2, 2, 2, 2))

p_pr <- mk_tile(d_pr, "PR total TPM")
p_cz <- mk_tile(d_cz, "CAZy family")

# genus labels italic (nomenclature); "Candidatus " prefixes stay roman,
# assembled with plotmath expressions
gn_lvls <- levels(d_gn$feature)
gn_lbls <- do.call(c, lapply(gn_lvls, function(g)
  if (startsWith(g, "Candidatus "))
    parse(text = sprintf('"Candidatus "*italic(%s)',
                         sub("^Candidatus ", "", g)))
  else
    parse(text = sprintf("italic(%s)", g))))
p_gn <- mk_tile(d_gn, "Genus (NR)", keep_legend = TRUE) +
  scale_x_discrete(labels = gn_lbls)

pA <- pHab + p_pr + p_cz + p_gn +
  plot_layout(widths = c(0.07, 0.9, 9, 7), guides = "collect")

# ============================================================================
# Step 5: Panel B - Spearman PR vs kelp-polysaccharide families
# ============================================================================
coup <- read.delim(file.path(out_dir, "pr_cazy_coupling.tsv"),
                   fileEncoding = "UTF-8-BOM", check.names = FALSE,
                   stringsAsFactors = FALSE)
topn <- coup[is.element(coup$family, fams),
             c("family", "spearman_prTpm_vs_family")]
stopifnot(nrow(topn) == length(fams))
names(topn)[2] <- "sig"
topn <- topn[order(topn$sig), ]
topn$family <- factor(topn$family, levels = topn$family)

pB <- ggplot(topn, aes(x = sig, y = family, fill = sig > 0)) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.35) +
  geom_col(width = 0.7, color = "grey25", linewidth = 0.15) +
  geom_text(aes(label = sprintf("%.2f", sig)),
            hjust = ifelse(topn$sig >= 0, -0.25, 1.25), size = 2.4) +
  scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"),
                    labels = c("FALSE" = "Negative", "TRUE" = "Positive"),
                    name = "Spearman rho") +
  scale_x_continuous(limits = c(-1.15, 1.15),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = expression("Spearman's "*rho*", PR total TPM vs family abundance"),
       y = "CAZy family",
       caption = "Kelp polysaccharide-degrading / laminarin-related families (n = 6 samples)",
       tag = "B") +
  # legend anchored in the empty band right of the short GH16/GH17 bars,
  # clear of the long PL18/CBM35 bars
  theme(axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 9),
        legend.position = c(0.85, 0.62),
        legend.background = element_rect(fill = "white", color = "grey85",
                                         linewidth = 0.2),
        legend.key.size = unit(0.32, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7.5),
        legend.margin = margin(2, 3, 2, 3),
        plot.tag = element_text(size = 11, face = "bold"),
        plot.caption = element_text(size = 6.5, color = "grey35"))

# ---- Step 6: assemble and export -------------------------------------------------
fig9 <- pA / pB +
  plot_layout(heights = c(1.1, 1))

print(fig9)

ggsave(file.path(out_dir, "Fig9_pr_cazy_coupling.tiff"), fig9,
       width = 190, height = 150, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Fig9_pr_cazy_coupling.pdf"), fig9,
       width = 190, height = 150, units = "mm", device = cairo_pdf)

cat("== Fig.9 done ->", file.path(out_dir, "Fig9_pr_cazy_coupling.pdf"), "==\n")
cat("Spearman rho (10 fams):\n"); print(topn)
