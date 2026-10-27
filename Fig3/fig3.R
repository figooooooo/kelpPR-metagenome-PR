# ============================================================================
# Fig.3  Alpha diversity (A, five-index boxplots, 190 mm wide) and beta
#        diversity (B, Bray-Curtis PCoA with exact PERMANOVA), assembled
#        into a two-row figure
# ----------------------------------------------------------------------------
# Project : kelpPR - metagenomic study of kelp aquaculture effects on
#           bacterioplankton and proteorhodopsin-based phototrophy
# Author  : Fei Zhang
# Date    : 2026-10-27
# Version : 1.0
# Input   : Fig3/alpha_index_stat.tsv   (6 samples x 7 alpha indices)
#           Fig3/pcoa_bray_coordinate.tsv, Fig3/pcoa_bray_explained.tsv
#           (vendor-delivered Bray-Curtis PCoA solution, used in the main
#           panel to match the manuscript's data-source wording)
#           data/group.txt, data/species_relabund.tsv (species-level relative
#           abundance: PERMANOVA distance input)
# Output  : Fig3/Fig3A_alpha_box.{tiff,pdf}, Fig3/Fig3B_pcoa.{tiff,pdf},
#           Fig3/Fig3_final.{tiff,pdf} (190 mm, 300 dpi + vector), plus
#           exact PERMANOVA results on the console
# Requires: ggplot2, dplyr, tidyr, vegan, patchwork
# ----------------------------------------------------------------------------

# ---- Step 1: package check (install missing packages from CRAN mirror) -----
pkgs <- c("ggplot2", "dplyr", "tidyr", "vegan", "patchwork")
need <- pkgs[!is.element(pkgs, rownames(installed.packages()))]
if (length(need) > 0)
  install.packages(need, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN",
                   type = "binary", quiet = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Step 2: locate working directories ------------------------------------
# Run from the project root (or from Fig3/: the script walks one level up).
# Figure-specific inputs live in Fig3/, shared tables in data/.
locate_dirs <- function(fig_id) {
  up <- if (dir.exists("data")) "." else if (dir.exists(file.path("..", "data"))) ".."
        else stop("Cannot find the project root. Set the working directory to ",
                  "the project root or to ", fig_id, "/")
  list(data_dir = file.path(up, "data"), out_dir = file.path(up, fig_id))
}
dirs <- locate_dirs("Fig3")
data_dir <- dirs$data_dir
out_dir  <- dirs$out_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Step 3: Panel A - alpha diversity boxplots (ACE / Chao1 / Shannon / Simpson
#   / Pielou). Wide table -> long table -> facet_wrap (small multiples);
#   exact Mann-Whitney per index (n = 3 vs 3); p labels drawn manually with
#   geom_text (no ggpubr dependency)
# ============================================================================

# ---- 3.1 read data -----------------------------------------------------------
alpha <- read.delim(file.path(out_dir, "alpha_index_stat.tsv"),
                    header = TRUE, check.names = FALSE)
grp <- read.delim(file.path(data_dir, "group.txt"),
                  header = TRUE, check.names = FALSE)
stopifnot(nrow(alpha) == 6, all(is.element(alpha$sample, grp$sample)))

# ---- 3.2 cross-check against the manuscript (observed_species assertion) ----
obs_ck <- range(alpha$observed_species[grepl("^CK", alpha$sample)])
obs_hd <- range(alpha$observed_species[grepl("^HD", alpha$sample)])
cat("observed_species  CK:", obs_ck[1], "-", obs_ck[2],
    " HD:", obs_hd[1], "-", obs_hd[2], "\n")
stopifnot(obs_ck[1] == 15577, obs_ck[2] == 16150,     # manuscript: CK 15,577-16,150
          obs_hd[1] == 14868, obs_hd[2] == 15055)     # manuscript: HD 14,868-15,055

# ---- 3.3 wide -> long ---------------------------------------------------------
# "five indices" per the figure caption: ACE, Chao1 (richness), Shannon,
# Simpson (diversity), Pielou evenness; observed_species is given in the text
idx_keep <- c("Chao1", "ACE", "Shannon", "Simpson", "pielou_evenness")
idx_label <- c(Chao1 = "Chao1", ACE = "ACE", Shannon = "Shannon",
               Simpson = "Simpson", pielou_evenness = "Pielou evenness")

long <- alpha |>
  dplyr::select(sample, dplyr::all_of(idx_keep)) |>
  pivot_longer(cols = -sample, names_to = "index", values_to = "value") |>
  left_join(grp, by = "sample") |>
  mutate(index = factor(index, levels = idx_keep),
         treat1 = factor(treat1, levels = c("CK", "HD")))

# ---- 3.4 exact Mann-Whitney per index -------------------------------------------
pvals <- sapply(idx_keep, function(k) {
  wilcox.test(value ~ treat1, data = long[long$index == k, ],
              exact = TRUE)$p.value
})
p_lab <- ifelse(pvals < 0.001, "p < 0.001", sprintf("p = %.3f", pvals))
cat("\n== exact Mann-Whitney, five indices ==\n")
print(round(pvals, 3))
pval_df <- data.frame(index = factor(idx_keep, levels = idx_keep), p_lab = p_lab)

# ---- 3.5 point alignment + tiny fixed offsets for near-coincident points ------
# Rules (fixed, transparent, reproducible, no random seed):
#   within a facet and group, after sorting by y, adjacent points closer than
#   5.66% of the facet data range (= 24 px at 300 dpi: point diameter
#   1.6 mm ~ 19 px + gap; 513 px panel height x 1.21 y-expand) form a
#   near-coincident cluster; 2-point cluster -> offsets -0.10/+0.10, 3-point
#   -> -0.10/0/+0.10; singletons are not offset. In this dataset only
#   Chao1 CK-1/CK-3 (delta = 18.7) and ACE CK-1/CK-3 (delta = 21.5) trigger,
#   2 pairs / 4 points; the other 26 points sit exactly on the group centres.
#   Companion caption sentence: Points with near-identical values within a
#   group were horizontally offset by 0.10 x-axis units for visibility.
OFFSET_W <- 0.10   # horizontal offset for near-coincident points (x units)
VIS_THR  <- 0.0566 # near-coincidence threshold (fraction of facet data range)

long <- long |>
  dplyr::group_by(index) |>
  dplyr::mutate(facet_rng = max(value) - min(value)) |>
  dplyr::ungroup() |>
  dplyr::arrange(index, treat1, value) |>
  dplyr::group_by(index, treat1) |>
  dplyr::mutate(
    new_cl = c(TRUE, diff(value) >= VIS_THR * facet_rng[1]),  # TRUE = new cluster
    cl = cumsum(new_cl)
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(index, treat1, cl) |>
  dplyr::mutate(
    k = dplyr::n(),
    x_off = if (k[1] == 1) 0
            else (dplyr::row_number() - (k[1] + 1) / 2) * (2 * OFFSET_W / (k[1] - 1))
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(sample)   # restore sample order for checking

n_off <- sum(long$x_off != 0)
cat("\n== near-coincident point offsets ==\n")
print(long[long$x_off != 0, c("sample", "index", "value", "x_off")])
cat("offset points:", n_off, "(expected 4: Chao1/ACE, CK-1/CK-3 each)\n")
stopifnot(n_off == 4)  # guards against data drift: only the two known pairs

# continuous x coordinates (CK = 1, HD = 2) - equivalent to the discrete
# layout (expand add 0.6 reproduces the ggplot2 discrete default)
long$x_pos <- as.numeric(long$treat1)

# ---- 3.6 plot --------------------------------------------------------------------
pA <- ggplot(long, aes(x = x_pos, y = value, fill = treat1)) +
  geom_boxplot(width = 0.5, colour = "grey25", linewidth = 0.35,
               outlier.shape = NA, alpha = 0.85) +
  # n = 3 shown per point, aligned in columns: non-coincident points sit
  # exactly on the group centre; near-coincident pairs get the fixed
  # +/-0.10 offset from 3.5; point size 1.6 mm
  geom_point(aes(x = x_pos + x_off, colour = treat1),
             size = 1.6, show.legend = FALSE) +
  geom_text(data = pval_df, aes(label = p_lab),
            x = 1.5, y = Inf, vjust = 1.5, size = 2.9,
            inherit.aes = FALSE) +
  # top headroom 16% + vjust 1.5 keeps the p label clear of the highest
  # box / point under free scales
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
  # continuous x axis (points need numeric offsets), add = 0.6 reproduces
  # discrete-layout panel margins
  scale_x_continuous(breaks = c(1, 2), labels = c("CK", "HD"),
                     expand = expansion(add = 0.6)) +
  facet_wrap(~index, scales = "free", nrow = 1,
             labeller = labeller(index = idx_label)) +
  scale_fill_manual(values = c(CK = "#4C72B0", HD = "#C44E52")) +
  scale_colour_manual(values = c(CK = "#2E4A78", HD = "#8F3336")) +  # darker than boxes
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 9.5, face = "bold"),
        axis.text.x = element_text(size = 8.5))

# ---- 3.7 export panel (190 mm five-panel strip) -----------------------------------
ggsave(file.path(out_dir, "Fig3A_alpha_box.tiff"), pA,
       width = 190, height = 55, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "Fig3A_alpha_box.pdf"), pA,
       width = 190, height = 55, units = "mm")

cat("\n== Fig.3 panel A done ==\nOutput:", out_dir, "/Fig3A_alpha_box.{tiff,pdf}\n")

# ============================================================================
# Step 4: Panel B - beta diversity PCoA (Bray-Curtis) + exact PERMANOVA.
#   No random permutations: with n = 3 vs 3 only C(6,3) = 20 groupings exist,
#   so the test enumerates all of them (exact), matching the exact
#   Mann-Whitney statistics used elsewhere; complete separation gives
#   p = 2/20 = 0.100, the power floor of this design.
#   Three cross-checks guard against silent bugs: (1) manual pseudo-F vs
#   vegan::adonis2; (2) self-computed cmdscale coordinates vs the
#   vendor-delivered table (|r| = 1, PCoA axis signs are arbitrary);
#   (3) Bray-Curtis as the primary metric + binary Jaccard as a check.
# ============================================================================

# ---- 4.1 read data ---------------------------------------------------------------
grp <- read.delim(file.path(data_dir, "group.txt"),
                  check.names = FALSE, stringsAsFactors = FALSE)
sp <- read.delim(file.path(data_dir, "species_relabund.tsv"),
                 row.names = 1, check.names = FALSE)
sp <- as.data.frame(t(sp))          # rows = samples x columns = species (vegan)

coord_bray <- read.delim(file.path(out_dir, "pcoa_bray_coordinate.tsv"),
                         row.names = 1, check.names = FALSE)
exp_bray   <- read.delim(file.path(out_dir, "pcoa_bray_explained.tsv"),
                         check.names = FALSE)
# consistency assertion: sample order must match, misalignment would be hard
# to spot downstream
stopifnot(identical(rownames(sp), grp$sample),
          identical(rownames(coord_bray), grp$sample))
group <- factor(grp$treat1, levels = c("CK", "HD"))

# ---- 4.2 distance matrices (same source as the vendor run: species-level) --------
d_bray <- vegdist(sp, method = "bray")                    # primary metric
d_jac  <- vegdist(sp, method = "jaccard", binary = TRUE)  # check metric

# ---- 4.3 PERMANOVA: complete enumeration exact test ------------------------------
# pseudo-F decomposition (Anderson 2001): SS_T = (1/n) sum d_ij^2;
# SS_W = sum_g (1/n_g) sum d_ij^2; F = (SS_A/(g-1)) / (SS_W/(n-g));
# R2 = SS_A/SS_T. Notes: both sum-of-squares terms must square d; the return
# must contain F, R2 and p (sprintf silently drops NULLs in annotations)
pseudo_f <- function(d, grp) {
  n  <- length(grp)
  dm <- as.matrix(d)
  ss_t <- sum(dm[upper.tri(dm)]^2) / n
  ss_w <- 0
  for (g in unique(grp)) {
    idx <- which(grp == g)
    if (length(idx) > 1) {
      sub <- dm[idx, idx]
      ss_w <- ss_w + sum(sub[upper.tri(sub)]^2) / length(idx)
    }
  }
  ss_a <- ss_t - ss_w
  g <- length(unique(grp))
  list(F  = (ss_a / (g - 1)) / (ss_w / (n - g)),
       R2 = ss_a / ss_t)
}
# all 20 groupings (C(6,3)), F computed for each - exact, no sampling
perm_exact <- function(d, labels) {
  n <- length(labels)
  combos <- combn(n, n / 2)
  obs   <- pseudo_f(d, labels)
  f_obs <- obs$F
  f_perm <- apply(combos, 2, function(idx) {
    lab <- rep("CK", n)
    lab[is.element(seq_len(n), idx)] <- "HD"
    pseudo_f(d, lab)$F
  })
  list(F = f_obs, R2 = obs$R2,
       p = sum(f_perm >= f_obs - 1e-12) / ncol(combos),
       n_perm = ncol(combos))
}

set.seed(42)
res_bray <- perm_exact(d_bray, group)
res_jac  <- perm_exact(d_jac,  group)
chk_bray <- adonis2(d_bray ~ group, permutations = 199)

cat("\n== PERMANOVA (exact enumeration, C(6,3) = 20 groupings) ==\n")
cat(sprintf("Bray-Curtis    : F = %.3f, R2 = %.4f, p = %.3f (n_perm = %d)\n",
            res_bray$F, res_bray$R2, res_bray$p, res_bray$n_perm))
cat(sprintf("binary Jaccard : F = %.3f, R2 = %.4f, p = %.3f (n_perm = %d)\n",
            res_jac$F, res_jac$R2, res_jac$p, res_jac$n_perm))
cat(sprintf("cross-check adonis2: F = %.3f, R2 = %.4f (must match manual Bray-Curtis)\n",
            unname(chk_bray[1, "F"]), unname(chk_bray[1, "R2"])))
stopifnot(abs(unname(chk_bray[1, "F"]) - res_bray$F) < 1e-6,
          abs(res_bray$p - 0.100) < 1e-9)

# ---- 4.4 PCoA coordinate cross-check (cmdscale vs vendor table) -------------------
# (1) PCoA axis signs are arbitrary, unify signs before comparing;
# (2) cmdscale needs eig = TRUE for the explained fractions
pco_own <- cmdscale(d_bray, k = 5, eig = TRUE)
cor_pc1 <- cor(pco_own$points[, 1], coord_bray[, 1])
cor_pc1 <- sign(cor_pc1) * abs(cor_pc1)
cat(sprintf("self-computed vs vendor PC1 coordinates: |r| = %.6f; PC1 explains %.2f%%\n",
            abs(cor_pc1), exp_bray$explained[1]))
stopifnot(abs(abs(cor_pc1) - 1) < 1e-4)

# ---- 4.5 plot (main panel uses the vendor coordinates, matching the caption) ------
df <- data.frame(sample = rownames(coord_bray),
                 group  = group,
                 PC1 = coord_bray[, 1],
                 PC2 = coord_bray[, 2])
lab1 <- sprintf("PC1 (%.1f%%)", exp_bray$explained[1])
lab2 <- sprintf("PC2 (%.1f%%)", exp_bray$explained[2])

# legend in the top-right corner (top-left holds the PERMANOVA annotation:
# together they form a top information band with zero point overlap)
pB <- ggplot(df, aes(x = PC1, y = PC2, colour = group)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey70", linewidth = 0.3) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sample), size = 2.1, show.legend = FALSE,
            vjust = 1.9, colour = "grey30") +
  annotate("text", x = min(df$PC1), y = max(df$PC2) + 0.004, hjust = 0, vjust = 1,
           label = sprintf("PERMANOVA:  R\u00b2 = %.2f,  p = %.3f",
                           res_bray$R2, res_bray$p),
           size = 2.6, colour = "grey20") +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.18))) +
  scale_colour_manual(values = c(CK = "#4C72B0", HD = "#C44E52"), name = NULL) +
  labs(x = lab1, y = lab2) +
  theme_classic(base_size = 8) +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        legend.key.size = unit(0.25, "cm"), legend.text = element_text(size = 7))

ggsave(file.path(out_dir, "Fig3B_pcoa.tiff"), pB,
       width = 90, height = 85, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Fig3B_pcoa.pdf"), pB,
       width = 90, height = 85, units = "mm")

cat("\n== Fig.3 panel B done ==\nOutput:", out_dir, "/Fig3B_pcoa.{tiff,pdf}\n")
cat(sprintf("values for the manuscript: PERMANOVA (Bray-Curtis) R2 = %.2f, F = %.2f, p = 0.100 (exact)\n",
            res_bray$R2, res_bray$F))

# ============================================================================
# Step 5: assemble the final figure
#   A is a 190 mm five-panel strip, B a 90 mm square scatter. Stacked directly,
#   B would be stretched - plot_spacer() pads both sides of row 2 so B keeps
#   its aspect ratio and stays centred.
# ============================================================================

fig3 <- pA /
        (plot_spacer() + pB + plot_spacer() +
           plot_layout(widths = c(0.5, 1, 0.5))) +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold"))

print(fig3)

ggsave(file.path(out_dir, "Fig3_final.tiff"), fig3,
       width = 190, height = 160, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Fig3_final.pdf"), fig3,
       width = 190, height = 160, units = "mm", device = cairo_pdf)

cat("\n== Fig.3 final figure done ==\nOutput:", out_dir, "/Fig3_final.{tiff,pdf}\n")
