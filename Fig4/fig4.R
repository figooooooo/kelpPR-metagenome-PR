# ============================================================================
# Fig.4  LEfSe analysis of zone-discriminating taxa:
#        (A) circular cladogram of the differential-lineage backbone,
#        (B) two-directional LDA score bar plot (HD right / CK left)
# ----------------------------------------------------------------------------
# Project : kelpPR - metagenomic study of kelp aquaculture effects on
#           bacterioplankton and proteorhodopsin-based phototrophy
# Author  : Fei Zhang
# Date    : 2026-10-27
# Version : 1.0
# Input   : data/lefse_rep.tsv (full LEfSe feature table),
#           data/{phylum,class,order,family,genus,species}_relabund.tsv
#           (used for the abundance-direction check)
# Output  : Fig4/Fig4A_cladogram.{tiff,pdf}, Fig4/Fig4B_lda.{tiff,pdf},
#           Fig4/Fig4_final.{tiff,pdf} (190 mm, 300 dpi + vector)
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

# ---- Step 2: locate working directories ------------------------------------
locate_dirs <- function(fig_id) {
  up <- if (dir.exists("data")) "." else if (dir.exists(file.path("..", "data"))) ".."
        else stop("Cannot find the project root. Set the working directory to ",
                  "the project root or to ", fig_id, "/")
  list(data_dir = file.path(up, "data"), out_dir = file.path(up, fig_id))
}
dirs <- locate_dirs("Fig4")
data_dir <- dirs$data_dir
out_dir  <- dirs$out_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

samples <- c("CK-1", "CK-2", "CK-3", "HD-1", "HD-2", "HD-3")
COL_HD <- "#C44E52"; COL_CK <- "#4C72B0"

# ============================================================================
# Step 3: LEfSe differential-taxon cleaning (identical rules as the Fig.2C
#   heatmap, so both figures show the same 21-taxa set):
#   keep labelled rows only; drop kingdom rank; deduplicate (keep max LDA);
#   verify each LEfSe label against the measured abundance direction and
#   remove artefacts (e.g. p__Proteobacteria labelled both HD and CK).
#   Numbering: sorted by (Groups == "HD", LDA), numbers run 1-21 top-down
#   (largest HD LDA = 1, smallest CK = 21); panel A shows the numbers on the
#   ring, panel B carries the full names prefixed with the same numbers.
# ============================================================================
lefse_raw <- read.delim(file.path(data_dir, "lefse_rep.tsv"),
                        header = TRUE, check.names = FALSE, colClasses = "character")
lefse <- lefse_raw[is.element(lefse_raw$Groups, c("HD", "CK")), ]
lefse$LDA <- as.numeric(lefse$LDA_value)
tok <- strsplit(lefse$Biomarker_names, ".", fixed = TRUE)
lefse$rank  <- vapply(tok, function(x) substr(x[length(x)], 1, 1), "")
lefse$label <- vapply(tok, function(x) substr(x[length(x)], 4, 100), "")
lefse <- lefse[lefse$rank != "k", ]
lefse <- lefse[order(lefse$Biomarker_names, -lefse$LDA), ]
lefse <- lefse[!duplicated(lefse$Biomarker_names), ]

# direction check: verify LEfSe labels against measured abundances,
# drop artefact rows
rank_files <- c(p = "phylum_relabund.tsv", c = "class_relabund.tsv",
                o = "order_relabund.tsv",  f = "family_relabund.tsv",
                g = "genus_relabund.tsv",  s = "species_relabund.tsv")
match_row <- function(tab, name) {
  if (is.element(name, rownames(tab))) return(name)
  alt <- gsub("__", "._", name)
  if (is.element(alt, rownames(tab))) return(alt)
  NA_character_
}
dirs_ck <- vapply(seq_len(nrow(lefse)), function(i) {
  tab  <- read.delim(file.path(data_dir, rank_files[[lefse$rank[i]]]),
                     header = TRUE, row.names = 1, check.names = FALSE)
  nm   <- match_row(tab, lefse$label[i])
  if (is.na(nm)) return(lefse$Groups[i])
  v    <- as.numeric(tab[nm, samples])
  ifelse(mean(v[4:6]) > mean(v[1:3]), "HD", "CK")
}, "")
drop <- dirs_ck != lefse$Groups
for (i in which(drop))
  cat(sprintf("  direction conflict removed: %s (LEfSe label %s)\n",
              lefse$label[i], lefse$Groups[i]))
lefse <- lefse[!drop, ]
stopifnot(nrow(lefse) == 21)                     # consistency with Fig.2C
cat("final differential taxa:", nrow(lefse), "(HD", sum(lefse$Groups == "HD"),
    "/ CK", sum(lefse$Groups == "CK"), ")\n")

# numbering shared by panels A and B (same sort key, numbers top-down 1-21)
lefse$no <- NA_integer_
lefse$no[order(lefse$Groups == "HD", lefse$LDA)] <- rev(seq_len(nrow(lefse)))
stopifnot(identical(sort(lefse$no), seq_len(nrow(lefse))))
cat("\n== numbering map (must match row-by-row between panels A and B) ==\n")
for (i in order(lefse$no))
  cat(sprintf("  %2d. [%s] %s\n", lefse$no[i], lefse$Groups[i], lefse$label[i]))

# ============================================================================
# Step 4: Panel A - circular cladogram of the differential-lineage backbone.
#   The full LEfSe table holds ~28k lineages (unreadable as one cladogram),
#   so the tree is the union of the 21 differential lineages' paths
#   (~70 nodes). Layout: textbook equal-angle fan (no ggtree needed) - each
#   leaf gets 2*pi/n of angle, internal nodes sit at the middle of their
#   subtree's share, radii are layered by rank (k/p/c/o/f/g/s); edges are a
#   parent arc plus a radial line. Labels are spread with an order-preserving
#   minimum angular distance, then re-centred.
# ============================================================================

# ---- 4.1 build the tree: node union of the 21 lineage paths ----------------
paths <- strsplit(lefse$Biomarker_names, ".", fixed = TRUE)
node_set <- unique(unlist(lapply(paths, function(x)
  vapply(seq_along(x), function(i) paste(x[1:i], collapse = "."), ""))))
node_rank <- vapply(strsplit(node_set, ".", fixed = TRUE),
                    function(x) substr(x[length(x)], 1, 1), "")
node_name <- vapply(strsplit(node_set, ".", fixed = TRUE),
                    function(x) substr(x[length(x)], 4, 100), "")
parent_of <- setNames(
  vapply(strsplit(node_set, ".", fixed = TRUE), function(x)
    if (length(x) == 1) NA_character_ else paste(x[-length(x)], collapse = "."), ""),
  node_set)

# ---- 4.2 circular layout: equal-angle fan ------------------------------------
RAD <- c(k = 0.10, p = 0.27, c = 0.44, o = 0.61, f = 0.78, g = 0.95, s = 1.12)
children_of <- function(nd) node_set[parent_of == nd & !is.na(parent_of)]
angle <- setNames(rep(NA_real_, length(node_set)), node_set)
assign_angles <- function(nd, a0, a1) {
  kids <- sort(children_of(nd))
  if (length(kids) == 0) { angle[nd] <<- (a0 + a1) / 2; return(invisible()) }
  span <- (a1 - a0) / length(kids)
  for (j in seq_along(kids))
    assign_angles(kids[j], a0 + (j - 1) * span, a0 + j * span)
  angle[nd] <<- (a0 + a1) / 2
}
assign_angles(node_set[1], -pi / 2, 3 * pi / 2)   # root starts at the bottom, clockwise
stopifnot(!any(is.na(angle)))

nd_rank <- node_rank[match(names(angle), node_set)]
xy <- data.frame(node = names(angle), rank = nd_rank,
                 name = node_name[match(names(angle), node_set)],
                 theta = as.numeric(angle), r = RAD[nd_rank],
                 stringsAsFactors = FALSE)
xy$x <- xy$r * cos(xy$theta); xy$y <- xy$r * sin(xy$theta)

# ---- 4.3 edge geometry (parent arc + radial line) ------------------------------
edge_draw <- list()
for (nd in node_set[-1]) {
  pa  <- parent_of[nd]
  th_c <- xy$theta[xy$node == nd]; th_p <- xy$theta[xy$node == pa]
  r_p  <- xy$r[xy$node == pa];    r_c  <- xy$r[xy$node == nd]
  arc  <- seq(th_p, th_c, length.out = 24)
  edge_draw[[length(edge_draw) + 1]] <- data.frame(
    x = r_p * cos(arc), y = r_p * sin(arc), node = nd, stringsAsFactors = FALSE)
  edge_draw[[length(edge_draw) + 1]] <- data.frame(
    x = c(r_p * cos(th_c), r_c * cos(th_c)),
    y = c(r_p * sin(th_c), r_c * sin(th_c)), node = nd, stringsAsFactors = FALSE)
}
edges <- do.call(rbind, edge_draw)
edges$side <- lefse$Groups[match(edges$node, lefse$Biomarker_names)]
edges$side[is.na(edges$side)] <- "bg"

# ---- 4.4 node / label data frames -----------------------------------------------
nodes <- xy
nodes$side <- lefse$Groups[match(nodes$node, lefse$Biomarker_names)]
nodes$LDA  <- lefse$LDA[match(nodes$node, lefse$Biomarker_names)]
nodes$is_diff <- !is.na(nodes$side)

# display names (console cross-check only; on-figure labels are numbers,
# full names are carried by panel B)
rank_full <- c(p = "phylum", c = "class", o = "order", f = "family",
               g = "genus", s = "species", k = "kingdom")
disp <- ifelse(nodes$name == "Unclassified",
               paste0("Unclassified (", rank_full[nodes$rank], ")"), nodes$name)
disp <- gsub("sp__", "sp. ", disp); disp <- gsub("_", " ", disp)
# outer-ring label layout: slots start at node angles, overlapping labels are
# pushed apart along the circle (order-preserving minimum-distance
# constraint), then the whole set is re-centred to cancel the one-sided
# cumulative drift. Numbers are narrow (~0.1 units << DELTA_MIN 0.26 rad),
# so spacing is comfortable
lab_df <- nodes[nodes$is_diff, ]
lab_df$disp <- disp[nodes$is_diff]
lab_df$no <- lefse$no[match(lab_df$node, lefse$Biomarker_names)]
stopifnot(!any(is.na(lab_df$no)))
lab_df <- lab_df[order(lab_df$theta), ]
n_lab <- nrow(lab_df)
slot <- lab_df$theta
for (i in 2:n_lab) while (slot[i] < slot[i - 1] - 1e-9) slot[i] <- slot[i] + 2 * pi
DELTA_MIN <- 0.26   # minimum angular distance between adjacent labels (rad, ~14.9 deg)
for (pass in 1:3) {
  for (i in 2:n_lab)
    if (slot[i] - slot[i - 1] < DELTA_MIN) slot[i] <- slot[i - 1] + DELTA_MIN
}
slot <- slot - mean(slot - lab_df$theta)   # re-centre
lab_df$slot <- slot
lab_df$slot_wr <- atan2(sin(slot), cos(slot))
R_TXT <- 1.24   # label ring radius, widened together with the label font size
lab_df$lab_x <- R_TXT * cos(lab_df$slot_wr); lab_df$lab_y <- R_TXT * sin(lab_df$slot_wr)
lab_df$hjust <- ifelse(cos(lab_df$slot_wr) >= 0, 0, 1)

# ---- 4.5 plot ----------------------------------------------------------------------
ring_df <- data.frame(r = RAD[c("p", "c", "o", "f", "g", "s")])
rank_lbl <- data.frame(r = RAD[c("p", "c", "o", "f", "g", "s")],
                       lab = c("phylum", "class", "order", "family", "genus", "species"))
th_axis <- pi / 4
rank_lbl$x <- (rank_lbl$r + 0.055) * cos(th_axis)
rank_lbl$y <- (rank_lbl$r + 0.055) * sin(th_axis)

pA <- ggplot() +
  geom_path(data = do.call(rbind, lapply(ring_df$r, function(rr) {
    th <- seq(0, 2 * pi, length.out = 181)
    data.frame(x = rr * cos(th), y = rr * sin(th), grp = rr)
  })), aes(x, y, group = grp), colour = "grey88", linewidth = 0.25) +
  geom_path(data = edges[edges$side == "bg", ], aes(x, y, group = node),
            colour = "grey70", linewidth = 0.4) +
  geom_path(data = edges[edges$side == "HD", ], aes(x, y, group = node),
            colour = COL_HD, linewidth = 0.7, alpha = 0.85) +
  geom_path(data = edges[edges$side == "CK", ], aes(x, y, group = node),
            colour = COL_CK, linewidth = 0.7, alpha = 0.85) +
  geom_point(data = nodes[!nodes$is_diff, ], aes(x, y), size = 0.7, colour = "grey60") +
  geom_point(data = nodes[nodes$is_diff, ],
             aes(x, y, size = LDA, fill = side), shape = 21,
             colour = "grey25", stroke = 0.25) +
  scale_fill_manual(values = c(HD = COL_HD, CK = COL_CK), guide = "none") +
  scale_size_continuous(range = c(2.2, 5.2), guide = "none") +
  # corner notes removed: "dot size proportional to LDA" / "numbered taxa are
  # listed in panel B" live in the manuscript caption, canvas narrowed
  geom_segment(data = lab_df,
               aes(x = r * cos(theta), y = r * sin(theta),
                   xend = (R_TXT - 0.02) * cos(slot_wr), yend = (R_TXT - 0.02) * sin(slot_wr)),
               colour = "grey55", linewidth = 0.25) +
  geom_text(data = lab_df, aes(x = lab_x, y = lab_y, label = no,
                               colour = side, hjust = hjust),
            parse = FALSE, size = 4.8, show.legend = FALSE) +
  scale_colour_manual(values = c(HD = COL_HD, CK = COL_CK)) +
  geom_text(data = rank_lbl, aes(x, y, label = lab), size = 3.5,
            colour = "grey45", hjust = 0.5, vjust = -0.4) +
  coord_fixed(xlim = c(-1.42, 1.42), ylim = c(-1.32, 1.32)) +
  theme_void() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "Fig4A_cladogram.tiff"), pA,
       width = 150, height = 140, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Fig4A_cladogram.pdf"), pA,
       width = 150, height = 140, units = "mm", device = cairo_pdf)

cat("\n== Fig.4 panel A done ==\n")
cat("total nodes:", nrow(nodes), "(differential", sum(nodes$is_diff), ")\n")
cat("output:", out_dir, "/Fig4A_cladogram.{tiff,pdf}\n")

# ============================================================================
# Step 5: Panel B - two-directional LDA bar plot (HD right / CK left, LDA
#   decreasing within blocks; y order matches the Fig.2C heatmap rows).
#   y labels are plotmath expressions (genus/species italic, same rules as
#   Fig.2C), prefixed with the numbers from Step 3. LDA values are printed
#   at the bar ends with two decimals for eyeball checking (an earlier
#   manuscript draft had misread the Log10_value column as 5.26/5.17;
#   correct values are 4.18/4.38).
# ============================================================================

# ---- 5.1 display names and italic expressions ----------------------------------
rank_full <- c(p = "phylum", c = "class", o = "order", f = "family",
               g = "genus", s = "species")
disp <- ifelse(lefse$label == "Unclassified",
               paste0("Unclassified (", rank_full[lefse$rank], ")"), lefse$label)
disp <- gsub("sp__", "sp. ", disp); disp <- gsub("_", " ", disp)
make_label <- function(nm, rk) {
  if (grepl("^Unclassified", nm)) return(sprintf('"%s"', nm))
  if (grepl(" bacterium$", nm))   return(sprintf('"%s"', nm))
  if (grepl("^Candidatus ", nm))
    return(sprintf('italic("Candidatus")~"%s"', sub("^Candidatus ", "", nm)))
  if (rk == "g") return(sprintf('italic("%s")', nm))
  if (rk == "s") {
    if (grepl(" sp\\. ", nm)) {
      parts <- strsplit(nm, " ", fixed = TRUE)[[1]]
      return(sprintf('italic("%s")~"%s"', parts[1], paste(parts[-1], collapse = " ")))
    }
    return(sprintf('italic("%s")', nm))
  }
  sprintf('"%s"', nm)
}
lefse$disp <- disp
lefse$expr <- vapply(seq_len(nrow(lefse)),
                     function(i) make_label(lefse$disp[i], lefse$rank[i]), "")

# ---- 5.2 row order and number prefixes -------------------------------------------
# HD block on top (LDA decreasing), CK block below; y axis bottom-up
lefse <- lefse[order(lefse$Groups == "HD", lefse$LDA), ]
# numbers were assigned in Step 3 with the same sort key - verify, do not recompute
stopifnot(identical(lefse$no, rev(seq_len(nrow(lefse)))))
lefse$expr <- sprintf('"%d. " ~ %s', lefse$no, lefse$expr)
lefse$disp_f <- factor(lefse$disp, levels = lefse$disp)
lefse$x <- ifelse(lefse$Groups == "HD", lefse$LDA, -lefse$LDA)

# ---- 5.3 plot ------------------------------------------------------------------------
pB <- ggplot(lefse, aes(x = x, y = disp_f, fill = Groups)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.15) +
  geom_text(aes(x = x + ifelse(Groups == "HD", 0.06, -0.06),
                label = sprintf("%.2f", LDA),
                hjust = ifelse(Groups == "HD", 0, 1)),
            size = 2.1, colour = "grey30") +
  scale_y_discrete(labels = setNames(parse(text = lefse$expr), lefse$disp_f)) +
  scale_fill_manual(values = c(HD = COL_HD, CK = COL_CK),
                    labels = c(HD = "HD-enriched", CK = "CK-enriched"),
                    name = NULL) +
  scale_x_continuous(limits = c(-6.1, 6.1),
                     breaks = seq(-6, 6, 2), labels = abs(seq(-6, 6, 2))) +
  labs(x = "LDA score (log10)", y = NULL) +
  annotate("text", x = 5.9, y = nrow(lefse) + 1.2, label = "HD",
           colour = COL_HD, fontface = "bold", size = 3.2, hjust = 1) +
  annotate("text", x = -5.9, y = 0.2, label = "CK",
           colour = COL_CK, fontface = "bold", size = 3.2, hjust = 0) +
  coord_cartesian(ylim = c(0.2, nrow(lefse) + 1.6), clip = "off") +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_text(size = 7, colour = "black"),
        axis.text.x = element_text(size = 7, colour = "black"),
        axis.title.x = element_text(size = 8),
        legend.position = "none",
        axis.line.y = element_blank())

ggsave(file.path(out_dir, "Fig4B_lda.tiff"), pB,
       width = 90, height = 120, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Fig4B_lda.pdf"), pB,
       width = 90, height = 120, units = "mm", device = cairo_pdf)

cat("\n== Fig.4 panel B done ==\nOutput:", out_dir, "/Fig4B_lda.{tiff,pdf}\n")

# ============================================================================
# Step 6: assemble the final figure
#   Panel A canvas 150 x 140 mm (xlim reserves space for the left/bottom
#   label ring), panel B 90 x 120 mm. Column widths 1.3 : 1 give panel A
#   ~107 mm of the 190 mm canvas (scale 0.72: ring 76 mm, numbers ~6.5 pt
#   at publication size); panel B stretches freely without coord_fixed, so
#   its text keeps its physical size (y labels 7 pt, LDA values ~6 pt).
# ============================================================================

fig4 <- pA + pB +
  plot_layout(widths = c(1.3, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold"))

print(fig4)

ggsave(file.path(out_dir, "Fig4_final.tiff"), fig4,
       width = 190, height = 138, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Fig4_final.pdf"), fig4,
       width = 190, height = 138, units = "mm", device = cairo_pdf)

cat("\n== Fig.4 final figure done ==\nOutput:", out_dir, "/Fig4_final.{tiff,pdf}\n")
