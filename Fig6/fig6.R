# ============================================================================
# Fig.6  Maximum-likelihood phylogeny of proteorhodopsins
#        (137 kelp-bay OTUs + 185 global references), circular layout
# ----------------------------------------------------------------------------
# Project : kelpPR - metagenomic study of kelp aquaculture effects on
#           bacterioplankton and proteorhodopsin-based phototrophy
# Author  : Fei Zhang
# Date    : 2026-10-27
# Version : 1.0
# Input   : Fig6/kelpPR_pr.contree (IQ-TREE consensus tree, 1000 ultrafast
#           bootstraps; tip labels follow "KelpOTU###|type|enriched-side|n"),
#           Fig6/pr_otu_table.tsv (137 OTU metadata: spectral type / counts)
# Output  : Fig6/Fig6_pr_phylogeny.{tiff,pdf} (175 x 165 mm, 300 dpi + vector)
# Requires: ape, phangorn, ggplot2
# ----------------------------------------------------------------------------
# Design choices for a 322-tip tree: kelp OTUs are point-coloured by enriched
# side (HD red / CK blue), references in grey; spectral types are shown as an
# outer ring of green/blue squares; only the four key OTUs named in the
# manuscript are annotated (OTU059 without cultured relatives; OTU055+056,
# Legionella-deep sister pair; OTU110, blue-light type ~10x enriched in HD);
# internal nodes with UFBoot >= 95 get small black dots.
# Tree-handling pitfalls handled in code:
#   (a) the IQ-TREE contree is unrooted - midpoint-root before computing
#       depths;
#   (b) ape::node.depth(method = 1) returns descendant tip counts, not
#       root-distance depths, so edge-count depths are computed by a manual
#       root-down traversal;
#   (c) internal-node angle = midpoint of the descendant tip angular span
#       (not a parent-child iterative mean, which misweights subtrees).
# ----------------------------------------------------------------------------

# ---- Step 1: package check (install missing packages from CRAN mirror) -----
need_pkgs <- c("ape", "phangorn", "ggplot2")
miss_pkgs <- need_pkgs[!is.element(need_pkgs, rownames(installed.packages()))]
if (length(miss_pkgs) > 0) {
  install.packages(miss_pkgs, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN",
                   type = "binary", quiet = TRUE)
}
suppressMessages(library(ape))
suppressMessages(library(phangorn))
library(ggplot2)

# ---- Step 2: locate working directories ------------------------------------
locate_dirs <- function(fig_id) {
  up <- if (dir.exists("data")) "." else if (dir.exists(file.path("..", "data"))) ".."
        else stop("Cannot find the project root. Set the working directory to ",
                  "the project root or to ", fig_id, "/")
  list(data_dir = file.path(up, "data"), out_dir = file.path(up, fig_id))
}
dirs <- locate_dirs("Fig6")
data_dir <- dirs$data_dir
out_dir  <- dirs$out_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

COL_HD <- "#C44E52"; COL_CK <- "#4C72B0"
COL_GREEN <- "#3C9A5F"; COL_BLUE <- "#3C5488"

# ---- Step 3: read tree + OTU metadata ----------------------------------------
tr <- read.tree(file.path(out_dir, "kelpPR_pr.contree"))
tr <- midpoint(tr)                       # pitfall (a): unrooted tree, midpoint-root first
stopifnot(is.rooted(tr))
tr <- ladderize(tr)                      # ladderised for a stable visual order
otu <- read.delim(file.path(out_dir, "pr_otu_table.tsv"),
                  header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(length(tr$tip.label) == 322, nrow(otu) == 137)
cat("spectral types (137 OTUs):", paste(names(table(otu$consensus_spectral)),
                                   table(otu$consensus_spectral), collapse = ", "), "\n")

# ---- Step 4: tip metadata parsing ----------------------------------------------
tips <- tr$tip.label
is_kelp <- grepl("^KelpOTU", tips)
kelp_grp <- vapply(strsplit(tips[is_kelp], "|", fixed = TRUE),
                   function(x) x[3], "")                   # HD / CK
tip_meta <- data.frame(
  tip = tips, kelp = is_kelp,
  group = ifelse(is_kelp, kelp_grp[match(seq_along(tips)[is_kelp],
                                         which(is_kelp))], NA),
  stringsAsFactors = FALSE)
# spectral type: strip the "|..." suffix first, then the "Kelp" prefix
otu_idx <- match(sub("Kelp", "", sub("\\|.*$", "", tip_meta$tip)), otu$OTU_id)
tip_meta$spectral <- ifelse(tip_meta$kelp, otu$consensus_spectral[otu_idx], NA)
stopifnot(!any(is.na(tip_meta$group[tip_meta$kelp])),
          !any(is.na(tip_meta$spectral[tip_meta$kelp])))

# ---- Step 5: circular layout ------------------------------------------------------
n_tip <- length(tips)
ang_tip <- seq(0, 2 * pi * (1 - 1 / n_tip), length.out = n_tip)
# pitfall (c): internal-node angle = midpoint of the descendant tip span
# (subtrees of a ladderised tree are contiguous in tip order); memoised to
# avoid exponential recursion
kids <- split(tr$edge[, 2], tr$edge[, 1])
span_env <- new.env()
span_of <- function(nd) {
  key <- as.character(nd)
  if (exists(key, envir = span_env)) return(get(key, envir = span_env))
  sp <- if (nd <= n_tip) rep(ang_tip[nd], 2)
        else {
          s <- vapply(kids[[key]], span_of, numeric(2))
          c(min(s[1, ]), max(s[2, ]))
        }
  assign(key, sp, envir = span_env)
  sp
}
all_nodes <- c(seq_len(n_tip), n_tip + seq_len(tr$Nnode))
spans <- vapply(all_nodes, span_of, numeric(2))
ang <- setNames((spans[1, ] + spans[2, ]) / 2, all_nodes)
stopifnot(!any(is.na(ang)))

# radii: pitfall (b) - edge-count depths computed by a manual root-down
# traversal (cladogram; the caption requires no scale bar, "deep branches"
# are conveyed by topology)
n_all <- n_tip + tr$Nnode
dep <- rep(NA_integer_, n_all)
root_nd <- setdiff(n_tip + seq_len(tr$Nnode), tr$edge[, 2])
stopifnot(length(root_nd) == 1)
dep[root_nd] <- 0L
for (i in seq_len(nrow(tr$edge)))
  dep[tr$edge[i, 2]] <- dep[tr$edge[i, 1]] + 1L
stopifnot(!any(is.na(dep)))
r <- dep / max(dep[seq_len(n_tip)])
xy <- data.frame(node = seq_along(ang), theta = as.numeric(ang), r = r,
                 x = r * cos(as.numeric(ang)), y = r * sin(as.numeric(ang)),
                 stringsAsFactors = FALSE)

# ---- Step 6: edge geometry (parent arc + radial segment) --------------------------
edge_parts <- lapply(seq_len(nrow(tr$edge)), function(i) {
  p <- tr$edge[i, 1]; c <- tr$edge[i, 2]
  arc <- seq(xy$theta[p], xy$theta[c], length.out = 20)
  rbind(data.frame(x = xy$r[p] * cos(arc), y = xy$r[p] * sin(arc),
                   grp = i, kid = c),
        data.frame(x = c(xy$r[p] * cos(xy$theta[c]), xy$r[c] * cos(xy$theta[c])),
                   y = c(xy$r[p] * sin(xy$theta[c]), xy$r[c] * sin(xy$theta[c])),
                   grp = i, kid = c))
})
edges <- do.call(rbind, edge_parts)
kid_grp <- tip_meta$group[match(edges$kid, seq_len(n_tip))]
edges$kind <- ifelse(!is.na(kid_grp), kid_grp, "bg")

# ---- Step 7: node / tip / annotation layer data --------------------------------------
tip_df <- data.frame(tip = tips, theta = ang[seq_len(n_tip)],
                     r = r[seq_len(n_tip)],
                     x = xy$x[seq_len(n_tip)], y = xy$y[seq_len(n_tip)],
                     kelp = tip_meta$kelp, group = tip_meta$group,
                     spectral = tip_meta$spectral, stringsAsFactors = FALSE)

# UFBoot support (contree internal node labels)
boot <- suppressWarnings(as.numeric(tr$node.label))
boot_df <- data.frame(x = xy$x[n_tip + seq_len(tr$Nnode)],
                      y = xy$y[n_tip + seq_len(tr$Nnode)], boot = boot)
boot_df <- boot_df[is.finite(boot_df$boot) & boot_df$boot >= 95, ]
cat("internal nodes with UFBoot >= 95:", nrow(boot_df), "/", tr$Nnode, "\n")

# key-OTU annotations (spread with minimum angular distance, same approach
# as the Fig.4A label ring). On-figure only bold numbers 1-4; the full
# descriptions go into the manuscript caption
key_otus <- c("OTU059", "OTU055", "OTU056", "OTU110")
key_no <- setNames(as.character(seq_along(key_otus)), key_otus)
key_note <- c(OTU059 = "1", OTU055 = "2", OTU056 = "3", OTU110 = "4")
tip_df$otu_id <- sub("Kelp", "", sub("\\|.*$", "", tip_df$tip))
key_df <- tip_df[is.element(tip_df$otu_id, key_otus), ]
stopifnot(nrow(key_df) == length(key_otus))
key_df$lab <- key_note[key_df$otu_id]
key_df <- key_df[order(key_df$theta), ]
slot <- key_df$theta
if (nrow(key_df) > 1) {
  for (i in 2:nrow(key_df))
    if (slot[i] - slot[i - 1] < 0.16) slot[i] <- slot[i - 1] + 0.16
}
slot <- slot - mean(slot - key_df$theta)
key_df$slot <- atan2(sin(slot), cos(slot))
# numbers sit close to their tips (r + 0.43, past the spectral squares at
# r + 0.045), hjust centred; leaders end at r + 0.40 (~6.5 mm)
key_df$lab_r <- key_df$r + 0.43
key_df$lab_x <- key_df$lab_r * cos(key_df$slot)
key_df$lab_y <- key_df$lab_r * sin(key_df$slot)
key_df$hjust <- 0.5

# ---- Step 8: plot and export ---------------------------------------------------------
p <- ggplot() +
  geom_path(data = edges[edges$kind == "bg", ], aes(x, y, group = grp),
            colour = "grey72", linewidth = 0.3) +
  geom_path(data = edges[edges$kind == "HD", ], aes(x, y, group = grp),
            colour = COL_HD, linewidth = 0.5, alpha = 0.8) +
  geom_path(data = edges[edges$kind == "CK", ], aes(x, y, group = grp),
            colour = COL_CK, linewidth = 0.5, alpha = 0.8) +
  geom_point(data = tip_df[!tip_df$kelp, ], aes(x, y), size = 0.35,
             colour = "grey55") +
  geom_point(data = tip_df[tip_df$kelp, ], aes(x, y, colour = group),
             size = 1.15) +
  # spectral squares sit just outside each tip's own radius (+0.045); under a
  # cladogram tips have different depths, a fixed radius would float the
  # squares off the crown
  geom_point(data = tip_df[tip_df$kelp, ],
             aes(x = (r + 0.045) * cos(theta), y = (r + 0.045) * sin(theta),
                 colour = spectral), size = 0.9, shape = 15) +
  geom_point(data = boot_df, aes(x, y), size = 0.35, colour = "black") +
  geom_segment(data = key_df,
               aes(x = x, y = y,
                   xend = (r + 0.40) * cos(slot), yend = (r + 0.40) * sin(slot)),
               colour = "grey45", linewidth = 0.25) +
  geom_text(data = key_df, aes(x = lab_x, y = lab_y, label = lab,
                               hjust = hjust),
            size = 3.2, colour = "black", fontface = "bold",
            show.legend = FALSE) +
  scale_colour_manual(values = c(HD = COL_HD, CK = COL_CK,
                                 green = COL_GREEN, blue = COL_BLUE),
                      breaks = c("HD", "CK", "green", "blue"),
                      labels = c("kelp OTU, HD-enriched",
                                 "kelp OTU, CK-enriched",
                                 "green-light PR", "blue-light PR"),
                      name = NULL) +
  guides(colour = guide_legend(override.aes = list(
    size = 1.6, shape = c(16, 16, 15, 15)))) +
  coord_fixed(xlim = c(-1.10, 1.10), ylim = c(-1.10, 1.10), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(3, 3, 3, 3, unit = "mm"),
        legend.position = c(0.005, 0.01), legend.justification = c(0, 0),
        legend.background = element_blank(), legend.key = element_blank(),
        legend.text = element_text(size = 6.3),
        legend.spacing.y = unit(0.05, "cm"))

ggsave(file.path(out_dir, "Fig6_pr_phylogeny.tiff"), p,
       width = 175, height = 165, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Fig6_pr_phylogeny.pdf"), p,
       width = 175, height = 165, units = "mm", device = cairo_pdf)

cat("\n== Fig.6 done ==\n")
cat("kelp OTU tips:", sum(tip_df$kelp), "(HD", sum(tip_df$group == "HD", na.rm = TRUE),
    "/ CK", sum(tip_df$group == "CK", na.rm = TRUE), "); reference tips:",
    sum(!tip_df$kelp), "\n")
cat("output:", out_dir, "/Fig6_pr_phylogeny.{tiff,pdf}\n")
cat("\n== number map (on-figure number -> OTU) ==\n")
for (o in key_otus) cat("  ", key_no[o], "=", o, "\n")
cat("\ncaption sentence for the manuscript (Fig. 6):\n",
    "Key OTUs are numbered 1-4: 1, OTU059 (no cultured relative); ",
    "2, OTU055 and 3, OTU056 (Legionella-deep sister pair); ",
    "4, OTU110 (blue-light, ~10x enriched in HD).\n")
