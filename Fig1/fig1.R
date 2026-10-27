# ============================================================================
# Fig.1  Sampling-station map of the study area (standalone figure)
# ----------------------------------------------------------------------------
# Project : kelpPR - metagenomic study of kelp aquaculture effects on
#           bacterioplankton and proteorhodopsin-based phototrophy
# Author  : Fei Zhang
# Date    : 2026-10-27
# Version : 1.0
# Input   : Fig1/coastline_ne10m.gpkg (Natural Earth 10m coastline cache,
#           public domain)
# Output  : Fig1/Fig1_station_map.tiff (300 dpi) and Fig1/Fig1_station_map.pdf
#           (vector)
# Requires: sf, ggspatial, ggplot2, patchwork
# Note    : Natural Earth natural coastline only, no borders or administrative
#           units; inset extent <= 119.3 degE so the frame shows only the
#           coastal seas of the study region
# ----------------------------------------------------------------------------

# ---- Step 1: package check (install missing packages from CRAN mirror) -----
need_pkgs <- c("sf", "ggspatial", "ggplot2", "patchwork")
miss_pkgs <- need_pkgs[!is.element(need_pkgs, rownames(installed.packages()))]
if (length(miss_pkgs) > 0) {
  install.packages(miss_pkgs, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN",
                   type = "binary", quiet = TRUE)
}
library(sf)
library(ggplot2)
library(ggspatial)
library(patchwork)

# ---- Step 2: locate working directories ------------------------------------
# Run from the project root (or from this script's folder: the script walks
# one level up to find the root). Figure-specific inputs live in Fig1/.
locate_dirs <- function(fig_id) {
  up <- if (dir.exists("data")) "." else if (dir.exists(file.path("..", "data"))) ".."
        else stop("Cannot find the project root. Set the working directory to ",
                  "the project root or to ", fig_id, "/")
  list(data_dir = file.path(up, "data"), out_dir = file.path(up, fig_id))
}
dirs <- locate_dirs("Fig1")
data_dir <- dirs$data_dir
out_dir  <- dirs$out_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Step 3: load coastline (Natural Earth 10m, local cache) ---------------
gpkg <- file.path(out_dir, "coastline_ne10m.gpkg")
stopifnot(file.exists(gpkg))
land <- st_read(gpkg, quiet = TRUE)
stopifnot(st_crs(land)$epsg == 4326)   # must be WGS84 to match station coords

# ---- Step 4: station table (DMS coords of the paper -> decimal degrees) ----
# lab_x/lab_y are label anchors verified at pixel level: labels sit close to
# the stations without overlapping the points
stn <- data.frame(
  site = c("HD", "CK"),
  desc = c("kelp farm", "control"),
  lon  = c(117 + 32/60, 117 + 19/60),   # HD 117d32'E, CK 117d19'E
  lat  = c( 23 + 48/60,  23 + 39/60),   # HD  23d48'N, CK  23d39'N
  lab  = c("HD\n(kelp farm)", "CK\n(control)"),
  lab_x = c(117.533, 117.317),
  lab_y = c(23.772,  23.627)
)

# ---- Step 5: crop coastline for main map and inset --------------------------
# Main-map xmax = 117.74 keeps the rightmost tick label (117.7 degE) clear of
# the panel edge; combined with the wider right plot margin in Step 6 this
# prevents the last character of the longitude label from being clipped
ext_main  <- c(xmin = 117.10, xmax = 117.74, ymin = 23.35, ymax = 23.95)
ext_inset <- c(xmin = 116.30, xmax = 119.30, ymin = 22.70, ymax = 25.30)
# st_make_valid first fixes self-intersecting polygons (common in raw
# Natural Earth data), then crop
land_main  <- st_crop(st_make_valid(land), st_bbox(ext_main,  crs = st_crs(4326)))
land_inset <- st_crop(st_make_valid(land), st_bbox(ext_inset, crs = st_crs(4326)))

# ---- Step 6: main map -------------------------------------------------------
p_main <- ggplot() +
  # (1) ocean base colour: fill the whole panel with sea colour, land on top
  geom_rect(aes(xmin = ext_main["xmin"], xmax = ext_main["xmax"],
                ymin = ext_main["ymin"], ymax = ext_main["ymax"]),
            fill = "#D6E6F2", colour = NA) +
  # (2) land polygons
  geom_sf(data = land_main, fill = "#F2EFE9", colour = "grey35", linewidth = 0.35) +
  # (3) water / place names (cartographic convention: water italic, land roman)
  annotate("text", x = 117.5465, y = 23.835, label = "Dongshan\nBay",
           fontface = "italic", colour = "#33668C", size = 2.8, lineheight = 0.85) +
  annotate("text", x = 117.425, y = 23.700, label = "Dongshan Is.",
           fontface = "italic", colour = "grey35", size = 2.8) +
  annotate("text", x = 117.230, y = 23.578, label = "Zhao'an Bay",
           fontface = "italic", colour = "#33668C", size = 3.2) +
  # (4) stations: points + two-line labels; HD red / CK blue (figure-wide palette)
  geom_point(data = stn, aes(x = lon, y = lat, colour = site), size = 2.8) +
  geom_text(data = stn, aes(x = lab_x, y = lab_y, label = lab, colour = site),
            fontface = "bold", size = 2.2, lineheight = 0.85, show.legend = FALSE) +
  scale_colour_manual(values = c(CK = "#4C72B0", HD = "#C44E52"), guide = "none") +
  # (5) north arrow (top left) + scale bar (bottom left)
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_minimal(text_size = 9),
                         pad_x = unit(0.14, "in"), pad_y = unit(0.16, "in")) +
  annotation_scale(location = "bl", width_hint = 0.30, text_cex = 0.75,
                   pad_x = unit(0.14, "in"), pad_y = unit(0.14, "in")) +
  # lon/lat ticks every 0.1 degrees
  scale_x_continuous(breaks = 117.1 + 0.1 * (0:6)) +
  scale_y_continuous(breaks = 23.4 + 0.1 * (0:5)) +
  coord_sf(xlim = c(ext_main["xmin"], ext_main["xmax"]),
           ylim = c(ext_main["ymin"], ext_main["ymax"]),
           expand = FALSE, datum = st_crs(4326)) +
  labs(x = "Longitude (\u00b0E)", y = "Latitude (\u00b0N)") +
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_line(colour = "grey75", linewidth = 0.25,
                                        linetype = "dashed"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#D6E6F2"),
        axis.text = element_text(size = 8, colour = "grey20"),
        axis.title = element_text(size = 9),
        # right margin widened (16 pt): the rightmost longitude label
        # "117.7 degE" extends ~0.19 in right of its tick, the default
        # 5.5 pt (0.076 in) margin would clip its last character
        plot.margin = margin(5.5, 16, 5.5, 5.5, "pt"))

# ---- Step 7: regional inset (Fujian-Guangdong coast, red frame = extent) ----
p_inset <- ggplot() +
  geom_rect(aes(xmin = ext_inset["xmin"], xmax = ext_inset["xmax"],
                ymin = ext_inset["ymin"], ymax = ext_inset["ymax"]),
            fill = "#D6E6F2", colour = NA) +
  geom_sf(data = land_inset, fill = "#F2EFE9", colour = "grey55", linewidth = 0.25) +
  annotate("rect", xmin = ext_main["xmin"], xmax = ext_main["xmax"],
           ymin = ext_main["ymin"], ymax = ext_main["ymax"],
           fill = NA, colour = "#C44E52", linewidth = 0.7) +
  annotate("text", x = 117.95, y = 24.90, label = "Fujian",
           size = 2.6, colour = "grey25") +
  annotate("text", x = 118.90, y = 23.75, label = "Taiwan\nStrait",
           size = 2.3, colour = "#33668C", fontface = "italic", lineheight = 0.85) +
  annotate("text", x = 118.45, y = 22.80, label = "South China Sea",
           size = 2.3, colour = "#33668C", fontface = "italic") +
  coord_sf(xlim = c(ext_inset["xmin"], ext_inset["xmax"]),
           ylim = c(ext_inset["ymin"], ext_inset["ymax"]), expand = FALSE) +
  theme_void() +
  theme(panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.5))

# ---- Step 8: embed inset into main map and export both formats --------------
# inset_element(left/bottom/right/top): relative panel coordinates (0-1)
p_full <- p_main +
  inset_element(p_inset, left = 0.55, bottom = 0.01, right = 0.995, top = 0.355,
                align_to = "panel", on_top = TRUE)
print(p_full)

ggsave(file.path(out_dir, "Fig1_station_map.tiff"), p_full,
       width = 4.25, height = 4.9, units = "in", dpi = 300)
ggsave(file.path(out_dir, "Fig1_station_map.pdf"), p_full,
       width = 4.25, height = 4.9, units = "in")

cat("\n== Fig.1 done ==\nOutput: ", out_dir, "/Fig1_station_map.{tiff,pdf}\n")
