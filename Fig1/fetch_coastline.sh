#!/usr/bin/env bash
# Fetch the Natural Earth 1:10m coastline basemap used by Fig1/fig1.R.
# The geopackage is third-party data (public domain) and is therefore not redistributed
# inside this repository.
#
# Usage:  bash Fig1/fetch_coastline.sh [output_dir]
set -euo pipefail

OUTDIR="${1:-$(dirname "$0")}"
URL="https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_coastline.zip"
TARGET="${OUTDIR}/coastline_ne10m.gpkg"

if [ -f "${TARGET}" ]; then
  echo "Already present: ${TARGET}"
  exit 0
fi

TMP="$(mktemp -d)"
trap 'rm -rf "${TMP}"' EXIT

echo "Downloading ${URL}"
curl -L --fail -o "${TMP}/ne.zip" "${URL}"
unzip -q -o "${TMP}/ne.zip" -d "${TMP}/ne"

# Convert the bundled shapefile to GeoPackage (requires GDAL >= 3, e.g. via sf::st_write in R).
SHP="$(find "${TMP}/ne" -name 'ne_10m_coastline.shp' | head -n 1)"
if [ -z "${SHP}" ]; then
  echo "Shapefile not found in archive" >&2
  exit 1
fi

if command -v ogr2ogr >/dev/null 2>&1; then
  ogr2ogr -f GPKG "${TARGET}" "${SHP}"
else
  Rscript -e "sf::st_write(sf::st_read('${SHP}', quiet = TRUE), '${TARGET}', delete_dsn = TRUE, quiet = TRUE)"
fi

echo "Written: ${TARGET}"
echo "Source: Natural Earth (https://www.naturalearthdata.com/), public domain."
