# Code and data for a metagenomic study of kelp aquaculture effects on bacterioplankton and proteorhodopsin-based phototrophy

Scripts and input tables for the analyses reported in the associated study. Six surface-water metagenomes were collected in Dongshan Bay, Fujian, China (14 January 2025) from a kelp-farming zone (HD-1 to HD-3) and an adjacent non-farmed zone (CK-1 to CK-3).

## Contents

| Directory | Contents |
|---|---|
| `01_pr_analysis/` | PR identification, key-residue filtering, functional-site annotation, OTU clustering, phylogeny, read-recruitment abundance (TPM, GPM) |
| `02_mag/` | six-sample co-assembly, binning and hybrid bin refinement (131 MAGs) |
| `03_pr_mag/` | assignment of PR genes to MAGs by sequence-level mapping |
| `04_mag_abundance/` | MAG quantification and zone contrasts (TPM, log2 fold change) |
| `05_mag_pr/` | independent rhodopsin screen inside MAGs (Prodigal + PF01036) |
| `06_functional/` | KEGG and CAZy profiling, PR-CAZy coupling |
| `Fig1/` to `Fig9/` | one R script per figure plus the input tables used only by that figure |
| `data/` | input tables shared by several figures (UTF-8 tab-separated text) and the data dictionary |

## Requirements

- R 4.5.1: `ggplot2`, `dplyr`, `tidyr`, `ape`, `pheatmap`, `patchwork`, `cowplot`, `scales`, `ggtext`
- Python 3.9+: `pandas`, `Biopython`
- HMMER 3.4, seqkit 2.10.0, CD-HIT 4.8.1, MAFFT 7.526, trimAl 1.5.1, IQ-TREE 3.1.3, bowtie2 2.5.1, samtools 1.20, MEGAHIT 1.2.9, metaWRAP 1.3.2, CheckM 1.2.4, GTDB-Tk 2.5.2 (GTDB r226), BLAST (megablast), Prodigal, salmon

`ENVIRONMENT.md` lists the versions used. Conda environments are configured by the user; activate the required environment before each stage.

## Running

Stages run in numeric order from the local clone. Each script writes to a project directory taken from the `PROJECT_DIR` environment variable (default `~/projects/kelpPR`) and skips outputs that already exist.

```bash
bash   01_pr_analysis/step1_pr_search.sh
bash   01_pr_analysis/step2_pipeline.sh
python 01_pr_analysis/step3_otu_table.py
bash   01_pr_analysis/step4_tree.sh
bash   01_pr_analysis/step5_pr_abundance.sh
python 01_pr_analysis/step5_tpm.py
bash   02_mag/step6_binning.sh
bash   03_pr_mag/step7_pr_mag_map.sh
python 03_pr_mag/step7_join.py
bash   04_mag_abundance/step8_mag_abundance.sh
python 04_mag_abundance/step8_compare.py
bash   05_mag_pr/step9_mag_pr_screen.sh
python 05_mag_pr/step9_summarize.py
python 06_functional/step10_functional_analysis.py
```

Each figure is produced by a single script, `FigN/figN.R`, which builds every panel of that figure in numbered steps:

```bash
Rscript Fig1/fig1.R
Rscript Fig3/fig3.R
Rscript Fig6/fig6.R
```

Scripts may be run from the repository root or from their own `FigN` directory; they resolve the shared `data/` directory and their figure-specific inputs automatically, and write TIFF (300 dpi) and PDF versions of the figure to the `FigN` directory.

## Data sources

Raw metagenomic reads are available from NCBI under BioProject PRJNA1532396 (BioSamples SAMN63273868 to SAMN63273873; runs SRR40761233 to SRR40761238). All tables in `data/` and in the `FigN` directories derive from those reads; their column definitions are given in `data/README_data_dictionary.md`. The coastline basemap used in the first figure is not bundled and is retrieved with `bash Fig1/fetch_coastline.sh`.

File integrity can be checked against `checksums_sha256.txt`.

## Licence

Code: MIT, see `LICENSE`. Data and documentation: CC BY 4.0.

## Citation

Cite this archive through `CITATION.cff`, which GitHub and Zenodo read automatically.
