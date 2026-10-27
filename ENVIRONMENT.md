# Software environment

Versions used to produce the reported results.

## Tool versions

| Software | Version | Purpose |
|---|---|---|
| HMMER | 3.4 | PF01036 HMM search |
| seqkit | 2.10.0 | sequence extraction, length filtering |
| MAFFT | 7.526 | multiple sequence alignment |
| trimAl | 1.5.1 | alignment trimming |
| IQ-TREE | 3.1.3 | phylogeny (ModelFinder + UFBoot 1000) |
| CD-HIT | 4.8.1 | OTU clustering at 95% amino-acid identity |
| bowtie2 | 2.5.1 | read recruitment to gene catalogue |
| samtools | 1.20 | BAM processing and counting |
| MEGAHIT | 1.2.9 | six-sample co-assembly |
| metaWRAP | 1.3.2 | binning, refinement, quantification |
| MetaBAT2 | 2.12.1 | binner 1 |
| MaxBin2 | bundled with metaWRAP 1.3.2 | binner 2 |
| CONCOCT | 1.0.0 | binner 3 |
| CheckM | 1.2.4 | completeness / contamination |
| GTDB-Tk | 2.5.2 | taxonomy, GTDB r226 |
| BLAST | 2.x (megablast) | PR nucleotide mapping to co-assembly |
| Prodigal | (metagenomic mode) | gene prediction in MAGs |
| salmon | bundled with metaWRAP | MAG abundance (TPM) |
| Python | 3.9 or later | filtering and statistics scripts |
| R | 4.5.1 | figures (ggplot2, ape, patchwork, cowplot, dplyr, tidyr, scales) |

## Databases

| Database | Release / version | Used for |
|---|---|---|
| Pfam | 38.2 | bacteriorhodopsin-like family **PF01036.25** (HMM, E-value <= 1e-10) |
| GTDB | r226 | MAG taxonomy |
| CheckM lineage data | bundled | MAG quality |
| Curated type-I rhodopsin reference set | 54 reviewed type-I rhodopsins + 295 proteorhodopsins with documented green/blue spectral types | spectral tuning-site calibration |

Vendor annotation (KEGG, CAZy, NR, eggNOG) was performed by Biomarker Technologies on releases accessed March 2025; version numbers were not supplied in the vendor report.

## Notes

- GTDB-Tk requires `GTDBTK_DATA_PATH` to point to the extracted GTDB r226 release.
- Exact Mann-Whitney tests are exact enumerations and therefore deterministic; no random seed affects the reported p-values. IQ-TREE bootstrap values reproduce within stochastic variation.
