# Data dictionary

All tables are tab-separated UTF-8 text (no BOM). Unless stated otherwise, abundances are
relative abundances (proportions, summing to 1 per sample) or TPM.
Sample identifiers: `CK-1/2/3` = control zone, `HD-1/2/3` = kelp-farming zone.
In vendor-derived tables the read-ID equivalents are `S01_R0001–R0003` (CK) and `S01_R0004–R0006` (HD).

Tables used by more than one figure are kept in `data/`. Tables used by a single figure live
next to that figure's script, in the corresponding `FigN/` directory.

## Shared tables (`data/`)

| File | Columns | Notes |
|---|---|---|
| `group.txt` | `sample, treat1` | sample-to-group assignment used by several figure scripts |
| `lefse_rep.tsv` | `Biomarker_names, Log10_value, Groups, LDA_value, p_value` | LEfSe output; taxon string formatted `k__…p__…c__…o__…f__…g__…s__` (Fig. 2C, Fig. 4) |
| `phylum_relabund.tsv`, `class_relabund.tsv`, `order_relabund.tsv`, `family_relabund.tsv`, `genus_relabund.tsv`, `species_relabund.tsv` | first column = taxon name (unnamed header), then `CK-1 … HD-3` | relative abundance per rank; first row `Unassigned` (Fig. 2, Fig. 3, Fig. 4) |
| `pr_tpm.tsv` | `gene_id, CK-1 … HD-3, sum` | TPM per PR gene (Fig. 7, Fig. 8, Fig. 9) |

## Figure-specific tables

| File | Columns | Notes |
|---|---|---|
| `Fig3/alpha_index_stat.tsv` | `sample, observed_species, ACE, Chao1, Shannon, Simpson, goods_Coverage, pielou_evenness` | alpha-diversity indices per sample (Fig. 3A) |
| `Fig3/pcoa_bray_coordinate.tsv` | `sampleID, PC1 … PC6` | Bray–Curtis PCoA sample scores (Fig. 3B) |
| `Fig3/pcoa_bray_explained.tsv` | `axes, explained` | percentage variance explained per axis |
| `Fig5/kegg_l3_diff.tsv` | `pathway, kegg_pathway_id, KEGG_pathway level 2, KEGG_pathway level 1, mean_CK, mean_HD, ratio_HD_CK, log2FC, direction, p_exactMW, q_BH` | KEGG level-3 pathway contrast; `p_exactMW` = exact Mann–Whitney p, `q_BH` = Benjamini–Hochberg adjusted (Fig. 5A) |
| `Fig5/kelp_polysaccharide_families.tsv` | `family, activity, mean_CK, mean_HD, ratio_HD_CK, p_exactMW, direction` | kelp-polysaccharide-related CAZy families (Fig. 5B) |
| `Fig6/pr_otu_table.tsv` | `OTU_id, centroid, size, CK-1 … HD-3, consensus_spectral, consensus_motif` | 137 OTUs clustered at 95% amino-acid identity; counts per sample |
| `Fig6/kelpPR_pr.contree` | — | IQ-TREE consensus tree in Newick format, PR OTUs plus curated type-I rhodopsin references |
| `Fig7/pr_annotation_all.tsv` | `gene_id, sample, group, len, Evalue, K231, res105, spectral, motif_97_101_108, model_cov, filter` | all 139 PR genes; `K231` retinal-binding lysine, `res105` tuning-site residue, `spectral` = green/blue/unclassified, `filter` = PASS/FAIL |
| `Fig7/pr_gpm_summary.tsv` | `sample, group, PR_TPM_sum, PR_GPM, recA_meanTPM, rpoB_meanTPM, rplB_meanTPM` | PR abundance normalised to single-copy housekeeping genes |
| `Fig8/bin_abundance_table.tab` | `Genomic bins, S01_R0001 … S01_R0006` | per-MAG abundance (TPM) from metaWRAP `quant_bins` |
| `Fig8/mag_hd_vs_ck.tsv` | `MAG, GTDB_taxonomy, Completeness, Contamination, meanTPM_CK, meanTPM_HD, log2FC_HD_CK, exact_MW_p` | zone contrast for all 131 MAGs (Fig. 8A/B) |
| `Fig8/mag_pr_summary.tsv` | `MAG, n_PR, PR_proteins, mean_len_aa, GTDB_taxonomy, step7_crossval` | rhodopsin genes found in MAG proteomes; `step7_crossval` = host assignment confirmed by independent read mapping (Fig. 8C/D) |
| `Fig9/cazy_family_relabundance.tsv` | `Family, activities, CK-1 … HD-3` | CAZy family relative abundance with activity annotation |
| `Fig9/nr_genus_top.tsv` | `genus, phylum, CK-1 … HD-3, mean_CK, mean_HD, p_exactMW` | NR-annotated genus table; counts, not relative abundance |
| `Fig9/pr_cazy_coupling.tsv` | `family, spearman_prTpm_vs_family, activity` | Spearman correlation of CAZy family abundance with total PR TPM |

## Not included

Natural Earth 1:10m coastline (`coastline_ne10m.gpkg`, 7.3 MB) is third-party data and is not
redistributed here. Fetch it into `Fig1/` before running `Fig1/fig1.R`:

```bash
bash Fig1/fetch_coastline.sh
```

Two Bray/Jaccard PCoA coordinate tables produced during exploratory analysis are not used by
any figure script and are therefore not included.
