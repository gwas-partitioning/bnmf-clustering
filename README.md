## T2D Multi-ancestry Partitioned Polygenic Scores
Cluster weights are available from the [Smith, Deutsch et al Nature Medicine 2024](https://www.nature.com/articles/s41591-024-02865-3) paper in the "Smith_Deutsch_NatureMedicine_2024" folder. Partitioned polygenic scores (pPS) can be generated using the enclosed variant cluster weights. We have included weights for ancestry-specific and multi-ancestry clusters:

* In the weights files, the "Effect_Allele" column denotes the T2D risk-increasing allele.
* When generating the pPS, all genotypes should be aligned to this allele!
* Weights have been provided in hg38, however a liftover map (hg19 to hg38) is included in each subfolder.

---

## Pipeline for GWAS clustering using Bayesian non-negative matrix factorization (bNMF)

The bNMF procedure, as applied here, detects clusters of GWAS variants for some outcome of interest based on the associations of those variants with a set of additional traits. The pipeline covers variant selection and QC, proxy variant search, z-score matrix preparation, clustering, and post-hoc summarization.

**Genome build:** hg19/GRCh37 for primary analysis; hg38 conversion is handled in the post-hoc step via liftover.

### Quick start

1. Copy `scripts/main_script_example.R` to your project folder and open it.
2. Fill in the **USER CONFIGURATION** section at the top (paths, thresholds, LDlink token).
3. Populate a copy of `example_data/clustering_data_source_example.xlsx` with paths to your own summary statistics files.
4. Run sections sequentially. Steps 2–3 (LD pruning via LDlink) are slow — run once and use `save.image()` / `load()` checkpoints to resume.

Sample input files are provided in `example_data/` to illustrate the expected formats.

### Requirements

**R packages:** `data.table`, `dplyr`, `tidyverse`, `magrittr`, `readxl`, `LDlinkR`, `furrr`, `parallelly`, `softImpute`, `strex`, `GenomicRanges`, `rtracklayer`, `Rsamtools`, `openxlsx`, `Homo.sapiens`, `vroom`

**LDlink API token:** Several steps use the [LDlinkR](https://cran.r-project.org/package=LDlinkR) package for LD-based operations. Request a free token at https://ldlink.nih.gov/?tab=apiaccess and store it as the environment variable `LDLINK_TOKEN`.

**rsID map:** The proxy search step requires a per-chromosome file mapping `VAR_ID → rsID`. Generate it from 1000 Genomes Phase 3 using `scripts/generate_varid_to_rsid_map_file.R`, which downloads per-chromosome VCFs from Ensembl and writes one map file per chromosome to a local directory. Pass the directory path as `my_rsid_map_dir` in the USER CONFIGURATION section.

### Input data format

**Variant IDs:** All summary statistics files must use `VAR_ID` in `CHR_POS_REF_ALT` format (e.g. `1_123456_A_G`), with alleles aligned consistently across files. Chromosome should be numeric (no `chr` prefix).

**Summary statistics columns** (required): `VAR_ID`, `Effect_Allele_PH`, `P_VALUE`, and either `BETA`+`SE` or `ODDS_RATIO`. An `N_PH` (sample size) column is optional — sample sizes can also be supplied via the manifest Excel file.

**GWAS manifest (Excel):** The pipeline is configured via a two-sheet Excel workbook:
* **`main_gwas` sheet** — one row per GWAS study for the primary trait being clustered. Required columns: `study`, `trait`, `population` (used together to build a per-study ID and to group studies for multi-ancestry LD pruning), `full_path` (path to summary stats file), `largest` (`Yes` for the study used to define the variant set; `No` otherwise).
* **`trait_gwas` sheet** — one row per additional trait used to characterize clusters. Required columns: `trait`, `full_path`, `sample_size`. Z-scores for each variant are computed from each trait's summary stats.

Only `full_path` is read for the file location. The `root_path` and `file` helper columns in the example are optional bookkeeping — in the example workbook `full_path` is simply their concatenation, written out as a plain string so it can be read without Excel recalculating any formulas.

See `example_data/clustering_data_source_example.xlsx` for a working example of this format.

---

### Pipeline scripts

#### `choose_variants_2025.R` — Variant selection and proxy search
* **`get_sig_snps()`** — Extract genome-wide significant variants from GWAS summary statistics
* **`get_biggest_gwas()`** — Retrieve variant info from the primary (largest) GWAS file
* **`snp_clump()`** — Position-based clumping of variants within a window
* **`ld_pruning_SNP.clip()`** — Multi-population LD pruning via LDlink SNPclip API
* **`window_to_sentinels()`** — Restrict a variant pool to within a window of LD-pruned sentinels (keeps proxy candidate fetch tractable for polygenic traits)
* **`find_variants_needing_proxies()`** — Flag strand-ambiguous, multiallelic, high-missingness, or non-TOPMed variants for proxy replacement
* **`choose_proxies()`** — Search for LD proxies via LDlinkR and select the best candidate per variant

#### `prep_bNMF_2025.R` — Summary statistics fetch and z-score matrix preparation
* **`fetch_summary_stats()`** — Retrieve z-scores and sample sizes across all traits for the variant set
* **`prep_z_matrix()`** — Final trait filters, sample size adjustment, and construction of the non-negative bNMF input matrix (N variants × 2M traits, with separate columns for positive and negative associations)

#### `run_bNMF_2025.R` — bNMF clustering
* **`BayesNMF.L2EU()`** — Core L2-Euclidean Bayesian NMF with ARD prior (auto-selects number of clusters K)
* **`run_bNMF_parallel()`** — Run multiple bNMF replicates in parallel using `furrr`
* **`summarize_bNMF()`** — Summarize replicate results and generate W/H heatmaps

#### `post_bNMF_2025.R` — Post-hoc analysis
* **`calculate_cutoff()`** — Determine optimal cluster activity weight cutoff (elbow method)
* **`v2g()`** — Variant-to-gene mapping using Open Targets V2G scores (requires local V2G database files)
* **`do_post_analysis()`** — Full post-hoc pipeline: nearest-gene annotation, Excel workbook export, hg19→hg38 liftover, PRS score file generation, and optional V2G mapping

#### `format_bNMF_results.Rmd` — Results visualization
Generates an HTML report with heatmaps, cluster summary tables, and an optimal weight threshold calculation.

---

### Outputs
Most pipeline steps print progress messages. Key output files written to the results directory:

| File | Description |
|---|---|
| `run_summary.txt` | Chosen K and negative log-likelihood for each bNMF replicate |
| `z_score_mat.rds` | Final N × M z-score matrix after all preprocessing |
| `L2EU.W.mat.K[K].txt` | Feature (trait) weight matrix for cluster count K |
| `L2EU.H.mat.K[K].txt` | Variant weight matrix for cluster count K |
| `W_plot_K[K].pdf` | Heatmap of trait contributions to clusters |
| `H_plot_K[K].pdf` | Heatmap of variant contributions to clusters |
| `quality_control_report.txt` | QC summary across all pipeline steps |
| `rsID_map.txt` | Final variant set with VAR_ID → rsID mapping |
| `alignment_GWAS_summStats.csv` | Aligned GWAS summary statistics for final variant set |

---

### Contributors
* Claire Kim (design and code)
* Kenny Westerman (design and code)
* Kirk Smith (code)
* Jaegil Kim (code)
* Marcin von Grotthuss (code)
* Miriam Udler (design and supervision)
