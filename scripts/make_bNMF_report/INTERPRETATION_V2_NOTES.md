# bNMF interpretation report

An HTML report built from the files a run writes to its results directory. It
makes no web requests and needs no annotation databases.

## Files

- `format_bNMF_results_2025_INTERPRETATION_V2.Rmd`: report template.
- `make_cluster_figures_INTERPRETATION_V2.R`: plots and interpretation helpers.
- `render_bNMF_report_INTERPRETATION_V2.R`: command-line renderer.

## Render after a pipeline run

The results directory must contain the output of `do_post_analysis()`
(`cluster_post_hoc.rds`) plus `zmat_preImputed.rds`, `zmat_imputed.rds`,
`rsID_map.txt`, `run_summary.txt` and `df_traits.csv`.

```sh
Rscript scripts/make_bNMF_report/render_bNMF_report_INTERPRETATION_V2.R /path/to/CURRENT_results
```

To compare signed trait profiles with a prior iteration:

```sh
Rscript scripts/make_bNMF_report/render_bNMF_report_INTERPRETATION_V2.R \
  /path/to/CURRENT_results /path/to/PREVIOUS_results
```

The renderer determines K from `cluster_post_hoc.rds`. The Rmd rejects a
conflicting `params$k`, so a report cannot show one K while using another K's
W/H matrices.

## What it shows

Run summary and K frequency, trait and variant processing (with the run's
`QC_report.html` when present), convergence, the W-only weight cutoff and
variant membership, trait core membership, and one section per cluster: defining
signed traits, defining loci with soft-membership specificity, trait-domain
composition, a circular weight plot and a z-score heatmap. Gene labels are the
nearest gene to each variant.

## Interpretation semantics

- Trait direction is derived from paired `_pos` and `_neg` H components and refers to effects aligned to the primary-GWAS risk allele.
- Locus specificity is the weight in one cluster divided by total weight across clusters.
- "Effective" component counts use inverse concentration (`1 / sum(p^2)`).
- Provisional signatures come from transparent pattern-based trait domains in `classify_trait_domain()`; review domain assignments when adding new trait naming systems. A `Group` column in the trait manifest (written to `trait_groups.csv`) overrides them.
- Dominance flags trigger when one locus or trait contributes at least 40% of its above-cutoff weight. "Mostly shared" triggers when at least 50% of loci also exceed the cutoff elsewhere.

## Current limitation

The pipeline saves one summarized W/H solution per K, not aligned W/H matrices from every repetition. The report can therefore show K frequency and evidence variation but cannot estimate per-cluster reproducibility across repetitions.
