# Z-score alignment audit against the GWAS Catalog — design

Date: 2026-11-13

> **Revision (2026-11-13, during implementation):** the REST API v2 proved too
> slow for whole-trait queries (about 45 s per page of 500; BMI alone about
> 40 minutes). The audit instead reads the Catalog's bulk file
> `gwas-catalog-associations_ontology-annotated-full.zip` (about 70 MB, about
> 6 s with `data.table::fread`), passed with `--catalog`. Consequences:
> `httr2` and the cache directory are gone; `LOOKUP_ERROR` is replaced by an
> "EFO ID(s) not in the Catalog file" note and report flag; rows have no
> association ID, so comparisons carry `catalog_row`; direction is "increase"
> or "decrease" in the CI text, else an unsigned value above 1 read as an
> odds ratio, else none; merged rsIDs also match through `SNP_ID_CURRENT`;
> term IDs include OBA and MONDO (pattern `^[A-Za-z]+_[A-Za-z0-9]+$`); draft
> mode suggests terms from the file's own annotations (with abbreviation
> expansion and word-set matching) and needs no network. Matrix rows named
> `chr1:123` match alignment rows named `1:123`. The sections below describe
> the original API design except where this note overrides them.

## Purpose

Check that every trait column in the pre-imputed bNMF z-score matrix is
oriented to the T2D risk allele. The most likely source of an allele bug is the
user's own preprocessing of trait summary statistics (for example, choosing the
wrong effect-allele column), which flips or scrambles a whole trait column.

The check must be independent of the user's summary-statistics files: any
comparison against those files (or the same GWAS via OpenGWAS) can reproduce
the same preprocessing error. The reference is therefore the NHGRI-EBI GWAS
Catalog's curated associations.

This replaces `T2DGGI_2025/R/audit_preimputed_zscores.R`, which compared
against the user's own files, required exact study accessions, compared Z
magnitudes reconstructed from rounded P-values, and judged a fixed sample of
~18 cells instead of each trait.

## Constraints and findings

- `gwasrapidd` (0.99.18, latest on CRAN) is unusable: the Catalog retired REST
  API v1 in August 2026 and every v1 endpoint returns HTTP 410. The maintainer
  plans a successor (gwasrapidd2) with no date.
- REST API v2 (`https://www.ebi.ac.uk/gwas/rest/api/v2/`) works and supports
  `associations?efo_id=...` (paged), `associations?efo_id=...&rs_id=...`, and
  `efo-traits?trait=...` (name search). Each association returns
  `snp_allele` (rs_id, effect_allele), `p_value`, `beta` text
  ("0.018 unit decrease"), `or_per_copy_num`, `accession_id`, `efo_traits`,
  `reported_trait`, `multi_snp_haplotype`, `snp_interaction`.
- Matrix rows are hg19 `CHR:POS`; the Catalog is GRCh38. Matching is by rsID
  from `alignment_GWAS_summStats.csv`, so the build difference does not matter.

## Approach

Sign concordance per trait. For each trait, take genome-wide-significant
Catalog associations for that trait's EFO term(s) from any study, match them to
matrix SNPs by rsID, derive the expected sign relative to the T2D risk allele,
and compare with the sign of the matrix z-score. A correctly aligned trait
agrees near 100%; a flipped column near 0%; a scrambled one near 50%.

Data access: page through `associations?efo_id=` per EFO term and cache the
result on disk. Rejected alternatives: per-rsID queries (~986 × 94 requests)
and the Catalog's bulk download (large, loosely structured columns).

## Files

- `scripts/audit_zscore_alignment.R` — sourceable functions plus a command-line
  entry point (`if (sys.nframe() == 0L)`), following the other scripts' style.
  Dependencies: `httr2` and base R. HTML is built with `qc_html_escape()` and
  `qc_html_table()` sourced from `scripts/qc_report.R`. The plot is
  hand-written inline SVG (no `svglite`).
- `tests/test-zscore-audit.R` — runs on saved Catalog JSON fixtures with no
  network access; picked up by `tests/run_tests.R`.
- README: a short section describing both modes.

## Inputs

- `zmat_preImputed.rds` — numeric matrix, rownames `CHR:POS`, colnames traits.
- `alignment_GWAS_summStats.csv` — `ChrPos`, `rsID`, `REF`, `ALT`,
  `Risk_Allele`, one row per matrix row.
- Audit TSV (separate from the pipeline's clustering xlsx), columns:
  - `trait` — matrix column name (required).
  - `efo_id` — one or more EFO IDs separated by `;`. Blank means skip the trait.
  - `accession` — optional GCST accession of the study the clustering used.
  - `reported_trait` — optional regular expression; only associations whose
    reported trait matches are kept (e.g. to separate BMI-adjusted traits).
  - `suggestions` — written by draft mode; ignored by audit mode.

## Modes

Draft: `Rscript scripts/audit_zscore_alignment.R --draft zmat_preImputed.rds audit.tsv`

Searches `efo-traits?trait=` with each matrix column name and writes one row
per trait: `efo_id` set to the top hit, `accession` and `reported_trait`
blank, `suggestions` holding the top three hits as `EFO_id (label)`. Traits
with no hit get a blank `efo_id`. Refuses to overwrite an existing file. The
user reviews and edits the TSV before auditing.

Audit: `Rscript scripts/audit_zscore_alignment.R zmat_preImputed.rds alignment_GWAS_summStats.csv audit.tsv <output_prefix> [--min-abs-z 2]`

## Comparison rules

1. Keep Catalog associations that are single-SNP (not haplotype or
   interaction), have exactly one rsID/effect-allele pair, have
   P ≤ 5×10⁻⁸, and have a determinable direction: `beta` text containing
   "increase"/"decrease", or `or_per_copy_num` ≠ 1. Apply the `reported_trait`
   filter when given.
2. Match to matrix rows by rsID. Drop A/T and C/G SNPs (by matrix REF/ALT),
   effect alleles that are `?` or not a single base, and effect alleles that
   match neither REF/ALT nor their strand complement.
3. Expected sign = Catalog direction × (+1 if the aligned effect allele is the
   T2D risk allele, −1 otherwise).
4. One call per SNP–trait pair. The pool is the associations from the
   trait's `accession` when any exist, otherwise all retained associations;
   the call uses the smallest-P association in the pool. If associations in
   the pool disagree in expected sign, mark the pair `catalog_conflict` and
   exclude it from concordance. (Other studies do not override the
   clustering study.) Note that for a non-palindromic SNP every single base
   matches REF/ALT on one strand or the other, so the strand rule in step 2
   only ever reads an allele on the opposite strand; it never drops one.
5. A pair is informative when the matrix z is finite and |z| ≥ `min_abs_z`
   (default 2). Non-informative pairs are reported but not counted.

## Verdicts (per trait)

Over informative pairs, with n = count and c = fraction concordant:

- `INSUFFICIENT` — n < 5.
- `PASS` — c ≥ 0.8.
- `FLIPPED` — c ≤ 0.2.
- `REVIEW` — otherwise.

Each trait also gets a two-sided binomial P against 0.5, concordance among
same-study pairs only (when `accession` is given), and a `study_differs` flag
when no informative pair came from the clustering study. Traits with a blank
`efo_id` are listed as `NOT_CHECKED`; Catalog lookup failures as
`LOOKUP_ERROR` with the message.

## Outputs

- `<prefix>.csv` — one row per SNP–trait comparison: trait, rsID, ChrPos,
  REF, ALT, Risk_Allele, association_id, accession, same_study, reported_trait,
  p_value, catalog effect allele, strand_complemented, catalog direction,
  expected sign, matrix z, informative, concordant, status.
- `<prefix>_traits.csv` — one row per trait: the manifest columns, verdict,
  counts, concordance, binomial P, same-study figures, study_differs, note.
- `<prefix>.html` — self-contained report: run summary; per-trait table sorted
  FLIPPED, REVIEW, LOOKUP_ERROR, INSUFFICIENT, PASS, NOT_CHECKED; inline SVG
  dot plot of concordance by trait (point size by n, shaded ≤ 0.2 and
  0.2–0.8 bands); a collapsible `<details>` per trait with its SNP rows;
  a methods and limitations note.
- `<prefix>_catalog_cache/` — one RDS per EFO ID of the raw, parsed
  associations. Reused on rerun; delete to refresh.

## Error handling

- Validate inputs up front: matrix shape and unique names, alignment columns
  and one row per matrix row, audit TSV columns, traits present in the matrix,
  EFO IDs matching `^[A-Za-z]+_[0-9]+$`.
- HTTP requests use a timeout and retries (`httr2::req_retry`). A failure for
  one EFO term marks the affected traits `LOOKUP_ERROR` and the run continues.
- The script exits non-zero only on invalid inputs, not on FLIPPED/REVIEW
  verdicts; verdict counts are printed at the end.

## Testing

Fixture JSON pages imitate v2 responses. Tests cover: direction parsing
(beta text, OR above/below 1, missing); filtering (P threshold, haplotypes,
`?` alleles, reported_trait regex); allele alignment (direct, complement,
palindromic dropped, mismatched dropped); expected-sign flip when the effect
allele is the non-risk allele; same-study preference and `catalog_conflict`;
verdict thresholds including a synthetic fully flipped trait and a scrambled
trait; paging and cache reuse (no second request); draft-mode output and
refusal to overwrite; HTML written and containing each trait.

## Limitations (stated in the report)

- EFO terms can pool related phenotypes (adjusted vs unadjusted, sex-specific),
  whose directions can legitimately differ; hence `REVIEW` rather than a fail,
  and the `reported_trait` filter.
- Catalog hits are curated top associations, so traits with few significant
  loci among the matrix SNPs may be `INSUFFICIENT`.
- Only the sign is checked, not the magnitude.
