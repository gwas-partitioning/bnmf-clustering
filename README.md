## Pipeline for GWAS clustering using Bayesian NMF

### Choose the set of variants to be clustered (choose_variants.R)
**ld_pruning**: LD-based pruning of the input variant set  
**count_traits_per_variant**: Assess the fraction of traits missing each variant of interest  
**find_variants_needing_proxies**: Determine which variants need proxies (allele considerations, missingness, etc.) 
**choose_proxies**: Select proxies for the necessary variants and output the final variant set for clustering

### Prepare the final z-score matrix (prep_bNMF.R)
**fetch_summary_stats**: Retrieve z-scores and sample sizes across all traits for the final variant set  
**prep_z_matrix**: Final trait filters, sample size adjustment, and creation of non-negative input matrix

### Run bNMF and summarize the results (run_bNMF.R)
**run_bNMF**: Run the bNMF procedure (over multiple iterations)  
**summarize_bNMF**: Summarize the results and create heatmaps for visualization

## To-do:
* Allele alignment issues
