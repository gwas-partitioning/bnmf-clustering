1_proxy_prep.R: Prepare dictionary of possible proxies for primary index variants.
* Collect index variants
* Retrieve all potential proxies for each index variant (to be run remotely on HPC)
* Filter out potential proxies in low LD with associated index variant (<0.6) and strand-ambigouous variants.

2_gather_sumstats.R: Collect and preprocess summary statistics for all traits of interest.
* Collect rsID, A1, A2, z-score, and GWAS samples size.
* Include index variants and all potential proxies.

3_process_sumstats.Rmd: Prepare summary statistics for clustering.
* Assign proxies to variants that are strand-ambiguous or missing for >20% of traits (prioritize variants in order of LD with index variant).
* Prune away traits with no nominally-significant (p < 0.05) variants.
* Flip z-score signs as necessary to align to the trait-increasing alleles.
* Adjust z-scores for GWAS samples size: z_adj = z / sqrt(N).
* Replace remaining missing z-scores with zero.
* Augment the summary statistic matrix (traits x variants) to include separate positively- and negatively-associated rows for each trait (with the remaining set of variants set to zero).
* Visualize the input matrix in various ways (basic heatmap, trait correlation heatmap, etc.).
