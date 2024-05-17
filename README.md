## T2D Multi-ancestry Partitioned Polygenic Scores
Cluster weights are available from the [Smith, Deutsch et al Nature Medicine 2024](https://www.nature.com/articles/s41591-024-02865-3) paper in the "Smith_Deutsch_NatureMedicine_2024" folder. Partitioned polygenic scores (pPS) can be generated using the enclosed variant cluster weights. We have included weights for ancestry-specific and multi-ancestry clusters:

* In the weights files, the "Effect_Allele" column denotes the T2D risk-increasing allele.
* When generating the pPS, all genotypes should be aligned to this allele! 
* Weights have been provided in hg38, however a liftover map (hg19 to hg38) is included in each subfolder.

## Pipeline for GWAS clustering using Bayesian non-negative matrix factorization (bNMF)

The bNMF procedure, as applied here, is used to detect clusters of GWAS variants for some outcome of interest based on the associations of those variants with a set of additional traits. This pipeline includes pre-processing steps (such as quality control of variants and traits and the choice of proxy variants), preparation of the z-score matrix, clustering, and summarization of results.

**Important:** The current pipeline makes certain assumptions and uses some hard-coded filenames, including:

* "VAR_ID"s for GWAS and trait-specific summary statistics are in a specific format (CHR_POS_REF_ALT), with alleles aligned in a consistent way across traits (i.e. variant matching is performed using a simple string match).
* The variant reference file linking VAR_IDs to rsIDs is based on the VAR_ID format above, and points to a file available on the Broad Institute compute cluster. The VAR_ID to rsID pairs can be generated using the script **"generate_varid_to_rsid_map_file.R"**

### Requirments
The pipeline relies APIs for performing several linkage disequilibrium (LD)-based operations, including LD-pruning and proxy variant searches. The user must either acquire an personal token for LDlinkR or download the topLD API, depending on how they wish to perform these steps:

* **LDlinkR**: request token from https://ldlink.nih.gov/?tab=apiaccess
  * LDlinkR is used within the ld_prune and choose_proxies functions. Alternatives include performing position-based clumping using the snp_clump function, and using the topLD API option for choose_proxies
* **topLD**: download API to your project folder (instructions here: https://github.com/linnabrown/topld_api)
  * Note that topLD is currently not compatible with macOS. 

### 1. Choose the set of variants to be clustered (choose_variants.R)
**ld_prune**: LD-based pruning of the input variant set  
**snp_clip**: Alternative to ld_prune; uses chromosomal position to prune variant set
**count_traits_per_variant**: Assess the fraction of traits missing each variant of interest  
**find_variants_needing_proxies**: Determine which variants need proxies (allele considerations, missingness, etc.) 
**choose_proxies**: Select proxies for the necessary variants and output the final variant set for clustering

### 2. Prepare the final z-score matrix (prep_bNMF.R)
**fetch_summary_stats**: Retrieve z-scores and sample sizes across all traits for the final variant set 
**fill_missing_zcores**: Fill missing values in the variant-trait association matrix
**prep_z_matrix**: Final trait filters, sample size adjustment, and creation of non-negative input matrix (N (variants) x 2M (traits), with separate columns for positive trait associations (zero otherwise) and negative trait associations)

### 3. Run bNMF and summarize the results (run_bNMF.R)
**run_bNMF**: Run the bNMF procedure (over multiple iterations)  
**summarize_bNMF**: Summarize the results and create heatmaps for visualization

### 4. Further summarize and visualize results
**format_bNMF_results.Rmd**: Generates a HTML file which includes several plots and tables which summarize the cluster results

  * also included is a calculation for finding an optimal weight threshold for defining cluster-defining variants and traits

### Outputs
Most steps of the pipeline will print messages with details of the procedure. In addition, the following outputs will be written to the working directory.  

* no_proxies_found.txt: A list of variants that were excluded and for which no acceptable proxies were found.  
* run_summary.txt: A table listing the chosen K (# of clusters) and negative log-likelihood for each bNMF iteration.
* z_score_mat.rds: A binary R object containing the N x M z-score matrix after all preprocessing steps.
* L2EU.W.mat.K[]: The matrix of feature contributions to clusters for the K in question (one per K chosen in at least one iteration).
* L2EU.H.mat.K[]: The matrix of variant contributions to clusters for the K in question.
* W_plot_K[].pdf: A heatmap displaying feature contributions to clusters for the K in question.
* H_plot_K[].pdf: A heatmap displaying variant contributions to clusters for the K in question.

### Contributors
* Claire Kim (design and code)
* Kenny Westerman (design and code)
* Kirk Smith (code)
* Jaegil Kim (code)
* Marcin von Grotthuss (code)
* Miriam Udler (design and supervision)
