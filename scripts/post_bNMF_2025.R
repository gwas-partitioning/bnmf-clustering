# =============================================================================
# post_bNMF_2025.R
# =============================================================================
# Post-hoc analysis functions applied after bNMF clustering.
#
# Functions:
#   calculate_cutoff()  — Method M W elbow cutoff (top 0.2%, min n = 25)
#   calculate_h_core_traits() — H core traits: signed loading >= 50% of the cluster's
#                         largest single loading (scale-invariant; interpretation only)
#   do_post_analysis()  — Full post-hoc pipeline: load results, annotate nearest
#                         gene, export Excel workbook, liftover hg19→hg38,
#                         generate PRS score files
#
# Dependencies:
#   GenomicRanges, Homo.sapiens, openxlsx, rtracklayer (liftover), vroom
# =============================================================================

library(tidyverse)
library(magrittr)
library(rtracklayer)
library(vroom)
library(openxlsx)
library(GenomicRanges)
library(Homo.sapiens)

# post-hoc functions

calculate_cutoff <- function(w, h = NULL, output_dir = ".",
                             target = c("W", "H"),
                             upper_anchor_fraction = 0.002,
                             min_upper_anchor_n = 25L,
                             bottom_anchor_fraction = 0.80,
                             r2_warning = 0.40) {
  # Method M is used for W/SNP membership. H core traits are selected by
  # calculate_h_core_traits(), not by applying the W cutoff to H.
  target <- match.arg(target)
  if (target == "W") {
    cluster_cols <- grep("^X[0-9]+$", names(w), value = TRUE)
    if (!length(cluster_cols)) {
      stop("No X# cluster columns were found in W.")
    }
    factor_matrix <- as.matrix(w[, cluster_cols, drop = FALSE])
  } else {
    if (is.null(h)) stop("H must be supplied when target = 'H'.")
    numeric_cols <- vapply(h, is.numeric, logical(1))
    if (!any(numeric_cols)) stop("No numeric trait-loading columns were found in H.")
    factor_matrix <- as.matrix(h[, numeric_cols, drop = FALSE])
  }
  storage.mode(factor_matrix) <- "double"
  if (any(!is.finite(factor_matrix)) || any(factor_matrix < 0)) {
    stop("Method M requires finite, nonnegative ", target, " values.")
  }

  weight <- sort(as.numeric(factor_matrix), decreasing = TRUE)
  n_weights <- length(weight)
  if (!n_weights || max(weight) <= 0) {
    stop("Cannot calculate a cutoff from empty or all-zero ", target, " weights.")
  }

  upper_n_from_fraction <- as.integer(ceiling(
    upper_anchor_fraction * n_weights
  ))
  upper_n <- max(as.integer(min_upper_anchor_n), upper_n_from_fraction)
  if (upper_n >= n_weights) {
    stop("Method M upper anchor must be smaller than the pooled W vector.")
  }

  weights <- tibble(
    weight = weight,
    rank = seq_len(n_weights),
    x_scale = seq_len(n_weights) / n_weights,
    w_scale = weight / max(weight)
  )

  top_idx <- seq_len(upper_n)
  bottom_threshold <- as.numeric(quantile(
    weight,
    probs = bottom_anchor_fraction,
    names = FALSE,
    type = 7
  ))
  bottom_idx <- which(weight <= bottom_threshold)
  if (length(bottom_idx) < 2L) {
    stop("Method M bottom anchor contains fewer than two observations.")
  }

  upper_fit <- lm(w_scale ~ x_scale, data = weights[top_idx, , drop = FALSE])
  bottom_fit <- lm(w_scale ~ x_scale, data = weights[bottom_idx, , drop = FALSE])

  point_line_distance <- function(x, y, fit) {
    intercept <- unname(coef(fit)[1])
    slope <- unname(coef(fit)[2])
    abs(slope * x - y + intercept) / sqrt(slope^2 + 1)
  }

  weights$distance_to_upper <- point_line_distance(
    weights$x_scale, weights$w_scale, upper_fit
  )
  weights$distance_to_bottom <- point_line_distance(
    weights$x_scale, weights$w_scale, bottom_fit
  )
  weights$distance_difference <-
    weights$distance_to_upper - weights$distance_to_bottom

  cutoff_rank <- which(weights$distance_difference > 0)[1]
  if (is.na(cutoff_rank)) {
    stop("Method M reference-line distances never cross the observed ", target, " curve.")
  }

  cutoff <- weights$weight[cutoff_rank]
  upper_r2 <- summary(upper_fit)$r.squared
  upper_adjusted_r2 <- summary(upper_fit)$adj.r.squared
  bottom_r2 <- summary(bottom_fit)$r.squared
  bottom_adjusted_r2 <- summary(bottom_fit)$adj.r.squared

  diagnostics <- list(
    method = paste0("M_", tolower(target), "_only_fixed_upper_anchor"),
    method_label = paste0(
      "Method M: ", target, "-only elbow, top 0.2% anchor (min 25)"
    ),
    matrix_side = target,
    n_pooled_weights = n_weights,
    n_pooled_w_weights = if (target == "W") n_weights else NA_integer_,
    n_pooled_h_weights = if (target == "H") n_weights else NA_integer_,
    upper_anchor_fraction_requested = upper_anchor_fraction,
    upper_anchor_n_from_fraction = upper_n_from_fraction,
    minimum_upper_anchor_n = as.integer(min_upper_anchor_n),
    upper_anchor_n = upper_n,
    upper_anchor_floor_applied = upper_n > upper_n_from_fraction,
    upper_r2 = upper_r2,
    upper_adjusted_r2 = upper_adjusted_r2,
    bottom_anchor_fraction = bottom_anchor_fraction,
    bottom_anchor_n = length(bottom_idx),
    bottom_r2 = bottom_r2,
    bottom_adjusted_r2 = bottom_adjusted_r2,
    cutoff = cutoff,
    cutoff_rank = cutoff_rank,
    cutoff_rank_fraction = weights$x_scale[cutoff_rank],
    normalization = paste0("rank/N and ", target, "/max(", target, ")"),
    distance_rule = paste(
      "first ranked", target, "value whose perpendicular distance to the bottom",
      "reference line is smaller than its distance to the upper reference line"
    ),
    upper_intercept = unname(coef(upper_fit)[1]),
    upper_slope = unname(coef(upper_fit)[2]),
    bottom_intercept = unname(coef(bottom_fit)[1]),
    bottom_slope = unname(coef(bottom_fit)[2])
  )

  if (is.finite(upper_adjusted_r2) && upper_adjusted_r2 < r2_warning) {
    warning(
      "Method M ", target, " upper-anchor adjusted R2 = ",
      signif(upper_adjusted_r2, 3),
      " is below the diagnostic warning threshold of ",
      r2_warning,
      "."
    )
  }

  weights$group <- "Middle transition"
  weights$group[bottom_idx] <- "Bottom 80% anchor"
  weights$group[top_idx] <- "Top 0.2% anchor"
  weights$group <- factor(
    weights$group,
    levels = c("Top 0.2% anchor", "Middle transition", "Bottom 80% anchor")
  )
  highlight <- weights[cutoff_rank, , drop = FALSE]

  cutoff_plot <- ggplot(weights, aes(x_scale, w_scale)) +
    geom_point(aes(color = group), size = 0.75, alpha = 0.55) +
    geom_abline(
      intercept = coef(upper_fit)[1],
      slope = coef(upper_fit)[2],
      color = "#0072B2",
      linewidth = 1.0
    ) +
    geom_abline(
      intercept = coef(bottom_fit)[1],
      slope = coef(bottom_fit)[2],
      color = "#009E73",
      linewidth = 1.0
    ) +
    geom_vline(
      xintercept = max(weights$x_scale[top_idx]),
      color = "#0072B2",
      linewidth = 0.55,
      linetype = "dotted"
    ) +
    geom_point(
      data = highlight,
      shape = 21,
      size = 4.5,
      stroke = 1.0,
      fill = "#F0E442",
      color = "black"
    ) +
    annotate(
      "label",
      x = 0.98,
      y = 0.98,
      hjust = 1,
      vjust = 1,
      size = 3.4,
      label = sprintf(
        paste0(
          "Top 0.2%% anchor: n = %s\n",
          "Upper-line R² = %.3f (adjusted %.3f)\n",
          "Bottom 80%% anchor: n = %s\n",
          "Bottom-line R² = %.3f (adjusted %.3f)"
        ),
        format(upper_n, big.mark = ","),
        upper_r2,
        upper_adjusted_r2,
        format(length(bottom_idx), big.mark = ","),
        bottom_r2,
        bottom_adjusted_r2
      )
    ) +
    annotate(
      "label",
      x = min(0.96, highlight$x_scale + 0.08),
      y = min(0.90, highlight$w_scale + 0.08),
      hjust = 0,
      size = 3.2,
      label = sprintf(
        "Cutoff = %.4f\nrank = %s (%.2f%%)",
        cutoff,
        format(cutoff_rank, big.mark = ","),
        100 * highlight$x_scale
      )
    ) +
    scale_color_manual(values = c(
      "Top 0.2% anchor" = "#56B4E9",
      "Middle transition" = "grey55",
      "Bottom 80% anchor" = "#009E73"
    )) +
    coord_fixed(
      ratio = 1,
      xlim = c(-0.01, 1.02),
      ylim = c(-0.01, 1.02),
      clip = "on"
    ) +
    labs(
      title = paste0("Method M ", target, "-only elbow cutoff"),
      subtitle = paste0(
        "Upper reference: top 0.2% of pooled ", target,
        " weights (minimum n = ",
        min_upper_anchor_n,
        "); lower reference: bottom 80%"
      ),
      x = paste0("Normalized rank among pooled ", target, " weights"),
      y = paste0(target, " / maximum ", target),
      color = NULL,
      caption = paste(
        "Yellow point: first observed weight closer to the bottom reference",
        "line than to the upper reference line. Equal x/y coordinate scaling."
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      plot.caption = element_text(hjust = 0, size = 9)
    )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  figure_stem <- if (target == "W") {
    "bNMF_cutoff_calc_method_M"
  } else {
    "bNMF_cutoff_calc_method_M_H"
  }
  diagnostics_stem <- if (target == "W") {
    "bNMF_cutoff_diagnostics"
  } else {
    "bNMF_cutoff_diagnostics_H"
  }
  ggsave(
    file.path(output_dir, paste0(figure_stem, ".png")),
    cutoff_plot,
    width = 8.5,
    height = 8.5,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(output_dir, paste0(figure_stem, ".pdf")),
    cutoff_plot,
    width = 8.5,
    height = 8.5,
    units = "in",
    bg = "white"
  )

  # Preserve the historical W-figure filename used by the general HTML report.
  if (target == "W") {
    ggsave(
      file.path(output_dir, "bNMF_cutoff_calc.png"),
      cutoff_plot,
      width = 8.5,
      height = 8.5,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  }

  saveRDS(diagnostics, file.path(output_dir, paste0(diagnostics_stem, ".rds")))
  write.csv(
    as.data.frame(diagnostics, stringsAsFactors = FALSE),
    file.path(output_dir, paste0(diagnostics_stem, ".csv")),
    row.names = FALSE,
    quote = TRUE
  )

  cutoff_with_diagnostics <- cutoff
  attr(cutoff_with_diagnostics, "diagnostics") <- diagnostics
  cutoff_with_diagnostics
}


calculate_h_core_traits <- function(h, cluster_names,
                                    relative_max_threshold = 0.50) {
  # Select descriptive H "core traits" independently within each cluster,
  # thresholding the ORIGINAL signed H columns (trait_pos / trait_neg / unsigned):
  #
  #   1. H_max,k = the largest single loading in cluster k (any one column).
  #   2. T_k     = relative_max_threshold * H_max,k.
  #   3. Every column with H_kj >= T_k is selected.
  #   4. Selected columns are collapsed to base traits only AFTERWARD, for
  #      counts, PPI, and enrichment. A base trait with both its positive and
  #      negative column selected is labelled "bidirectional".
  #
  # Positive and negative columns are never combined before thresholding, so a
  # trait is selected only if one of its own loadings reaches T_k.
  #
  # The rule is invariant to component-wise NMF rescaling and is intended for
  # cluster interpretation, PPI, and enrichment, not for SNP/pPS membership.
  #
  # Returns:
  #   all_traits         one row per cluster x base trait (selected or not),
  #                      including full-H loading_energy (sum of squared
  #                      loadings) as a diagnostic that is NOT used for selection
  #   core_traits        the selected rows of all_traits (one row per base trait)
  #   core_signed_traits one row per selected signed column, for plots that
  #                      must show both signs of a bidirectional trait
  #   cluster_summary    per-cluster counts, membership-direction counts, and
  #                      full-H concentration diagnostics

  if (!is.numeric(relative_max_threshold) ||
      length(relative_max_threshold) != 1L ||
      !is.finite(relative_max_threshold) ||
      relative_max_threshold <= 0 || relative_max_threshold > 1) {
    stop("relative_max_threshold must be one number in (0, 1].")
  }
  if (nrow(h) != length(cluster_names)) {
    stop(
      "H has ", nrow(h), " rows but ", length(cluster_names),
      " cluster names were supplied."
    )
  }

  h <- h[, !(names(h) %in% c("cluster", "Cluster")), drop = FALSE]
  numeric_cols <- vapply(h, is.numeric, logical(1))
  if (!any(numeric_cols)) stop("No numeric trait-loading columns were found in H.")
  H <- as.matrix(h[, numeric_cols, drop = FALSE])
  storage.mode(H) <- "double"
  if (any(!is.finite(H)) || any(H < 0)) {
    stop("H core-trait selection requires finite, nonnegative H values.")
  }
  if (anyDuplicated(colnames(H))) {
    stop("H has duplicated trait column names; core-trait selection needs unique names.")
  }

  signed <- as.data.frame(H, check.names = FALSE) %>%
    mutate(cluster = cluster_names, .before = 1) %>%
    pivot_longer(-cluster, names_to = "signed_trait", values_to = "weight") %>%
    mutate(
      base_trait_id = stringr::str_remove(signed_trait, "_(pos|neg)$"),
      base_trait = base_trait_id %>%
        stringr::str_remove(stringr::fixed("protein..")) %>%
        stringr::str_remove(stringr::fixed("trait..")),
      explicit_trait_type = case_when(
        stringr::str_starts(base_trait_id, stringr::fixed("protein..")) ~ "protein",
        stringr::str_starts(base_trait_id, stringr::fixed("trait..")) ~ "phenotype",
        TRUE ~ NA_character_
      ),
      direction = case_when(
        stringr::str_ends(signed_trait, "_neg") ~ "negative",
        stringr::str_ends(signed_trait, "_pos") ~ "positive",
        TRUE ~ "unsigned"
      )
    ) %>%
    group_by(cluster) %>%
    mutate(
      maximum_loading = max(weight),
      threshold_loading = relative_max_threshold * maximum_loading,
      relative_to_maximum = if_else(maximum_loading > 0, weight / maximum_loading, 0),
      selected = maximum_loading > 0 & weight > 0 & weight >= threshold_loading
    ) %>%
    ungroup()

  all_traits <- signed %>%
    group_by(cluster, base_trait_id, base_trait, explicit_trait_type) %>%
    summarise(
      positive_loading = if (any(direction == "positive")) {
        sum(weight[direction == "positive"])
      } else {
        NA_real_
      },
      negative_loading = if (any(direction == "negative")) {
        sum(weight[direction == "negative"])
      } else {
        NA_real_
      },
      unsigned_loading = if (any(direction == "unsigned")) {
        sum(weight[direction == "unsigned"])
      } else {
        NA_real_
      },
      max_signed_loading = max(weight),
      dominant_signed_trait = signed_trait[which.max(weight)],
      dominant_direction = direction[which.max(weight)],
      loading_energy = sum(weight^2),
      positive_selected = any(selected & direction == "positive"),
      negative_selected = any(selected & direction == "negative"),
      unsigned_selected = any(selected & direction == "unsigned"),
      maximum_loading = dplyr::first(maximum_loading),
      threshold_loading = dplyr::first(threshold_loading),
      .groups = "drop"
    ) %>%
    mutate(
      selected = positive_selected | negative_selected | unsigned_selected,
      membership_direction = case_when(
        positive_selected & negative_selected ~ "bidirectional",
        positive_selected ~ "positive",
        negative_selected ~ "negative",
        unsigned_selected ~ "unsigned",
        TRUE ~ NA_character_
      ),
      relative_to_maximum = if_else(
        maximum_loading > 0,
        max_signed_loading / maximum_loading,
        0
      ),
      relative_max_threshold = relative_max_threshold
    ) %>%
    arrange(match(cluster, cluster_names), desc(max_signed_loading), base_trait_id)

  core_signed_traits <- signed %>%
    filter(selected) %>%
    left_join(
      all_traits %>% dplyr::select(cluster, base_trait_id, membership_direction),
      by = c("cluster", "base_trait_id")
    ) %>%
    dplyr::select(
      cluster, signed_trait, base_trait_id, base_trait, explicit_trait_type,
      direction, weight, relative_to_maximum, maximum_loading,
      threshold_loading, membership_direction
    ) %>%
    arrange(match(cluster, cluster_names), desc(weight), signed_trait)

  signed_summary <- core_signed_traits %>%
    group_by(cluster) %>%
    summarise(
      core_signed_columns = n(),
      boundary_loading = min(weight),
      .groups = "drop"
    )

  cluster_summary <- all_traits %>%
    group_by(cluster) %>%
    summarise(
      total_base_traits = n(),
      core_traits = sum(selected),
      positive_only = sum(membership_direction == "positive", na.rm = TRUE),
      negative_only = sum(membership_direction == "negative", na.rm = TRUE),
      bidirectional = sum(membership_direction == "bidirectional", na.rm = TRUE),
      unsigned = sum(membership_direction == "unsigned", na.rm = TRUE),
      maximum_loading = dplyr::first(maximum_loading),
      threshold_loading = dplyr::first(threshold_loading),
      # Full-H concentration: computed over ALL base traits, independent of
      # which traits were selected, so a one-trait core set is not by itself
      # evidence that a single trait drives the cluster.
      full_h_top_trait_energy_share = if (sum(loading_energy) > 0) {
        max(loading_energy) / sum(loading_energy)
      } else {
        NA_real_
      },
      full_h_effective_traits = if (sum(loading_energy) > 0) {
        sum(loading_energy)^2 / sum(loading_energy^2)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    left_join(signed_summary, by = "cluster") %>%
    mutate(core_signed_columns = tidyr::replace_na(core_signed_columns, 0L)) %>%
    dplyr::select(
      cluster, total_base_traits, core_traits, core_signed_columns,
      positive_only, negative_only, bidirectional, unsigned,
      maximum_loading, threshold_loading, boundary_loading,
      full_h_top_trait_energy_share, full_h_effective_traits
    ) %>%
    arrange(match(cluster, cluster_names))

  list(
    method = "H_signed_relative_to_cluster_maximum",
    method_label = paste0(
      "H core traits: signed loading >= ",
      scales::percent(relative_max_threshold, accuracy = 1),
      " of the cluster's largest single loading"
    ),
    relative_max_threshold = relative_max_threshold,
    all_traits = all_traits,
    core_traits = all_traits %>% filter(selected),
    core_signed_traits = core_signed_traits,
    cluster_summary = cluster_summary
  )
}


do_post_analysis <- function(main_dir,
                            my_chain,
                            df_gwas,
                            k              = NULL,
                            cluster_names  = NULL,
                            loci_file      = "query",
                            make_excel     = TRUE,
                            do_liftover    = TRUE,
                            existing_dir   = NULL) {
  # Full post-hoc pipeline: load bNMF outputs, annotate nearest genes,
  # create Excel workbook, perform hg19→hg38 liftover, and write PRS score files.
  #
  # Args:
  #   main_dir         : Output directory from summarize_bNMF()
  #   my_chain         : rtracklayer chain object for hg19→hg38 liftover
  #   df_gwas          : Aligned GWAS summary stats data frame (from Step 7)
  #   k                : Number of clusters (NULL = auto-select by majority vote)
  #   cluster_names    : Character vector of cluster labels (NULL = "X1","X2",...)
  #   loci_file        : "query" to annotate nearest gene; otherwise path to loci file
  #   make_excel       : Write per-cluster sorted weight Excel workbook
  #   do_liftover      : Perform hg19→hg38 coordinate liftover
  #   existing_dir     : Optional alternative directory for rsID_map.txt lookup
  
  # 1.) load data ---
  
  # find majority K from run summary
  df_run_summary <- fread(file.path(main_dir,"run_summary.txt"),
                          stringsAsFactors = F, data.table = F)
  
  k_counts <- df_run_summary %>%
    dplyr::count(K) %>%
    mutate(perc = n/nrow(df_run_summary))
  
  # ONLY calculate best K if user did not provide one
  if (is.null(k)) {
    k <- k_counts$K[which.max(k_counts$n)]
    message(sprintf("Automatically selected K=%d (majority vote)", k))
  } else {
    message(sprintf("Using manually requested K=%d", k))
  }
  
  if (is.null(cluster_names)) {
    cluster_names <- paste0("X",1:k)
  }
  
  w <- fread(sprintf("%s/L2EU.W.mat.%i.txt",main_dir, k),
             stringsAsFactors = FALSE, data.table = F) %>%
    dplyr::select(variant, all_of(cluster_names))
  n_variants <- nrow(w)
  
  h <- fread(sprintf("%s/L2EU.H.mat.%i.txt",main_dir,k),
             stringsAsFactors = FALSE, data.table = F) %>%
    dplyr::select(-any_of("cluster")) %>%
    rename_at(
      vars(starts_with("X2hr")), function(x) {gsub("X","",x) })
  
  # W defines SNP/pPS membership using Method M.
  cutoff_w_result <- calculate_cutoff(
    w, h, output_dir = main_dir, target = "W"
  )
  cutoff_diagnostics_w <- attr(cutoff_w_result, "diagnostics")
  cutoff_w <- as.numeric(cutoff_w_result)

  # H is used descriptively for cluster interpretation. Within each cluster the
  # original signed H columns (trait_pos / trait_neg / unsigned) are selected
  # when their own loading is at least 50% of that cluster's largest single
  # loading; selected columns are collapsed to base traits afterward.
  h_core_result <- calculate_h_core_traits(
    h = h,
    cluster_names = cluster_names,
    relative_max_threshold = 0.50
  )
  readr::write_csv(
    h_core_result$core_traits,
    file.path(main_dir, "H_core_traits_relative_to_cluster_maximum.csv")
  )
  readr::write_csv(
    h_core_result$core_signed_traits,
    file.path(main_dir, "H_core_signed_columns_relative_to_cluster_maximum.csv")
  )
  readr::write_csv(
    h_core_result$cluster_summary,
    file.path(main_dir, "H_core_trait_counts_by_cluster.csv")
  )

  # Backward-compatible aliases: historical callers interpret `cutoff` and
  # `cutoff_diagnostics` as the SNP/W threshold.
  cutoff <- cutoff_w
  cutoff_diagnostics <- cutoff_diagnostics_w
  print(sprintf(
    "Method M optimal W-only weight cutoff is %.4f (upper n=%d; upper R2=%.3f; bottom R2=%.3f).",
    cutoff_w,
    cutoff_diagnostics_w$upper_anchor_n,
    cutoff_diagnostics_w$upper_r2,
    cutoff_diagnostics_w$bottom_r2
  ))
  print(sprintf(
    paste0(
      "H interpretation selects signed columns with loading >= %.0f%% of the ",
      "cluster's largest single loading; core sizes range from %d to %d base ",
      "traits across clusters (%d bidirectional trait-cluster memberships)."
    ),
    100 * h_core_result$relative_max_threshold,
    min(h_core_result$cluster_summary$core_traits),
    max(h_core_result$cluster_summary$core_traits),
    sum(h_core_result$cluster_summary$bidirectional)
  ))

  # 2.) get nearest gene ---
  print("reading rsID map...")
  rsID_dir <- ifelse(is.null(existing_dir), main_dir, existing_dir)
  df_rsIDs <- fread(file.path(rsID_dir,"rsID_map.txt"),
                    data.table=F, stringsAsFactors=F)
  
  variants_have_prefix <- grepl("chr", w$variant[1])
  
  if (!'variant' %in% names(df_rsIDs)) {
    if (variants_have_prefix) {
      df_rsIDs <- df_rsIDs %>%
        separate(VAR_ID, into=c('CHR','POS','REF','ALT'),sep='_',remove = F) %>%
        mutate(variant = paste0('chr',CHR,":",POS))
    } else {
      df_rsIDs <- df_rsIDs %>%
        separate(VAR_ID, into=c('CHR','POS','REF','ALT'),sep='_',remove = F) %>%
        mutate(variant = paste0(CHR,":",POS))
    }
    
  }
  
  # IMPROVED VERSION - Handle variant filtering throughout the gene annotation pipeline
  
  print("Making excel file...")
  if (make_excel==T) {
    print("Making excel file...")
    
    if (loci_file=="query") {
      # Gene annotation logic with robust parsing
      geneRanges <- function(db, column="ENTREZID"){
        g <- genes(db, columns=column)
        col <- mcols(g)[[column]]
        genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
        mcols(genes)[[column]] <- as.character(unlist(col))
        genes
      }
      
      # Suppress the gene dropping messages
      suppressMessages({
        gns = geneRanges(Homo.sapiens, column="SYMBOL")
        df_entrez_data <- geneRanges(Homo.sapiens, column="ENTREZID")
      })
      
      df_gns <- data.frame(gns)
      gtf.gene <- df_gns %>%
        mutate(chr = gsub("chr","",seqnames)) %>%
        dplyr::select(chr, start, end, SYMBOL)
      
      df_entrez <- data.frame(df_entrez_data) %>%
        mutate(chr = gsub("chr","",seqnames)) %>%
        mutate(ChrPos = paste(chr,start,sep=":")) %>%
        dplyr::select(ChrPos, ENTREZID)
      
      range2GRanges <- function(df) {
        require(GenomicRanges)
        require(IRanges)
        gr <- GenomicRanges::GRanges(
          seqnames = df[,1],
          ranges=IRanges(start = df[,2], end = df[,3])
        )
        return(gr)
      }
      
      # ROBUST VARIANT PARSING
      snps <- w$variant
      message(sprintf("Annotating %d variants with nearest gene...", length(snps)))

      # SAFE PARSING - handles any format
      ranges_df <- data.frame(
        chr = character(length(snps)),
        start = numeric(length(snps)),
        end = numeric(length(snps)),
        original_variant = snps,  # Keep track of original variant names
        stringsAsFactors = FALSE
      )
      
      for (i in seq_along(snps)) {
        snp <- snps[i]
        
        tryCatch({
          if (grepl("^chr", snp)) {
            # Format: chr1:12345
            clean_snp <- gsub("^chr", "", snp)
            parts <- strsplit(clean_snp, ":")[[1]]
            ranges_df$chr[i] <- parts[1]
            ranges_df$start[i] <- as.numeric(parts[2])
            ranges_df$end[i] <- as.numeric(parts[2])
            
          } else if (grepl(":", snp)) {
            # Format: 1:12345
            parts <- strsplit(snp, ":")[[1]]
            ranges_df$chr[i] <- parts[1]
            ranges_df$start[i] <- as.numeric(parts[2])
            ranges_df$end[i] <- as.numeric(parts[2])
            
          } else if (grepl("_", snp)) {
            # Format: 1_12345_A_T
            parts <- strsplit(snp, "_")[[1]]
            ranges_df$chr[i] <- parts[1]
            ranges_df$start[i] <- as.numeric(parts[2])
            ranges_df$end[i] <- as.numeric(parts[2])
            
          } else {
            # Unknown format - skip
            ranges_df$chr[i] <- NA
            ranges_df$start[i] <- NA
            ranges_df$end[i] <- NA
          }
        }, error = function(e) {
          # Skip any problematic variants
          ranges_df$chr[i] <- NA
          ranges_df$start[i] <- NA
          ranges_df$end[i] <- NA
        })
      }
      
      # Remove failed parses
      valid_rows <- !is.na(ranges_df$chr) & !is.na(ranges_df$start) & !is.na(ranges_df$end)
      ranges_df_clean <- ranges_df[valid_rows, ]
      snps_clean <- ranges_df_clean$original_variant  # Use original variant names
      
      message(sprintf("  Parsed %d/%d variants successfully", nrow(ranges_df_clean), length(snps)))
      
      if (nrow(ranges_df_clean) == 0) {
        warning("No variants could be parsed for gene annotation - using variant IDs as gene names.")
        w_wLoci <- w %>%
          mutate(gene = variant)
      } else {
        
        # Create GRanges object with error checking
        snps.granges <- GenomicRanges::GRanges(
          seqnames = ranges_df_clean$chr,
          ranges = IRanges::IRanges(start = ranges_df_clean$start, end = ranges_df_clean$end)
        )
        names(snps.granges) <- snps_clean
        
        message(sprintf("  Created GRanges with %d variants", length(snps.granges)))
        
        # Convert genes to GRanges
        gtf.granges <- GenomicRanges::GRanges(
          seqnames = gtf.gene$chr,
          ranges = IRanges::IRanges(start = as.numeric(gtf.gene$start), end = as.numeric(gtf.gene$end))
        )
        names(gtf.granges) <- gtf.gene$SYMBOL
        
        # Find nearest genes
        hits <- GenomicRanges::nearest(snps.granges, gtf.granges)
        
        # Prepare gene annotation data
        df_gns_final <- df_gns %>%
          mutate(chr = gsub("chr","",seqnames)) %>%
          mutate(ChrPos = paste(chr,start,sep=":")) %>%
          dplyr::select(chr, start, end, ChrPos, SYMBOL) %>%
          merge(df_entrez, by="ChrPos")
        
        # Create hits data frame
        df_hits <- data.frame(
          gene = names(gtf.granges)[hits],
          variant = names(snps.granges),
          stringsAsFactors = FALSE
        ) %>%
          merge(df_gns_final[,c("SYMBOL","ENTREZID")], by.x="gene", by.y="SYMBOL", all.x=TRUE) %>%
          filter(!duplicated(variant))
        
        message(sprintf("  Gene annotation completed for %d variants", nrow(df_hits)))
        
        # IMPORTANT: Handle variants that couldn't be parsed
        # Include all original variants, with gene annotation where available
        w_wLoci <- w %>%
          left_join(df_hits, by="variant") %>%
          mutate(gene = dplyr::coalesce(gene, variant))  # Use variant name if no gene found
        
        message(sprintf("  Final annotation: %d variants", nrow(w_wLoci)))
      }
      
    } else {
      w_wLoci <- w %>%
        mutate(gene=variant)
    }
    
    if (exists("df_rsIDs") && nrow(df_rsIDs) > 0) {
      w_excel_ready <- w_wLoci %>%
        inner_join(df_rsIDs, by = "variant")
      n_missing <- nrow(w_wLoci) - nrow(w_excel_ready)
      if (n_missing > 0) {
        warning(sprintf("%d variants could not be matched to rsID map and will be dropped from Excel output.", n_missing))
      }
    } else {
      message("No rsID mapping available - using variant names as identifiers.")
      w_excel_ready <- w_wLoci %>%
        mutate(VAR_ID = variant, rsID = variant)
    }

    # Rename nearest gene column
    w_excel_ready <- w_excel_ready %>%
      dplyr::rename(Nearest_Gene = gene)

    # Make Excel with improved error handling
    tryCatch({
      h.exc <- t(h) %>%
        data.frame() %>%
        mutate(trait = row.names(.)) %>%
        relocate(trait)
      
      w.exc <- w_excel_ready %>%
        dplyr::select(-any_of(c('ENTREZID', 'variant'))) %>%
        relocate(any_of(c("VAR_ID", "rsID", "Nearest_Gene")))

      # Determine last ID column safely
      id_cols <- intersect(c("VAR_ID", "rsID", "Nearest_Gene"), colnames(w.exc))
      last_id_col <- length(id_cols)
      
      posStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
      
      wb <- openxlsx::createWorkbook()
      for (i in 1:length(cluster_names)){
        w_temp <- w.exc %>%
          dplyr::arrange(desc(!!as.symbol(cluster_names[i])))
        h_temp <- h.exc %>%
          arrange(desc(!!as.symbol(cluster_names[i])))
        
        cur_sheet = sprintf("cluster%i",i)
        addWorksheet(wb, cur_sheet)
        start_col1 <- 1
        start_col2 <- ncol(w_temp)+3
        
        writeData(wb=wb, sheet=cur_sheet, x=w_temp, startCol = start_col1, startRow = 1,
                  rowNames=F,colNames = T)
        writeData(wb=wb, sheet=cur_sheet, x=h_temp, startCol = start_col2, startRow = 1,
                  rowNames=F,colNames = T)
        
        # Safe conditional formatting
        if (last_id_col > 0 && last_id_col + i <= ncol(w_temp)) {
          conditionalFormatting(wb, cur_sheet,
                                cols = last_id_col + i, rows = 2:(nrow(w_temp)+1),
                                rule = sprintf("%s2>=%.7f",
                                               openxlsx::int2col(last_id_col + i), cutoff_w),
                                style = posStyle)
        }
        
        if (start_col2 + i <= ncol(h_temp) + start_col2) {
          selected_signed_traits <- h_core_result$core_signed_traits %>%
            filter(cluster == cluster_names[i]) %>%
            pull(signed_trait)
          selected_h_rows <- which(h_temp$trait %in% selected_signed_traits) + 1L
          if (length(selected_h_rows)) {
            addStyle(
              wb,
              sheet = cur_sheet,
              style = posStyle,
              rows = selected_h_rows,
              cols = start_col2 + i,
              gridExpand = FALSE,
              stack = TRUE
            )
          }
        }
      }

      core_trait_export <- h_core_result$core_traits %>%
        dplyr::select(
          cluster, base_trait, explicit_trait_type, membership_direction,
          positive_loading, negative_loading, max_signed_loading,
          relative_to_maximum, maximum_loading, threshold_loading
        )
      addWorksheet(wb, "H_core_traits")
      writeData(
        wb, "H_core_traits", core_trait_export,
        startCol = 1, startRow = 1, rowNames = FALSE, colNames = TRUE
      )
      addWorksheet(wb, "H_core_summary")
      writeData(
        wb, "H_core_summary", h_core_result$cluster_summary,
        startCol = 1, startRow = 1, rowNames = FALSE, colNames = TRUE
      )
      
      fp <- file.path(main_dir, sprintf("sorted_cluster_weights_K%i_rev.xlsx",k))
      openxlsx::saveWorkbook(wb, fp, overwrite = T)
      
      cat(sprintf("Excel file saved: %s\n", fp))
      
    }, error = function(e) {
      cat("Error creating Excel file:", e$message, "\n")
      cat("Continuing without Excel output...\n")
    })
  }  
  
  # 3.) hg19-to-hg38 liftover --
  
  if (do_liftover==T) {
    print("Doing liftover...")
    if (is.character(my_chain)) my_chain <- rtracklayer::import.chain(my_chain)
    
    df_weights <- w_wLoci %>%
      inner_join(df_rsIDs, by=c('variant')) %>%
      # dplyr::select(VAR_ID, rsID, gene) %>%
      separate(VAR_ID, into=c("CHR", "POS_hg19", "REF", "ALT"), sep = "_", remove = F) %>%
      mutate(ChrPos = paste(CHR, POS_hg19, sep = ":")) %>%
      mutate(group = seq.int(nrow(.))) %>%
      mutate_at(.vars = c('CHR','POS_hg19'),as.integer)
    
    gr <- GRanges(seqnames = paste0("chr",df_weights$CHR),
                  strand = "*",
                  ranges = IRanges(start = df_weights$POS_hg19,
                                   width = 1))
    
    df_liftover <- as.data.frame(rtracklayer::liftOver(x = gr, chain = my_chain)) %>%
      dplyr::select(group, seqnames, Pos_hg38=start) %>%
      mutate(ChrPos_hg38=paste(seqnames, Pos_hg38, sep=":"))
    
    missing <- df_weights %>%
      filter(!group %in% df_liftover$group)
    if (nrow(missing)>0) {
      print("Variants not found in liftover...")
      print(missing$variant)
    }
    
    # merge weights with hg38 SNPs using rsID
    #   new file should have VAR_ID_hg38, rsID, Effect_Allele, BETA and the cluster columns (X1, X2, etc.)
    input_snps38 <- df_weights %>%
      inner_join(df_liftover[, c("group","Pos_hg38")], by="group") %>%
      mutate(VAR_ID_hg38 = paste(CHR, Pos_hg38, REF, ALT, sep = "_")) %>%
      inner_join(df_gwas[,c("ChrPos","Risk_Allele","BETA")], by="ChrPos")
    
    
    liftover_map <- input_snps38 %>%
      dplyr::select(VAR_ID_hg19=VAR_ID, Risk_Allele, rsID, VAR_ID_hg38)
    print("Writing liftover map to directory...")
    write_delim(liftover_map, file.path(main_dir,"liftover_map.txt"),
                delim = "\t",col_names = T)
    
    
    # make weight files
    input_snps38.1 <- input_snps38 %>%
      dplyr::select(VAR_ID_hg38, rsID, Risk_Allele, BETA, starts_with("X"))
    
    
    #---- CREATE WEIGHTS FILE FOR EACH CLUSTER-----
    
    
    cluster_weights <- input_snps38.1 %>%
      dplyr::select(VAR_ID_hg38, starts_with("X")) %>%
      column_to_rownames(var="VAR_ID_hg38")
    
    
    cluster_info = data.frame(weight=numeric(),VAR_ID_hg38=character(),cluster=character())
    
    for(i in 1:ncol(cluster_weights)){
      w_tmp = cluster_weights[i]
      names(w_tmp) = c("Weight")
      w_tmp$VAR_ID_hg38 = rownames(w_tmp)
      w_tmp = w_tmp[w_tmp$Weight>=cutoff_w,]
      w_tmp$cluster=cluster_names[i]
      cluster_info = rbind(cluster_info,w_tmp)
    }
    
    
    # for PRS, need following columns:
    # Chr, Pos, Ref, Alt, RSID, Effect_Allele, Weight (Effect_Allele=Risk_Allele here)
    #   Ex: 11, 92673828, C, T, rs1387153, T, 4.658804585
    input_df <- input_snps38 %>%
      merge(cluster_info, by="VAR_ID_hg38") %>%
      separate(VAR_ID_hg38, into=c("Chr","Pos","Ref","Alt"),sep="_",remove = F) %>%
      arrange(as.integer(Chr), Pos) %>%
      dplyr::select(Chr, Pos, Ref, Alt, RSID=rsID, Risk_Allele, Weight, cluster)
    # check no issues with mapping
    sum(is.na(input_df))
    
    # create file with snps and weights for each cluster
    print("Writing score files...")
    dir.create(file.path(main_dir,"score_info"))
    for (c in cluster_names) {
      df <- input_df %>%
        subset(cluster==c) %>%
        dplyr::select(-c(cluster))
      print(nrow(df))
      write.csv(x=df,
                file=file.path(main_dir,"score_info",sprintf("%s.csv", c)),
                row.names=F, quote=F)
    }
    
    # repeat for all SNPs (for total PRS)
    all_snps <- input_snps38 %>%
      separate(VAR_ID_hg38, into=c("Chr","Pos","Ref","Alt"),sep="_",remove = F) %>%
      arrange(as.integer(Chr), Pos) %>%
      mutate(Weight=abs(as.numeric(BETA))) %>%
      dplyr::select(Chr, Pos, Ref, Alt, RSID=rsID, Risk_Allele, Weight)
    write.csv(all_snps,
              file=file.path(main_dir,"score_info","Total_GRS.csv"),
              row.names=F, quote=F)
    
    print("Writing BCFtools input file...")
    for_VCF <- input_snps38 %>%
      # filter(VAR_ID_hg38 %in% cluster_info$VAR_ID_hg38) %>%
      separate(VAR_ID_hg38, into=c("Chr","Pos","Ref","Alt"),sep="_",remove = F) %>%
      arrange(as.integer(Chr), as.integer(Pos)) %>%
      mutate(Chr = paste0("chr",Chr), start = Pos, end = Pos) %>%
      dplyr::select(Chr, start, end)
    write_delim(for_VCF, file.path(main_dir,"bcftools_input.txt"),
                delim = "\t",col_names = F)
  }
  
  
  return(list(
    df_run_summary = df_run_summary,
    k = k,
    cluster_names = cluster_names,
    w = w_wLoci,
    h = h,
    cutoff = cutoff,
    cutoff_w = cutoff_w,
    cutoff_diagnostics = cutoff_diagnostics,
    cutoff_diagnostics_w = cutoff_diagnostics_w,
    h_selection_method = h_core_result$method,
    h_selection_method_label = h_core_result$method_label,
    h_relative_max_threshold = h_core_result$relative_max_threshold,
    h_core_traits = h_core_result$core_traits,
    h_core_signed_traits = h_core_result$core_signed_traits,
    h_trait_loadings = h_core_result$all_traits,
    h_core_cluster_summary = h_core_result$cluster_summary
  ))
}
