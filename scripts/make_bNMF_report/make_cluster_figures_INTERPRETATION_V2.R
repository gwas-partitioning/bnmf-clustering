# make_cluster_figures_INTERPRETATION_V2.R
# Plots and deterministic cluster-interpretation summaries for the bNMF report.
# Uses only the pipeline's own outputs (W, H, z-score matrices, nearest-gene
# labels); it makes no web requests and needs no annotation databases.

library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(magrittr)
library(data.table)

# Bioconductor annotation packages export a different select() generic.
# Keep data-frame pipelines deterministic throughout this sourced helper.
select <- dplyr::select
rename <- dplyr::rename

# Core utility functions
sanitize_label <- function(x, wrap = NULL, truncate = NULL, ascii = TRUE, ellipsis = "...") {
  x <- as.character(x)
  is_na <- is.na(x)
  
  if (any(!is_na)) {
    x[!is_na] <- str_squish(str_replace_all(x[!is_na], "[^[:print:]]", ""))
    if (ascii) x[!is_na] <- iconv(x[!is_na], from = "", to = "ASCII//TRANSLIT", sub = "")
  }
  
  if (!is.null(truncate) && truncate > 0L) {
    over_idx <- which(!is_na & nchar(x) > truncate)
    if (length(over_idx)) {
      keep <- max(0L, truncate - nchar(ellipsis))
      x[over_idx] <- paste0(str_sub(x[over_idx], 1, keep), ellipsis)
    }
  }
  
  if (!is.null(wrap) && wrap > 0L) {
    ok <- !is_na & nzchar(x)
    if (any(ok)) x[ok] <- str_wrap(x[ok], width = wrap)
  }
  
  x
}

# Map raw UKB-IDP / cardiometabolic trait names -> short, clinician-readable
# labels for circle-plot bar labels (base R only; strips _pos/_neg; handles
# bilateral *_plus_* averages).
shorten_trait_label <- function(x) {
  core <- function(s) {
    if (grepl("^IDP_dMRI_TBSS_", s)) {
      body <- sub("^IDP_dMRI_TBSS_", "", s); metric <- sub("_.*$", "", body); rest <- sub("^[^_]+_", "", body)
      h <- if (grepl("_R$|_R_", body)) " (R)" else if (grepl("_L$|_L_", body)) " (L)" else ""
      rest <- sub("_(L|R)$", "", rest); return(paste0(gsub("_", " ", rest), " ", metric, h))
    }
    if (grepl("^IDP_dMRI_ProbtrackX_", s)) {
      body <- sub("^IDP_dMRI_ProbtrackX_", "", s)
      h <- if (grepl("_r$|_R$", body)) " (R)" else if (grepl("_l$|_L$", body)) " (L)" else ""
      body <- sub("_(l|r|L|R)$", "", body); metric <- sub("_.*$", "", body); tract <- sub("^[^_]+_", "", body)
      return(paste0(gsub("_", " ", tract), " ", metric, h))
    }
    if (grepl("^IDP_T1_FAST_ROIs_", s)) {
      h <- if (grepl("^IDP_T1_FAST_ROIs_L_", s)) " (L)" else if (grepl("^IDP_T1_FAST_ROIs_R_", s)) " (R)" else ""
      body <- sub("^IDP_T1_FAST_ROIs_[LR]_", "", s); return(paste0(gsub("_", " ", body), " GM", h))
    }
    if (grepl("^IDP_T1_FIRST_", s)) {
      h <- if (grepl("_left_", s)) " (L)" else if (grepl("_right_", s)) " (R)" else ""
      body <- sub("^IDP_T1_FIRST_(left|right)_", "", s); body <- sub("_volume$", "", body)
      return(paste0(gsub("_", " ", body), " vol", h))
    }
    if (grepl("^IDP_T1_SIENAX_", s)) {
      body <- sub("^IDP_T1_SIENAX_", "", s); body <- sub("_?normalised_volume$", "", body); body <- sub("_volume$", "", body)
      return(paste0("Total ", gsub("_", " ", body), " vol"))
    }
    if (grepl("WMH", s)) return("WMH volume")
    if (grepl("^IDP_SWI_T2star_", s)) {
      h <- if (grepl("_left_", s)) " (L)" else if (grepl("_right_", s)) " (R)" else ""
      body <- sub("^IDP_SWI_T2star_(left|right)_", "", s); return(paste0("T2* ", gsub("_", " ", body), h))
    }
    if (grepl("^IDP_tfMRI", s)) return(gsub("_", " ", sub("^IDP_", "", s)))
    if (grepl("^DKTatlas_", s)) {
      h <- if (grepl("^DKTatlas_lh_", s)) " (L)" else if (grepl("^DKTatlas_rh_", s)) " (R)" else ""
      body <- sub("^DKTatlas_(lh|rh)_", "", s); meas <- sub(".*_", "", body); reg <- sub("_[^_]*$", "", body)
      return(paste0(gsub("_", " ", reg), " ", meas, h))
    }
    if (grepl("^a2009s_", s)) {
      h <- if (grepl("^a2009s_lh_", s)) " (L)" else if (grepl("^a2009s_rh_", s)) " (R)" else ""
      body <- sub("^a2009s_(lh|rh)_", "", s); meas <- sub(".*_", "", body); reg <- sub("_[^_]*$", "", body)
      reg <- sub("^(G&S|G|S|Lat|Pole)_", "", reg); reg <- gsub("[_-]", " ", reg)
      return(paste0(trimws(reg), " ", meas, h))
    }
    if (grepl("^NODEamps", s)) {
      ica <- sub("^NODEamps(25|100)_.*", "\\1", s); node <- as.integer(sub("^NODEamps(25|100)_", "", s))
      return(sprintf("rfMRI amp node%d (ICA%s)", node, ica))
    }
    if (grepl("^volume_", s)) {
      body <- sub("^volume_", "", s)
      h <- if (grepl("^Left-|^lh", body)) " (L)" else if (grepl("^Right-|^rh", body)) " (R)" else ""
      body <- sub("^(Left|Right)-", "", body); body <- sub("^(lh|rh)", "", body)
      body <- gsub("-", " ", body); body <- gsub("([a-z])([A-Z])", "\\1 \\2", body)
      body <- sub("[ ]?Vol$", "", body, ignore.case = TRUE); return(paste0(trimws(body), " vol", h))
    }
    gsub("[._]", " ", sub("^IDP_", "", s))
  }
  one <- function(s) {
    s <- sub("_(pos|neg)$", "", as.character(s))
    bilat <- grepl("_plus_", s); if (bilat) s <- sub("_plus_.*$", "", s)
    lab <- core(s); if (bilat) lab <- sub("\\s*\\((L|R)\\)$", " (L+R)", lab)
    lab
  }
  vapply(as.character(x), one, character(1), USE.NAMES = FALSE)
}

# CORRECTED: Enhanced data preparation function (removed is.null checks per user)
prepare_cluster_data_enhanced <- function(w_mat_long, h_mat_long, cluster_name, cutoff,
                                         max_rows = 30, pos_multiplier = 1.5, neg_multiplier = -1.5, empty_bar = 1) {
  
  # Subset variants
  snp_data <- w_mat_long %>%
    filter(weight >= cutoff & cluster == cluster_name) %>%
    mutate(group = "Variant") %>%
    dplyr::rename(name = variant)
  
  # FIXED: Enhanced trait processing with modern stringr functions
  trait_data <- h_mat_long %>%
    filter(weight >= cutoff & cluster == cluster_name) %>%
    mutate(
      # Use str_ends instead of endsWith to avoid recursion
      pos_neg = case_when(
        is.na(trait) ~ "Unknown",
        str_ends(trait, "pos") ~ "pos",
        str_ends(trait, "neg") ~ "neg", 
        TRUE ~ "Unknown"
      ),
      # FIXED: Safer case_when logic (removed is.null checks per user)
      group = case_when(
        pos_neg == "pos" & trait_type == "protein" ~ "Pos_Protein",
        pos_neg == "neg" & trait_type == "protein" ~ "Neg_Protein",
        pos_neg == "pos" & trait_type == "phenotype" ~ "Pos_Phenotype",
        pos_neg == "neg" & trait_type == "phenotype" ~ "Neg_Phenotype",
        pos_neg == "pos" ~ "Pos_Trait",  # Fallback
        pos_neg == "neg" ~ "Neg_Trait",  # Fallback
        TRUE ~ "Unknown_Trait"
      )
    ) %>%
    dplyr::rename(name = trait)
  
  # Combine data
  combined_data <- bind_rows(snp_data, trait_data) %>%
    mutate(name = as.character(name))
  
  # Apply row limit
  if(nrow(combined_data) > max_rows){
    combined_data <- combined_data %>% 
      slice_max(weight, n = max_rows) %>% 
      arrange(group, desc(weight))
  } else {
    combined_data <- combined_data %>% 
      group_by(group) %>% 
      arrange(group, desc(weight)) %>% 
      ungroup()
  }
  
  # FIXED: Safer value calculation
  combined_data <- combined_data %>%
    mutate(
      value = case_when(
        str_starts(group, "Neg_") ~ weight * neg_multiplier,
        TRUE ~ weight * pos_multiplier
      ),
      # FIXED: Use str_replace_all instead of gsub
      name = str_replace_all(name, "_pos|_neg", ""),
      name = str_replace_all(name, "\\.|_", " ")
    )
  
  # Add empty rows for spacing
  if (nrow(combined_data) > 0) {
    groups <- unique(combined_data$group)
    groups <- groups[!is.na(groups)]
    
    if (length(groups) > 0 && empty_bar > 0) {
      na_row <- combined_data[1, , drop = FALSE]
      na_row[1, ] <- NA
      
      empty_list <- map(groups, ~ {
        out <- na_row[rep(1, empty_bar), , drop = FALSE]
        out$group <- .x
        out
      })
      
      empty_df <- bind_rows(empty_list)
      combined_data <- bind_rows(combined_data, empty_df) %>%
        arrange(group, desc(weight))
    }
  }
  
  combined_data$id <- seq_len(nrow(combined_data))
  return(combined_data)
}

# Component selection function
select_components_for_circle_enhanced <- function(df, cutoff, n_each = 5, n_remaining = 20) {
  
  df <- df %>% filter(is.finite(weight), weight >= cutoff)
  if (nrow(df) == 0) return(df)
  
  # Safe selection with proper error handling
  top_traits <- df %>% 
    filter(comp_type == "TRAIT") %>%
    arrange(desc(weight)) %>% 
    slice_head(n = min(n_each, nrow(.)))
  
  top_snps <- df %>% 
    filter(comp_type == "SNP") %>%
    arrange(desc(weight)) %>% 
    slice_head(n = min(n_each, nrow(.)))
  
  picked_keys <- c(top_traits$key, top_snps$key)
  
  remainder <- df %>% 
    filter(!(key %in% picked_keys)) %>%
    arrange(desc(weight)) %>% 
    slice_head(n = min(n_remaining, nrow(.)))
  
  bind_rows(top_traits, top_snps, remainder)
}

# Layout building function
build_circle_layout_enhanced <- function(df, pos_multiplier = 1.5, neg_multiplier = -1.5, empty_bar = 1) {
  
  if (nrow(df) == 0) return(df)
  
  df <- df %>%
    mutate(
      value = case_when(
        str_starts(group, "Neg_") ~ weight * neg_multiplier,
        TRUE ~ weight * pos_multiplier
      ),
      name = str_replace_all(name, "_pos|_neg", ""),
      name = str_replace_all(name, "[._]", " ")
    )
  
  # CRITICAL CHANGE: Create a simplified Sort_Group for ordering
  df <- df %>% 
    mutate(
      Sort_Group = case_when(
        group == "Variant" ~ "A_Variant",         
        str_starts(group, "Neg_") ~ "B_Neg_Trait", 
        str_starts(group, "Pos_") ~ "C_Pos_Trait", 
        TRUE ~ "D_Other"
      )
    )
  
  # Arrange first by the simplified group, and second by weight magnitude
  df <- df %>% arrange(Sort_Group, desc(weight))
  
  # 💡 REVISED LOGIC FOR SEPARATORS BETWEEN SORT_GROUPS
  if (empty_bar > 0 && nrow(df) > 0) {
    
    # 1. Define the desired order: Variants -> Positives -> Negatives
    #    This order matches the visual appearance in your plot (Green -> Red -> Blue)
    desired_order <- c("A_Variant", "C_Pos_Trait", "B_Neg_Trait", "D_Other")
    
    # Filter the order to only include groups actually present in the data
    present_sort_groups <- intersect(desired_order, unique(df$Sort_Group))
    
    empty_rows_list <- list()
    
    for (i in seq_along(present_sort_groups)) {
      current_sort_group <- present_sort_groups[i]
      group_data <- df %>% filter(Sort_Group == current_sort_group)
      
      # 1. If it's NOT the first group (A_Variant), insert a separator BEFORE the data
      if (i > 1) {
        # Create dummy row structure
        dummy_row <- group_data[1, , drop = FALSE]
        dummy_row[,] <- NA
        # Use a neutral/unique Sort_Group for ordering and transparent fill
        dummy_row$Sort_Group <- paste0(current_sort_group, "_SEP") 
        dummy_row$group <- paste0(current_sort_group, "_SEP") # Match the transparent fill in plot_circular_bar_enhanced
        dummy_row$name <- NA 
        dummy_row$value <- 0 
        
        # Add separator BEFORE the group data
        empty_rows_list[[length(empty_rows_list) + 1]] <- dummy_row[rep(1, empty_bar), , drop = FALSE]
      }
      
      # 2. Add the actual data for the current group
      empty_rows_list[[length(empty_rows_list) + 1]] <- group_data
    }
    
    # Reconstruct the dataframe. The sorting logic is now baked into the list order.
    df <- bind_rows(empty_rows_list)
    
  } else {
    # If no empty bars, just keep the arrangement
    df <- df %>% arrange(Sort_Group, desc(weight)) 
  }
  
  df <- df %>% mutate(id = dplyr::row_number())
  df$label_txt <- sanitize_label(df$name, truncate = 16, wrap = NULL)
  
  return(df)
}

# Enhanced plotting function with protein vs phenotype colors
plot_circular_bar_enhanced <- function(df, cluster_name, k, output_dir) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  label_data <- df %>% filter(!is.na(value) & !is.na(name))
  
  if (nrow(label_data) == 0) {
    message("No data to plot for cluster ", cluster_name)
    return(NULL)
  }
  
  # Calculate angles
  number_of_bar <- nrow(df)
  label_data$angle <- 90 - 360 * (label_data$id - 0.5) / number_of_bar
  label_data$hjust <- ifelse(label_data$angle < -90, 1, 0)
  label_data$angle <- ifelse(label_data$angle < -90, label_data$angle + 180, label_data$angle)
  
  # Get value range
  val_range <- range(abs(df$value), na.rm = TRUE)
  val_max <- max(val_range)
  
  # Enhanced color scheme to distinguish proteins from phenotypes
  enhanced_colors <- c(
    "Variant" = "green3",
    "Neg_Protein" = "navy",        # Dark blue for proteins
    "Pos_Protein" = "darkred",     # Dark red for proteins
    "Neg_Phenotype" = "lightblue", # Light blue for phenotypes
    "Pos_Phenotype" = "lightcoral",# Light coral for phenotypes
    "Neg_Trait" = "blue3",         # Fallback
    "Pos_Trait" = "red3"           # Fallback
  )
  
  # Create plot
  p <- ggplot(df, aes(x = id, y = abs(value), fill = group)) +
    geom_bar(stat = "identity", alpha = 0.7, width = 0.8, na.rm = TRUE) +
    scale_fill_manual(values = enhanced_colors, na.value = "transparent") +
    geom_text(data = label_data, 
              aes(x = id, y = abs(value) + 0.05 * val_max,
                  label = label_txt, hjust = hjust, angle = angle),
              color = "black", fontface = "bold", alpha = 0.8, 
              size = 2.8, inherit.aes = FALSE) +
    ylim(-0.3 * val_max, 1.2 * val_max) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    coord_polar() +
    ggtitle(sprintf("Cluster %s (Proteins vs Phenotypes)", cluster_name)) +
    guides(fill = guide_legend(title = "Component Type", ncol = 3))
  
  # Save plots
  tryCatch({
    ggsave(file.path(output_dir, sprintf("circlePlot_enhanced_k%i_cluster%s.pdf", k, cluster_name)),
           p, device = "pdf", width = 9, height = 9)
    ggsave(file.path(output_dir, sprintf("circlePlot_enhanced_k%i_cluster%s.png", k, cluster_name)),
           p, device = "png", width = 9, height = 9, dpi = 300)
  }, error = function(e) {
    message("Error saving plots: ", e$message)
  })
  
  return(p)
}



# KEEP: Data preparation function
prepare_cluster_data <- function(w_mat_long, h_mat_long, cluster_name, cutoff,
                                 max_rows = 30, pos_multiplier = 1.5, neg_multiplier = -1.5, empty_bar = 1) {
  
  # Subset loci (variants) and traits based on the cutoff for the given cluster
  snp_data <- w_mat_long %>%
    filter(weight >= cutoff & cluster == cluster_name) %>%
    mutate(group = "Variant") %>%
    dplyr::rename(name = variant)
  
  # Subset traits and ensure group is character
  trait_data <- h_mat_long %>%
    filter(weight >= cutoff & cluster == cluster_name) %>%
    mutate(group = case_when(
      is.na(trait) ~ "Unknown_Trait",  # Handle NA traits
      endsWith(trait, "pos") ~ "Pos_Trait",
      endsWith(trait, "neg") ~ "Neg_Trait",
      TRUE ~ "Unknown_Trait"  # Handle traits without _pos or _neg
    )) %>%
    dplyr::rename(name = trait)
  
  # Combine SNP and trait data
  combined_data <- bind_rows(snp_data, trait_data) %>%
    mutate(name = as.character(name))
  
  # If there are more than max_rows, keep only the top max_rows by weight
  if(nrow(combined_data) > max_rows){
    combined_data <- combined_data %>% top_n(max_rows, wt = weight) %>% arrange(group, desc(weight))
  } else {
    combined_data <- combined_data %>% group_by(group) %>% arrange(group, desc(weight)) %>% ungroup()
  }
  
  # Scale values: multiply positive group by pos_multiplier, negative group by neg_multiplier
  combined_data <- combined_data %>%
    mutate(value = ifelse(group == "Neg_Trait", weight * neg_multiplier, weight * pos_multiplier),
           name = gsub("_neg|_pos", "", name),
           name = gsub("\\.|\\_", " ", name))
  
  # Add empty rows for spacing between groups (preserve column classes)
  groups <- levels(as.factor(combined_data$group))
  na_row <- combined_data[1, , drop = FALSE]  # template with correct classes
  na_row[1, ] <- NA                            # set all fields to NA (keeps classes)
  
  empty_list <- lapply(groups, function(g) {
    out <- na_row[rep(1, empty_bar), , drop = FALSE]
    out$group <- g
    out
  })
  empty_df <- dplyr::bind_rows(empty_list)
  
  combined_data <- dplyr::bind_rows(combined_data, empty_df) %>%
    dplyr::arrange(group, dplyr::desc(weight))
  
  # Assign sequential ID (will be used for plotting)
  combined_data$id <- seq_len(nrow(combined_data))
  
  return(combined_data)
}

# KEEP: Component selection function
select_components_for_circle <- function(df, cutoff, n_each = 5, n_remaining = 20) {
  
  df <- df %>% filter(is.finite(weight), weight >= cutoff)
  
  # top N per class (handles n()<N)
  top_traits <- df %>% filter(comp_type == "TRAIT") %>%
    arrange(desc(weight)) %>% slice_head(n = n_each)
  
  top_snps <- df %>% filter(comp_type == "SNP") %>%
    arrange(desc(weight)) %>% slice_head(n = n_each)
  
  picked_keys <- c(top_traits$key, top_snps$key)
  
  remainder <- df %>% filter(!(key %in% picked_keys)) %>%
    arrange(desc(weight)) %>% slice_head(n = n_remaining)
  
  bind_rows(top_traits, top_snps, remainder)
}

# KEEP: Simplified layout building function 
build_circle_layout <- function(df, pos_multiplier = 1.5, neg_multiplier = -1.5,
                                empty_bar = 1) {
  # Expect columns: name, comp_type in {"TRAIT","SNP"}, group in {"Pos_Trait","Neg_Trait","Variant"}, weight
  
  # value for bar height (sign by group)
  df <- df %>%
    mutate(
      value = dplyr::case_when(
        group == "Neg_Trait" ~ weight * neg_multiplier,
        TRUE                  ~ weight * pos_multiplier
      ),
      # clean labels
      name = gsub("_neg|_pos", "", name),
      name = gsub("[._]", " ", name)
    )
  
  # order within groups
  group_levels <- c("Neg_Trait", "Pos_Trait", "Variant")
  df$group <- factor(df$group, levels = intersect(group_levels, unique(df$group)))
  
  df <- df %>%
    arrange(group, desc(weight))
  
  # add empty bars as separators between present groups
  if (empty_bar > 0) {
    present_groups <- levels(droplevels(df$group))
    sep_rows <- lapply(present_groups, function(g) {
      if (sum(df$group == g) == 0) return(NULL)
      out <- df[rep(which(df$group == g)[1], empty_bar), , drop = FALSE]
      out[,] <- NA
      out$group <- g
      out
    })
    sep_rows <- do.call(rbind, sep_rows)
    df <- bind_rows(df, sep_rows) %>% arrange(group, desc(weight))
  }
  
  # new sequential id for plotting
  df <- df %>%
    mutate(id = dplyr::row_number())
  
  # Create clean labels
  df$label_txt <- sanitize_label(df$name, truncate = 16, wrap = NULL)
  
  df
}

# KEEP: Final plotting function - styled like ggplot documentation

plot_circular_bar <- function(df, cluster_name, k, output_dir) {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get label data for bars with actual values (not NA spacers)
  label_data <- df %>% filter(!is.na(value) & !is.na(name))
  
  if (nrow(label_data) == 0) {
    message("No data to plot for cluster ", cluster_name)
    return(NULL)
  }
  
  # Calculate angles and hjust for labels (like the example)
  number_of_bar <- nrow(df)
  label_data$angle <- 90 - 360 * (label_data$id - 0.5) / number_of_bar
  label_data$hjust <- ifelse(label_data$angle < -90, 1, 0)
  label_data$angle <- ifelse(label_data$angle < -90, label_data$angle + 180, label_data$angle)
  
  # Prepare base line data for group separators
  base_data <- df %>% 
    filter(!is.na(group)) %>%
    group_by(group) %>% 
    summarize(start = min(id, na.rm = TRUE), 
              end = max(id, na.rm = TRUE),
              .groups = 'drop') %>% 
    rowwise() %>% 
    mutate(title_pos = mean(c(start, end))) # CHANGE HERE: Store the numeric position
  
  # Get value range for grid lines
  val_range <- range(abs(df$value), na.rm = TRUE)
  val_max <- max(val_range)
  
  # Create grid line values (like 25, 50, 75% of max)
  grid_values <- c(val_max * 0.25, val_max * 0.5, val_max * 0.75)
  
  # Prepare grid data for circular grid lines
  if (nrow(base_data) > 0 && length(grid_values) > 0) {
    # Create a proper expanded grid
    grid_data <- expand_grid(base_data, grid_y = grid_values)
  } else {
    grid_data <- data.frame()
  }
  
  # Create the plot
  p <- ggplot(df, aes(x = id, y = abs(value), fill = group)) +
    
    # Add grid lines first (so bars are on top)
    {if(nrow(grid_data) > 0) {
      geom_segment(data = grid_data, 
                   aes(x = start, y = grid_y, xend = end, yend = grid_y), 
                   colour = "grey80", alpha = 0.8, linewidth = 0.3, 
                   inherit.aes = FALSE)
    }} +
    
    # Add grid value labels
    {if(length(grid_values) > 0) {
      annotate("text", 
               x = rep(max(df$id), length(grid_values)), 
               y = grid_values, 
               label = round(grid_values, 1),
               color = "grey60", size = 2.5, angle = 0, 
               fontface = "bold", hjust = 1)
    }} +
    
    # Main bars
    geom_bar(stat = "identity", alpha = 0.7, width = 0.8, na.rm = TRUE) +
    
    # Color scheme
    scale_fill_manual(values = c("Variant" = "green3",
                                 "Neg_Trait" = "blue3", 
                                 "Pos_Trait" = "red3"),
                      na.value = "transparent") +
    
    # Labels at the end of bars (no leader lines)
    geom_text(data = label_data, 
              aes(x = id, y = abs(value) + 0.05 * val_max,
                  label = label_txt, hjust = hjust, angle = angle),
              color = "black", fontface = "bold", alpha = 0.8, 
              size = 2.8, inherit.aes = FALSE) +
    
    # Base lines for group separators
    {if(nrow(base_data) > 0) {
      geom_segment(data = base_data, 
                   aes(x = start, y = -0.1 * val_max,
                       xend = end, yend = -0.1 * val_max), 
                   colour = "black", alpha = 0.8, linewidth = 0.6, 
                   inherit.aes = FALSE)
    }} +
    
    # # Group labels
    # {if(nrow(base_data) > 0) {
    #   geom_text(data = base_data, 
    #             aes(x = title_pos, y = -0.2 * val_max, label = group), # CHANGE HERE: use the new numeric 'title_pos'
    #             colour = "black", alpha = 0.8, size = 3.5, 
    #             fontface = "bold", inherit.aes = FALSE)
    # }} +
    
    # Set limits to give space for labels
    ylim(-0.3 * val_max, 1.2 * val_max) +
    
    # Styling
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-0.5, 4), "cm"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    coord_polar() +
    ggtitle(sprintf("Cluster %s", cluster_name))
  
  # Save plots
  ggsave(file.path(output_dir, sprintf("circlePlot_k%i_cluster%s.pdf", k, cluster_name)),
         p, device = cairo_pdf, width = 8, height = 8)
  ggsave(file.path(output_dir, sprintf("circlePlot_k%i_cluster%s.png", k, cluster_name)),
         p, device = ragg::agg_png, width = 8, height = 8, dpi = 300)
  
  p
}

# FINALIZED: Main function
make_circle_plots <- function(w_mat_long_wGene, h_mat_long, output_dir, 
                              cutoff = 0, n_each = 5, n_remaining = 20) {
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get cluster names from the data
  cluster_names <- unique(w_mat_long_wGene$cluster)
  k <- length(cluster_names)
  
  circular_plots <- list()
  
  for (current_cluster in cluster_names) {
    message(sprintf("Processing cluster %s...", current_cluster))
    
    # 1) Pull all components for this cluster (no cutoff here)
    df_full <- prepare_cluster_data(
      w_mat_long = w_mat_long_wGene,
      h_mat_long = h_mat_long,
      cluster_name = current_cluster,
      cutoff = -Inf,          # Don't drop yet
      max_rows = 100,
      empty_bar = 0           # No separators here; we'll add later
    )
    
    # 2) Define comp_type and key for selection
    df_full <- df_full %>%
      mutate(
        comp_type = dplyr::case_when(
          grepl("Trait", group) ~ "TRAIT",
          group == "Variant"    ~ "SNP",
          TRUE                  ~ NA_character_
        ),
        # Use 'name' as the unique key
        key = name
      )
    
    # 3) Select top components (5 + 5 + 20 >= cutoff)
    df_sel <- select_components_for_circle(
      df = df_full,
      cutoff = cutoff,
      n_each = n_each,
      n_remaining = n_remaining
    )
    
    if (!nrow(df_sel)) {
      message(sprintf("Cluster %s: no components above cutoff; skipping.", current_cluster))
      next
    }
    
    # 4) Build circular layout after selection
    df_layout <- build_circle_layout(
      df_sel,
      pos_multiplier = 1.5,
      neg_multiplier = -1.5,
      empty_bar = 1            # Group separators visible
    )
    
    # 5) Create and save the plot
    p <- plot_circular_bar(df_layout, current_cluster, k, output_dir)
    
    # Store the plot in the list
    circular_plots[[current_cluster]] <- p
  }
  
  message(sprintf("Created %d circular plots in %s", length(circular_plots), output_dir))
  return(circular_plots)
}


# Filtered Z-matrix function
make_filtered_zmat <- function(zmat_imputed, zmat_final, w_mat_long, h_mat_long, cutoff) {
  
  # Use the post-imputation matrix for plotting
  my_zmat <- zmat_imputed
  names(my_zmat) <- gsub("\\-|\\.", " ", names(my_zmat))
  colnames(zmat_final) <- gsub("\\-|\\.", " ", names(zmat_final))
  
  # Find the global Z-score range for consistent coloring
  global_z_max_abs <- quantile(abs(my_zmat), 0.9, na.rm = TRUE)
  
  # Prepare trait names from h_mat_long:
  h_mat_long_rev <- h_mat_long %>%
    mutate(base_trait = gsub("_pos|_neg", "", trait),
           base_trait = gsub("\\.|\\-", " ", base_trait))
  
  # Get cluster names
  cluster_names <- unique(w_mat_long$cluster)
  
  zscore_heatmaps <- list()
  for (cluster_name in cluster_names) {
    message("Processing cluster: ", cluster_name)
    
    # Get the SNPs and traits passing the cutoff for this cluster
    df_cluster_snps <- w_mat_long %>%
      filter(cluster == cluster_name & weight >= cutoff) %>%
      arrange(desc(weight))
    
    # Check if gene column exists, if not create it
    if (!"gene" %in% colnames(df_cluster_snps)) {
      df_cluster_snps <- df_cluster_snps %>%
        mutate(gene = make.unique(variant))
    } else {
      df_cluster_snps <- df_cluster_snps %>%
        mutate(gene = make.unique(gene))
    }
    
    df_cluster_traits <- h_mat_long_rev %>%
      filter(cluster == cluster_name & weight >= cutoff) %>%
      distinct(base_trait, .keep_all = TRUE) %>%
      arrange(desc(weight))
    
    if (nrow(df_cluster_snps) == 0 || nrow(df_cluster_traits) == 0) {
      message(sprintf("Skipping cluster %s: No SNPs or traits above cutoff.", cluster_name))
      next
    }
    
    # Check which variants exist in the matrix
    available_variants <- intersect(df_cluster_snps$variant, rownames(my_zmat))
    available_traits <- intersect(df_cluster_traits$base_trait, colnames(my_zmat))
    
    if (length(available_variants) == 0 || length(available_traits) == 0) {
      message(sprintf("Skipping cluster %s: No matching variants or traits found in matrix.", cluster_name))
      next
    }
    
    # Update the dataframes to only include available items
    df_cluster_snps <- df_cluster_snps %>%
      filter(variant %in% available_variants)
    
    df_cluster_traits <- df_cluster_traits %>%
      filter(base_trait %in% available_traits)
    
    # Subset the matrices using only available variants and traits
    my_zmat_subset <- my_zmat[df_cluster_snps$variant, df_cluster_traits$base_trait, drop = FALSE]
    zmat_final_subset <- zmat_final[df_cluster_snps$variant, df_cluster_traits$base_trait, drop = FALSE]
    
    # Check for empty subsets
    if (is.null(my_zmat_subset) || is.null(zmat_final_subset) ||
        nrow(my_zmat_subset) == 0 || ncol(my_zmat_subset) == 0) {
      message(sprintf("Skipping cluster %s: No data after subsetting.", cluster_name))
      next
    }
    
    # Apply gene names to rows
    rownames(my_zmat_subset) <- df_cluster_snps$gene
    rownames(zmat_final_subset) <- df_cluster_snps$gene
    
    # Convert to a long data frame
    df_heat <- as.data.frame(my_zmat_subset) %>%
      mutate(variant = rownames(.)) %>%
      pivot_longer(-variant, names_to = "trait", values_to = "zscore")
    
    # Convert imputed_flag to a long data frame
    imputed_flag <- is.na(zmat_final_subset)
    df_flag <- as.data.frame(imputed_flag) %>%
      mutate(variant = rownames(.)) %>%
      pivot_longer(-variant, names_to = "trait", values_to = "imputed")
    
    # Merge the two data frames and factor for ordering
    df_plot <- inner_join(df_heat, df_flag, by = c("variant", "trait")) %>%
      mutate(trait = fct_relevel(trait, df_cluster_traits$base_trait),
             variant = fct_rev(fct_relevel(variant, df_cluster_snps$gene)))
    
    # Create the heatmap
    ht <- ggplot(df_plot, aes(x = trait, y = variant, fill = zscore)) +
      geom_tile() +
      scale_fill_gradient2(
        low = "navy", mid = "white", high = "red", midpoint = 0,
        limits = c(-global_z_max_abs, global_z_max_abs),
        oob = scales::oob_squish
      ) +
      geom_text(data = filter(df_plot, imputed), aes(label = "*"), color = "black", size = 5) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_blank()) +
      ggtitle(paste("Z-scores for Cluster", cluster_name))
    
    zscore_heatmaps[[cluster_name]] <- ht
  }
  
  return(zscore_heatmaps)
}

# =============================================================================
# Interpretation helpers (V2)
# =============================================================================

effective_number <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x > 0]
  if (!length(x) || sum(x) == 0) return(0)
  p <- x / sum(x)
  1 / sum(p^2)
}

classify_trait_domain <- function(trait) {
  x <- tolower(as.character(trait))
  dplyr::case_when(
    grepl("idp_|dktatlas|a2009s|nodeamps|aparc|aseg_|amyg|hippsub|thalam|brain|cort|white.?matter|grey.?matter|gray.?matter", x) ~ "Brain imaging",
    grepl("bmi|body.?mass|waist|hip.?circumference|whr|adipos|body.?fat|fat.?mass|obesity", x) ~ "Adiposity",
    grepl("glucose|insulin|hba1c|glyc|homa|diabetes|t2d|proinsulin|ogtt", x) ~ "Glycemic",
    grepl("hdl|ldl|trig|cholesterol|lipid|apolipoprotein|\\bapo[ab]", x) ~ "Lipids",
    grepl("nafld|liver|hepatic|alanine|aspartate|\\balt\\b|\\bast\\b|\\bggt\\b", x) ~ "Liver",
    grepl("blood.?pressure|systolic|diastolic|hypertension|\\bsbp\\b|\\bdbp\\b", x) ~ "Blood pressure",
    grepl("crp|inflamm|cytokine|interleukin|tumou?r.?necrosis|\\btnf\\b", x) ~ "Inflammation",
    grepl("creatinine|egfr|kidney|renal|albuminuria", x) ~ "Renal",
    grepl("coronary|cardiovascular|myocard|stroke|cad\\b|heart", x) ~ "Cardiovascular",
    TRUE ~ "Other"
  )
}

resolve_trait_domain <- function(trait, trait_groups = NULL,
                                 fallback = c("grep", "other")) {
  fallback <- match.arg(fallback)

  if (!is.null(trait_groups) && all(c("trait", "Group") %in% names(trait_groups))) {
    group_lookup <- trait_groups %>%
      transmute(
        trait_key = tolower(trimws(as.character(trait))),
        Group = na_if(trimws(as.character(Group)), "")
      ) %>%
      filter(!is.na(Group)) %>%
      distinct(trait_key, .keep_all = TRUE)
    group <- group_lookup$Group[match(tolower(trimws(as.character(trait))), group_lookup$trait_key)]
    if (!is.na(group)) return(group)
  }

  if (fallback == "grep") classify_trait_domain(trait) else "Other"
}

make_signed_trait_summary <- function(h_mat_long, cutoff = 0,
                                      trait_groups = NULL,
                                      domain_fallback = c("grep", "other")) {
  domain_fallback <- match.arg(domain_fallback)
  required <- c("trait", "cluster", "weight")
  if (!all(required %in% names(h_mat_long))) {
    stop("h_mat_long must contain: ", paste(required, collapse = ", "))
  }

  h_mat_long %>%
    mutate(
      component = case_when(
        grepl("_pos$", trait) ~ "positive",
        grepl("_neg$", trait) ~ "negative",
        TRUE ~ "unsigned"
      ),
      base_trait = sub("_(pos|neg)$", "", trait)
    ) %>%
    group_by(cluster, base_trait) %>%
    summarise(
      positive_weight = sum(weight[component == "positive"], na.rm = TRUE),
      negative_weight = sum(weight[component == "negative"], na.rm = TRUE),
      unsigned_weight = sum(weight[component == "unsigned"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      magnitude = pmax(positive_weight, negative_weight, unsigned_weight),
      net_weight = positive_weight - negative_weight,
      direction = case_when(
        unsigned_weight == magnitude & unsigned_weight > 0 ~ "Unsigned",
        positive_weight > 0 & negative_weight > 0 &
          positive_weight >= cutoff & negative_weight >= cutoff ~ "Bidirectional",
        net_weight > 0 ~ "Positive",
        net_weight < 0 ~ "Negative",
        TRUE ~ "Neutral"
      ),
      display_trait = shorten_trait_label(base_trait),
      domain = vapply(
        base_trait,
        resolve_trait_domain,
        character(1),
        trait_groups = trait_groups,
        fallback = domain_fallback
      ),
      above_cutoff = magnitude >= cutoff
    )
}

make_variant_membership_summary <- function(w_mat_long, cutoff) {
  required <- c("variant", "cluster", "weight")
  if (!all(required %in% names(w_mat_long))) {
    stop("w_mat_long must contain: ", paste(required, collapse = ", "))
  }

  totals <- w_mat_long %>%
    group_by(variant) %>%
    summarise(
      total_weight = sum(weight, na.rm = TRUE),
      max_weight = max(weight, na.rm = TRUE),
      primary_cluster = cluster[which.max(weight)][1],
      clusters_above_cutoff = sum(weight >= cutoff, na.rm = TRUE),
      .groups = "drop"
    )

  w_mat_long %>%
    left_join(totals, by = "variant") %>%
    mutate(
      specificity = if_else(total_weight > 0, weight / total_weight, 0),
      membership = case_when(
        weight < cutoff ~ "Below cutoff",
        clusters_above_cutoff > 1 ~ "Shared",
        TRUE ~ "Cluster-specific"
      )
    )
}

make_trait_domain_summary <- function(signed_traits, cutoff) {
  signed_traits %>%
    filter(magnitude >= cutoff) %>%
    group_by(cluster, domain) %>%
    summarise(domain_weight = sum(magnitude), n_traits = n(), .groups = "drop") %>%
    group_by(cluster) %>%
    mutate(weight_fraction = domain_weight / sum(domain_weight)) %>%
    ungroup() %>%
    arrange(cluster, desc(domain_weight))
}

make_cluster_interpretation <- function(w_mat_long, h_mat_long, cutoff,
                                        variant_annotations = NULL,
                                        top_n = 5,
                                        trait_groups = NULL,
                                        domain_fallback = c("grep", "other"),
                                        trait_concentration = NULL) {
  domain_fallback <- match.arg(domain_fallback)
  memberships <- make_variant_membership_summary(w_mat_long, cutoff)
  traits <- make_signed_trait_summary(
    h_mat_long,
    cutoff,
    trait_groups = trait_groups,
    domain_fallback = domain_fallback
  )
  domains <- make_trait_domain_summary(traits, cutoff)

  if (!is.null(variant_annotations) && "variant" %in% names(variant_annotations)) {
    ann <- variant_annotations %>% distinct(variant, .keep_all = TRUE)
    memberships <- memberships %>% left_join(ann, by = "variant")
  }

  cluster_ids <- unique(as.character(w_mat_long$cluster))
  overview <- vector("list", length(cluster_ids))
  narratives <- setNames(vector("list", length(cluster_ids)), cluster_ids)
  top_traits <- setNames(vector("list", length(cluster_ids)), cluster_ids)
  top_loci <- setNames(vector("list", length(cluster_ids)), cluster_ids)
  domain_tables <- setNames(vector("list", length(cluster_ids)), cluster_ids)

  for (i in seq_along(cluster_ids)) {
    cl <- cluster_ids[i]
    loci_cl <- memberships %>%
      filter(cluster == cl, weight >= cutoff) %>%
      arrange(desc(weight), desc(specificity))
    traits_cl <- traits %>%
      filter(cluster == cl, magnitude >= cutoff) %>%
      arrange(desc(magnitude))
    domains_cl <- domains %>% filter(cluster == cl) %>% arrange(desc(domain_weight))

    locus_weights <- loci_cl$weight
    trait_weights <- traits_cl$magnitude
    top_locus_share <- if (length(locus_weights) && sum(locus_weights) > 0) max(locus_weights) / sum(locus_weights) else NA_real_
    top_trait_share <- if (length(trait_weights) && sum(trait_weights) > 0) max(trait_weights) / sum(trait_weights) else NA_real_
    effective_traits_value <- effective_number(trait_weights)
    if (!is.null(trait_concentration)) {
      concentration_row <- trait_concentration %>%
        filter(cluster == cl) %>%
        slice_head(n = 1)
      if (nrow(concentration_row)) {
        top_trait_share <- concentration_row$top_trait_share[1]
        effective_traits_value <- concentration_row$effective_traits[1]
      }
    }
    shared_fraction <- if (nrow(loci_cl)) mean(loci_cl$clusters_above_cutoff > 1) else NA_real_

    label_domains <- domains_cl$domain[domains_cl$domain != "Other"]
    label_domains <- head(label_domains, 2)
    label_traits <- head(traits_cl$display_trait, 2)
    provisional_label <- if (length(label_domains)) {
      paste0(paste(label_domains, collapse = " / "), " signature")
    } else if (length(label_traits)) {
      paste0(paste(label_traits, collapse = " / "), " signature")
    } else {
      "No components above cutoff"
    }

    direction_text <- if (nrow(traits_cl)) {
      paste0(head(traits_cl$display_trait, 3), " (", tolower(head(traits_cl$direction, 3)), ")", collapse = ", ")
    } else "no traits above cutoff"

    gene_col <- intersect(c("gene", "Nearest_Gene"), names(loci_cl))[1]
    gene_text <- "no gene annotations available"
    if (length(gene_col) && !is.na(gene_col)) {
      genes <- unique(na.omit(as.character(loci_cl[[gene_col]])))
      genes <- genes[nzchar(genes)]
      if (length(genes)) gene_text <- paste(head(genes, 3), collapse = ", ")
    }

    flags <- character()
    if (is.finite(top_locus_share) && top_locus_share >= 0.4) flags <- c(flags, "single-locus dominated")
    if (is.finite(top_trait_share) && top_trait_share >= 0.4) flags <- c(flags, "single-trait dominated")
    if (is.finite(shared_fraction) && shared_fraction >= 0.5) flags <- c(flags, "mostly shared loci")
    flag_text <- if (length(flags)) paste(flags, collapse = "; ") else "no major concentration flag"

    shared_text <- if (is.finite(shared_fraction)) {
      sprintf("%.0f%% of loci above cutoff are shared with another cluster; ", 100 * shared_fraction)
    } else {
      "shared-locus burden is not estimable; "
    }

    narratives[[cl]] <- paste0(
      "This provisional ", provisional_label, " is defined most strongly by ",
      direction_text, ". Leading annotated genes are ", gene_text, ". ",
      shared_text,
      flag_text, ". Directions refer to effects aligned to the primary-GWAS risk allele."
    )

    overview[[i]] <- tibble(
      cluster = cl,
      provisional_label = provisional_label,
      n_loci = nrow(loci_cl),
      n_traits = nrow(traits_cl),
      effective_loci = effective_number(locus_weights),
      effective_traits = effective_traits_value,
      shared_loci_percent = 100 * shared_fraction,
      top_locus_share_percent = 100 * top_locus_share,
      top_trait_share_percent = 100 * top_trait_share,
      interpretation_flag = flag_text
    )

    top_traits[[cl]] <- traits_cl %>%
      select(Trait = display_trait, Domain = domain, Direction = direction,
             Positive = positive_weight, Negative = negative_weight,
             Net = net_weight, Magnitude = magnitude) %>%
      slice_head(n = top_n)

    keep_locus_cols <- intersect(
      c("variant", "rsID", "VAR_ID", "gene", "Nearest_Gene",
        "weight", "specificity", "clusters_above_cutoff", "membership", "P_VALUE"),
      names(loci_cl)
    )
    top_loci[[cl]] <- loci_cl %>%
      select(all_of(keep_locus_cols)) %>%
      slice_head(n = top_n)
    domain_tables[[cl]] <- domains_cl
  }

  list(
    overview = bind_rows(overview),
    narratives = narratives,
    top_traits = top_traits,
    top_loci = top_loci,
    domains = domain_tables,
    signed_traits = traits,
    memberships = memberships
  )
}

make_cluster_overlap_matrix <- function(w_mat_long, cutoff) {
  clusters <- unique(as.character(w_mat_long$cluster))
  sets <- setNames(lapply(clusters, function(cl) {
    unique(w_mat_long$variant[w_mat_long$cluster == cl & w_mat_long$weight >= cutoff])
  }), clusters)

  tidyr::expand_grid(cluster_1 = clusters, cluster_2 = clusters) %>%
    rowwise() %>%
    mutate(
      intersection = length(intersect(sets[[cluster_1]], sets[[cluster_2]])),
      union = length(union(sets[[cluster_1]], sets[[cluster_2]])),
      jaccard = if_else(union > 0, intersection / union, 0)
    ) %>%
    ungroup()
}

plot_cluster_overlap <- function(overlap_df) {
  ggplot(overlap_df, aes(cluster_1, cluster_2, fill = jaccard)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", jaccard)), size = 3) +
    scale_fill_gradient(low = "white", high = "#2274CC", limits = c(0, 1)) +
    coord_equal() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
    labs(x = NULL, y = NULL, fill = "Jaccard", title = "Locus membership overlap")
}

compare_cluster_trait_profiles <- function(current_h_long, previous_h_long) {
  current <- make_signed_trait_summary(current_h_long, cutoff = -Inf) %>%
    select(cluster, base_trait, net_weight)
  previous <- make_signed_trait_summary(previous_h_long, cutoff = -Inf) %>%
    select(cluster, base_trait, net_weight)

  pairs <- tidyr::expand_grid(
    current_cluster = unique(current$cluster),
    previous_cluster = unique(previous$cluster)
  )

  scores <- pairs %>%
    rowwise() %>%
    mutate(cosine_similarity = {
      a <- current %>% filter(cluster == current_cluster) %>% select(base_trait, a = net_weight)
      b <- previous %>% filter(cluster == previous_cluster) %>% select(base_trait, b = net_weight)
      ab <- inner_join(a, b, by = "base_trait")
      denom <- sqrt(sum(ab$a^2)) * sqrt(sum(ab$b^2))
      if (nrow(ab) && denom > 0) sum(ab$a * ab$b) / denom else NA_real_
    }) %>%
    ungroup()

  scores %>%
    group_by(current_cluster) %>%
    slice_max(cosine_similarity, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(current_cluster)
}
