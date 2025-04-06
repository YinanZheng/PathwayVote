#' @importFrom future plan
#' @importFrom furrr future_map furrr_options
#' @importFrom parallelly availableCores
#' @importFrom ReactomePA enrichPathway
#' @importFrom clusterProfiler enrichGO enrichKEGG setReadable
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @import dplyr

#' @title Filter eQTM Data by Statistics and Distance
#'
#' @description Filters eQTM data based on statistics and distance thresholds.
#' @param eQTM An eQTM object.
#' @param stat_threshold Numeric, the minimum threshold for the statistics value.
#' @param distance Numeric, the maximum distance threshold.
#' @return A list with filtered entrez IDs and p-values.
#' @export
#' @name gene_filter
gene_filter <- function(eQTM, stat_threshold, distance) {
  if (!inherits(eQTM, "eQTM")) {
    stop("eQTM must be an eQTM object")
  }
  eQTM_data <- getData(eQTM)

  # 通用过滤逻辑：对所有支持的 statistics_type 都取绝对值
  eQTM_filtered <- eQTM_data[abs(eQTM_data$statistics) >= stat_threshold & eQTM_data$distance <= distance, ]

  return(list(entrez = unique(na.omit(eQTM_filtered$entrez)), p_values = eQTM_filtered$p_value))
}

#' @title Filter eQTM Gene Lists with Voting Parameters
#'
#' @description Generate gene lists from an eQTM object using thresholds on statistics and distance.
#' @param eQTM An eQTM object.
#' @param stat_grid A vector of statistic thresholds.
#' @param distance_grid A vector of distance thresholds.
#' @param overlap_threshold Proportion threshold for keeping gene lists.
#' @param verbose Logical, whether to show messages.
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{gene_lists}: A list of gene vectors.
#'     \item \code{cpg_gene_maps}: A list of CpG → gene mappings.
#'     \item \code{params}: A list of parameter combinations used.
#'   }
#' @export
filter_gene_lists <- function(eQTM, stat_grid, distance_grid, overlap_threshold = 0.5, verbose = FALSE) {
  if (!inherits(eQTM, "eQTM")) stop("Input must be an eQTM object.")
  data <- getData(eQTM)

  gene_lists <- list()
  cpg_gene_maps <- list()
  params <- list()

  for (r in stat_grid) {
    for (d in distance_grid) {
      subset <- data[abs(data$statistics) > r & data$distance < d, ]
      if (nrow(subset) == 0) next

      map <- split(subset$entrez, subset$cpg)
      map <- lapply(map, function(ids) unique(na.omit(as.numeric(ids))))
      map <- map[lengths(map) > 0]

      if (length(map) == 0) next

      gene_list <- unique(unlist(map))
      gene_list <- gene_list[!is.na(gene_list)]
      if (length(gene_list) == 0) next

      # 过滤重复 gene 列表（防止 vote bias）
      gene_list_key <- paste(sort(gene_list), collapse = "_")
      if (any(vapply(gene_lists, function(g) identical(sort(g), sort(gene_list)), logical(1)))) next

      if (verbose) {
        message(sprintf("  Accepted: %d genes, %d CpGs (stat > %.2f, dist < %d)",
                        length(gene_list), length(map), r, d))
      }

      gene_lists[[length(gene_lists) + 1]] <- gene_list
      cpg_gene_maps[[length(cpg_gene_maps) + 1]] <- map
      params[[length(params) + 1]] <- c(stat = r, d = d)
    }
  }

  return(list(
    gene_lists = gene_lists,
    cpg_gene_maps = cpg_gene_maps,
    params = params
  ))
}

#' @title Run Enrichment Analysis on Gene Lists
#'
#' @description Performs pathway enrichment analysis on a gene list using specified databases.
#' @param gene_list Character vector of Entrez IDs.
#' @param databases Character vector of databases (e.g., "Reactome", "GO", "KEGG").
#' @param r Numeric, statistics threshold used for filtering.
#' @param d Numeric, distance threshold used for filtering.
#' @param verbose Logical, whether to print progress messages.
#' @return A list of enrichment results for each database.
#' @export
run_enrichment <- function(gene_list, databases, r, d, verbose = FALSE) {
  if (verbose) message(sprintf("Running enrichment for statistics threshold = %.2f, d = %d", r, d))
  enrich_results <- list()
  if ("Reactome" %in% databases) {
    if (verbose) message("  Analyzing Reactome...")
    enrich_results$Reactome <- as.data.frame(
      enrichPathway(gene = gene_list, organism = "human", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", readable = TRUE)
    )
  }
  if ("GO" %in% databases) {
    if (verbose) message("  Analyzing GO...")
    enrich_results$GO <- as.data.frame(
      enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", readable = TRUE)
    )
  }
  if ("KEGG" %in% databases) {
    if (verbose) message("  Analyzing KEGG...")
    kegg_enrich <- enrichKEGG(gene = gene_list, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")
    kegg_enrich_readable <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    enrich_results$KEGG <- as.data.frame(kegg_enrich_readable)
  }
  if (verbose) message("  Enrichment completed.")
  return(enrich_results)
}

#' @title Combine Enrichment Results
#'
#' @description Combines enrichment results from multiple gene lists using harmonic mean p-value.
#' @param enrich_results List of enrichment results from run_enrichment.
#' @param databases Character vector of databases.
#' @param verbose Logical, whether to print progress messages.
#' @return A list of combined result tables for each database.
#' @export
combine_enrichment_results <- function(enrich_results, databases, verbose = FALSE) {
  if (verbose) message("Combining enrichment results...")

  all_pathways <- list()
  for (db in databases) {
    if (verbose) message(sprintf("  Extracting pathways for %s...", db))
    all_pathways[[db]] <- unique(unlist(lapply(enrich_results, function(x) x[[db]]$ID)))
  }

  p_value_matrices <- list()
  for (db in databases) {
    if (verbose) message(sprintf("  Building p-value matrix for %s...", db))
    p_value_matrices[[db]] <- matrix(NA, nrow = length(all_pathways[[db]]), ncol = length(enrich_results))
    rownames(p_value_matrices[[db]]) <- all_pathways[[db]]
    for (j in 1:length(enrich_results)) {
      current_df <- enrich_results[[j]][[db]]
      for (pathway in all_pathways[[db]]) {
        p_value_matrices[[db]][pathway, j] <- ifelse(pathway %in% current_df$ID, current_df$pvalue[current_df$ID == pathway], 1)
      }
    }
  }

  harmonic_mean_p <- function(p_values) {
    valid_p <- p_values[!is.na(p_values) & p_values >= 0 & p_values <= 1]
    if (length(valid_p) == 0) return(1)
    if (length(valid_p) == 1) return(valid_p[1])
    return(length(valid_p) / sum(1 / valid_p))
  }

  combined_p <- list()
  for (db in databases) {
    if (verbose) message(sprintf("  Combining p-values for %s using HMP...", db))
    combined_p[[db]] <- c()
    for (pathway in all_pathways[[db]]) {
      p_values <- p_value_matrices[[db]][pathway, ]
      combined_p[[db]][pathway] <- harmonic_mean_p(p_values)
    }
    combined_p[[db]] <- combined_p[[db]][order(combined_p[[db]])]
    combined_p[[db]] <- p.adjust(combined_p[[db]], method = "BH")
  }

  pathway_info <- list()
  for (db in databases) {
    if (verbose) message(sprintf("  Collecting pathway info for %s...", db))
    pathway_info[[db]] <- list()
    for (i in 1:length(enrich_results)) {
      if (nrow(enrich_results[[i]][[db]]) > 0) {
        df <- enrich_results[[i]][[db]][, c("ID", "Description", "geneID")]
        for (j in 1:nrow(df)) {
          pathway_id <- df$ID[j]
          if (is.null(pathway_info[[db]][[pathway_id]])) {
            pathway_info[[db]][[pathway_id]] <- list(Description = df$Description[j], geneIDs = character())
          }
          genes <- unlist(strsplit(df$geneID[j], "/"))
          pathway_info[[db]][[pathway_id]]$geneIDs <- union(pathway_info[[db]][[pathway_id]]$geneIDs, genes)
        }
      }
    }
  }

  result_tables <- list()
  for (db in databases) {
    if (verbose) message(sprintf("  Generating result table for %s...", db))
    pathway_df <- data.frame(
      ID = names(pathway_info[[db]]),
      Description = sapply(pathway_info[[db]], function(x) x$Description),
      geneID = sapply(pathway_info[[db]], function(x) paste(x$geneIDs, collapse = "/")),
      stringsAsFactors = FALSE
    )
    result_tables[[db]] <- data.frame(
      ID = names(combined_p[[db]]),
      p.adjust = combined_p[[db]],
      row.names = NULL
    )
    result_tables[[db]] <- merge(result_tables[[db]], pathway_df, by = "ID", all.x = TRUE)
    result_tables[[db]] <- result_tables[[db]][order(result_tables[[db]]$p.adjust), ]
  }

  if (verbose) message("Result combination completed.")
  return(result_tables)
}

#' @title Pathway Vote Algorithm for eQTM Data (Auto Parallel)
#'
#' @description Performs pathway enrichment analysis using a voting-based approach on eQTM data.
#' @param ewas_data A data.frame with columns: cpg and a ranking column (e.g., p_value, score).
#' @param eQTM An eQTM object containing eQTM data.
#' @param k_values A numeric vector of top k CpGs to select (e.g., c(10, 50, 100)).
#' @param stat_grid A numeric vector of statistics thresholds.
#' @param distance_grid A numeric vector of distance thresholds.
#' @param overlap_threshold A numeric value for gene list overlap threshold.
#' @param databases A character vector of pathway databases (e.g., "Reactome").
#' @param rank_column A character string indicating which column in `ewas_data` to use for ranking (default: "p_value").
#' @param rank_decreasing Logical. If TRUE (default), sorts CpGs from high to low based on `rank_column`.
#'                        If FALSE, sorts from low to high. Use FALSE for p-values; TRUE for absolute statistics.
#' @param use_abs Logical. Whether to apply `abs()` to the ranking column (e.g., p-value, correlation, score) before sorting CpGs.
#' @param min_vote_support Minimum number of enrichment combinations in which a pathway must appear to be retained. Default = 3.
#' @param min_genes Minimum number of genes (`Count`) a pathway must include in any enrichment result to be considered. Default = 3.
#' @param verbose Logical, whether to print progress messages.
#' @return A list containing enrichment results and CpG-gene mappings.
#' @export
pathway_vote <- function(ewas_data, eQTM, k_values, stat_grid, distance_grid,
                         overlap_threshold, databases = c("Reactome"),
                         rank_column = "p_value",
                         rank_decreasing = FALSE,
                         use_abs = FALSE,
                         min_vote_support = 3,
                         min_genes = 3,
                         verbose = FALSE) {

  # ---- Load and setup parallel environment ----
  required_pkgs <- c("PathwayVote", "furrr", "future", "ReactomePA", "clusterProfiler", "org.Hs.eg.db")
  lapply(required_pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required. Please install it first.")
    }
  })
  suppressMessages({
    lapply(required_pkgs, library, character.only = TRUE)
  })

  available_cores <- parallelly::availableCores(logical = TRUE)
  max_safe_workers <- floor(min(available_cores, parallelly::availableCores("system")) * 0.75)
  max_safe_workers <- max(1, max_safe_workers)

  if (!inherits(future::plan(), "multisession")) {
    future::plan(future::multisession, workers = max_safe_workers)
  }
  if (verbose) message("Using ", max_safe_workers, " parallel workers.")

  # ---- Input checks ----
  if (!inherits(eQTM, "eQTM")) stop("eQTM must be an eQTM object")
  if (!"cpg" %in% colnames(ewas_data)) stop("ewas_data must contain a 'cpg' column.")
  if (!rank_column %in% colnames(ewas_data)) stop("Column '", rank_column, "' not found in ewas_data")
  if (all(is.na(getData(eQTM)$entrez))) stop("Entrez IDs are required for pathway analysis")

  # ---- Sort EWAS data ----
  if (!rank_column %in% colnames(ewas_data)) stop("Column '", rank_column, "' not found in ewas_data")

  ranking_values <- if (use_abs) abs(ewas_data[[rank_column]]) else ewas_data[[rank_column]]
  ewas_data <- ewas_data[order(ranking_values, decreasing = rank_decreasing), ]

  # ---- Filter gene lists for each k ----
  if (verbose) message("Generating gene lists for all k values...")
  all_gene_sets <- list()
  valid_combination_count <- 0  # 统计有效组合数量

  for (k in k_values) {
    if (verbose) message(sprintf("Processing top %d CpGs...", k))
    selected_cpgs <- head(ewas_data$cpg, k)
    eQTM_subset <- new("eQTM", data = getData(eQTM)[getData(eQTM)$cpg %in% selected_cpgs, ],
                       metadata = getMetadata(eQTM))

    filtered_results <- filter_gene_lists(eQTM_subset, stat_grid, distance_grid, overlap_threshold, verbose = FALSE)

    for (i in seq_along(filtered_results$gene_lists)) {
      gene_list_i <- filtered_results$gene_lists[[i]]
      if (length(gene_list_i) == 0) next  # 跳过无效组合

      valid_combination_count <- valid_combination_count + 1

      all_gene_sets[[length(all_gene_sets) + 1]] <- list(
        gene_list = gene_list_i,
        cpg_gene_map = filtered_results$cpg_gene_maps[[i]],
        param = filtered_results$params[[i]],
        k = k
      )
    }
  }

  if (verbose) message(sprintf("Gene filtering completed. %d valid combinations retained.", valid_combination_count))

  # ---- Run enrichment analysis in parallel ----
  if (verbose) message("Running enrichment analysis (parallel)...")
  enrich_results <- furrr::future_map(
    all_gene_sets,
    function(x) {
      run_enrichment(
        gene_list = x$gene_list,
        databases = databases,
        r = x$param["stat"],
        d = x$param["d"],
        verbose = verbose
      )
    },
    .options = furrr::furrr_options(
      seed = TRUE,
      packages = c("PathwayVote", "ReactomePA", "clusterProfiler", "org.Hs.eg.db")
    )
  )

  # ---- Prune pathways by vote ----
  enrich_results <- prune_pathways_by_vote(
    enrich_results,
    min_vote_support = min_vote_support,
    min_genes = min_genes_per_hit,
    verbose = verbose
  )

  # ---- Combine and return results ----
  cpg_gene_maps <- lapply(all_gene_sets, function(x) x$cpg_gene_map)
  result_tables <- combine_enrichment_results(enrich_results, databases, verbose = verbose)
  return(list(result_tables = result_tables, cpg_gene_maps = cpg_gene_maps))
}

