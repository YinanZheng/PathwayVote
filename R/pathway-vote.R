#' @importFrom stats na.omit p.adjust
#' @importFrom utils head
#' @importFrom future plan
#' @importFrom furrr future_map furrr_options
#' @importFrom parallelly availableCores
#' @importFrom ReactomePA enrichPathway
#' @importFrom clusterProfiler enrichGO enrichKEGG setReadable
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @import dplyr

#' @title Run Enrichment Analysis on Gene Lists
#'
#' @description Performs pathway enrichment analysis on a gene list using specified databases.
#' @param gene_list Character vector of Entrez IDs.
#' @param databases Character vector of databases (e.g., "Reactome", "GO", "KEGG").
#' @param r Numeric, statistics threshold used for filtering.
#' @param d Numeric, distance threshold used for filtering.
#' @param readable Logical, whether to convert Entrez IDs to gene symbols in enrichment results. Default is FALSE to retain Entrez IDs (recommended for programmatic comparison).
#' @param verbose Logical, whether to print progress messages.
#' @return A list of enrichment results for each database.
#' @export
run_enrichment <- function(gene_list, databases, r, d, readable = FALSE, verbose = FALSE) {
  if (verbose) message(sprintf("Running enrichment for statistics threshold = %.2f, d = %d", r, d))
  enrich_results <- list()
  if ("Reactome" %in% databases) {
    if (verbose) message("  Analyzing Reactome...")
    enrich_results$Reactome <- as.data.frame(
      enrichPathway(gene = gene_list, organism = "human", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", readable = readable)
    )
  }
  if ("GO" %in% databases) {
    if (verbose) message("  Analyzing GO...")
    enrich_results$GO <- as.data.frame(
      enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", readable = readable)
    )
  }
  if ("KEGG" %in% databases) {
    if (verbose) message("  Analyzing KEGG...")
    kegg_enrich <- enrichKEGG(gene = gene_list, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")
    if (readable) {
      kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    }
    enrich_results$KEGG <- as.data.frame(kegg_enrich)
  }
  if (verbose) message("  Enrichment completed.")
  return(enrich_results)
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
#' @param min_genes_per_hit Minimum number of genes (`Count`) a pathway must include in any enrichment result to be considered. Default = 3.
#' @param readable Logical, whether to convert Entrez IDs to gene symbols in enrichment results. Default is FALSE to retain Entrez IDs (recommended for programmatic comparison).
#' @param verbose Logical, whether to print progress messages.
#' @return A list of enrichment results.
#' @export
pathway_vote <- function(ewas_data, eQTM, k_values, stat_grid, distance_grid,
                         overlap_threshold, databases = c("Reactome"),
                         rank_column = "p_value",
                         rank_decreasing = FALSE,
                         use_abs = FALSE,
                         min_vote_support = 3,
                         min_genes_per_hit = 3,
                         readable = FALSE,
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
      if (length(gene_list_i) == 0) next  # skip invalid combo

      valid_combination_count <- valid_combination_count + 1

      all_gene_sets[[length(all_gene_sets) + 1]] <- list(
        gene_list = gene_list_i,
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
        readable = readable,
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
  result_tables <- combine_enrichment_results(enrich_results, databases, verbose = verbose)
  return(result_tables)
}

