#' @importFrom stats na.omit p.adjust
#' @importFrom utils head
#' @importFrom future plan
#' @importFrom furrr future_map furrr_options
#' @importFrom parallelly availableCores
#' @importFrom ReactomePA enrichPathway
#' @importFrom clusterProfiler enrichGO enrichKEGG setReadable
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @import dplyr

#' @title Pathway Vote
#'
#' @description Performs methylation pathway enrichment analysis using a voting-based approach.
#'
#' @param ewas_data A data.frame with columns: cpg and a ranking column (e.g., p_value, score).
#' @param eQTM An eQTM object containing eQTM data.
#' @param databases A character vector of pathway databases (e.g., "Reactome").
#' @param rank_column A character string indicating which column in `ewas_data` to use for ranking.
#' @param rank_decreasing Logical. If TRUE (default), sorts CpGs from high to low based on `rank_column`.
#' @param use_abs Logical. Whether to apply `abs()` to the ranking column before sorting CpGs.
#' @param k_grid A numeric vector of top k CpGs to select. If NULL, inferred automatically.
#' @param stat_grid A numeric vector of statistics thresholds. If NULL, inferred automatically.
#' @param distance_grid A numeric vector of distance thresholds. If NULL, inferred automatically.
#' @param overlap_threshold A numeric value for gene list overlap threshold. If NULL, inferred automatically.
#' @param fixed_prune Integer or NULL. Minimum number of votes to retain a pathway. If NULL, will use cuberoot(N) where N is the number of enrichment runs.
#' @param grid_size Integer. Number of values in each grid when auto-generating. Default is 5.
#' @param min_genes_per_hit Minimum number of genes (`Count`) a pathway must include to be considered.
#' @param workers Optional integer. Number of parallel workers. If NULL, use 2 logical cores.
#' @param readable Logical. whether to convert Entrez IDs to gene symbols in enrichment results.
#' @param verbose Logical. whether to print progress messages.
#'
#' @return A named list of data.frames, each containing enrichment results (pathway ID, p.adjust, Description, geneID) for one database (e.g., Reactome, KEGG).
#'
#' @examples
#' set.seed(123)
#'
#' # Simulated EWAS result: a mix of signal and noise
#' ewas_data <- data.frame(
#'   cpg = paste0("cg", sprintf("%06d", 1:20)),
#'   p_value = c(runif(5, 1e-5, 0.001), runif(5, 0.01, 0.05), runif(10, 0.1, 1))
#' )
#'
#' # Corresponding eQTM mapping (some of these CpGs have gene links)
#' eqtm_data <- data.frame(
#'   cpg = ewas_data$cpg,
#'   statistics = rnorm(20, mean = 2, sd = 1),
#'   p_value = runif(20, 0.001, 0.05),
#'   distance = sample(1000:100000, 20),
#'   entrez = rep(c("5290", "673", "1956", "7157", "7422"), length.out = 20)
#' )
#' eqtm_obj <- create_eQTM(eqtm_data)
#'
#' results <- pathway_vote(
#'   ewas_data = data,
#'   eQTM = eqtm_obj,
#'   databases = c("KEGG"),
#'   rank_column = "p_value",
#'   rank_decreasing = FALSE,
#'   use_abs = FALSE,
#'   worker = 1, # If not specified, will use 2 cores by default
#'   verbose = TRUE
#' )
#'
#' @export
#'
pathway_vote <- function(ewas_data, eQTM,
                         databases = c("Reactome"),
                         rank_column = "p_value",
                         rank_decreasing = FALSE,
                         use_abs = FALSE,
                         k_grid = NULL,
                         stat_grid = NULL,
                         distance_grid = NULL,
                         overlap_threshold = NULL,
                         fixed_prune = NULL,
                         grid_size = 5,
                         min_genes_per_hit = 3,
                         readable = FALSE,
                         workers = NULL,
                         verbose = FALSE) {

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
  user_specified_workers <- !is.null(workers)

  if (is.null(workers)) {
    workers <- min(2, available_cores)
  } else {
    if (!is.numeric(workers) || length(workers) != 1 || workers < 1) {
      stop("`workers` must be a positive integer")
    }
    workers <- min(workers, available_cores)
  }

  current_plan <- future::plan()
  if (!inherits(current_plan, "multisession") || future::nbrOfWorkers() != workers) {
    future::plan(future::multisession, workers = workers)
  }

  if (verbose) {
    message("Using ", workers, ifelse(workers == 1, " worker", " workers"))
  }

  if (!inherits(eQTM, "eQTM")) stop("eQTM must be an eQTM object")
  if (!"cpg" %in% colnames(ewas_data)) stop("ewas_data must contain a 'cpg' column.")
  if (!rank_column %in% colnames(ewas_data)) stop("Column '", rank_column, "' not found in ewas_data")
  if (all(is.na(getData(eQTM)$entrez))) stop("Entrez IDs are required for pathway analysis")

  if (is.null(k_grid)) {
    k_grid <- auto_generate_k_grid_inflection(ewas_data, rank_column, rank_decreasing, grid_size, verbose)
  }

  if (is.null(stat_grid)) {
    stat_vals <- abs(getData(eQTM)$statistics)
    stat_grid <- round(quantile(stat_vals, probs = seq(0.2, 0.8, length.out = grid_size), na.rm = TRUE), 2)
  }

  if (is.null(distance_grid)) {
    dist_vals <- getData(eQTM)$distance
    distance_grid <- round(quantile(dist_vals, probs = seq(0.5, 0.9, length.out = grid_size), na.rm = TRUE), -3)
  }

  ranking_values <- if (use_abs) abs(ewas_data[[rank_column]]) else ewas_data[[rank_column]]
  ewas_data <- ewas_data[order(ranking_values, decreasing = rank_decreasing), ]

  if (verbose) message("Generating gene lists for all k values...")
  all_gene_sets <- list()
  valid_combination_count <- 0

  for (k in k_grid) {
    if (verbose) message(sprintf("Processing top %d CpGs...", k))
    selected_cpgs <- head(ewas_data$cpg, k)
    eQTM_subset <- new("eQTM", data = getData(eQTM)[getData(eQTM)$cpg %in% selected_cpgs, ],
                       metadata = getMetadata(eQTM))

    if (is.null(overlap_threshold)) {
      temp_lists <- filter_gene_lists(eQTM_subset, stat_grid, distance_grid, overlap_threshold = 1)
      overlap_threshold <- auto_overlap_threshold(temp_lists$gene_lists, quantile_level = 1 - 1/grid_size, verbose)
    }

    filtered_results <- filter_gene_lists(eQTM_subset, stat_grid, distance_grid, overlap_threshold, verbose = FALSE)

    for (i in seq_along(filtered_results$gene_lists)) {
      gene_list_i <- filtered_results$gene_lists[[i]]
      if (length(gene_list_i) == 0) next

      valid_combination_count <- valid_combination_count + 1

      all_gene_sets[[length(all_gene_sets) + 1]] <- list(
        gene_list = gene_list_i,
        param = filtered_results$params[[i]],
        k = k
      )
    }
  }

  if (verbose) message(sprintf("Gene filtering completed. %d valid combinations retained.", valid_combination_count))

  if (verbose) message("Running enrichment analysis...")

  enrich_results <- furrr::future_map(
    all_gene_sets,
    function(x) {
      tryCatch({
        run_enrichment(
          gene_list = x$gene_list,
          databases = databases,
          readable = readable,
          verbose = FALSE
        )
      }, error = function(e) {
        warning("Enrichment failed for one gene list: ", conditionMessage(e))
        return(NULL)
      })
    },
    .options = furrr::furrr_options(
      seed = TRUE,
      packages = c("PathwayVote", "ReactomePA", "clusterProfiler", "org.Hs.eg.db")
    )
  )

  enrich_results <- prune_pathways_by_vote(
    enrich_results,
    fixed_prune = fixed_prune,
    min_genes = min_genes_per_hit,
    verbose = verbose
  )

  result_tables <- combine_enrichment_results(enrich_results, databases, verbose = verbose)
  return(result_tables)
}


