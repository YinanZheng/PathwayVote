#' @importFrom stats na.omit p.adjust
#' @importFrom utils head
#' @importFrom future plan
#' @importFrom furrr future_map furrr_options
#' @importFrom parallelly availableCores
#' @importFrom ReactomePA enrichPathway
#' @importFrom clusterProfiler enrichGO enrichKEGG setReadable
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @import dplyr

#' @title Pathway Voting-Based Enrichment Analysis
#'
#' @description Performs pathway enrichment analysis on CpGs using a voting-based approach.
#'
#' @param ewas_data A data.frame with columns: `cpg` and a numeric ranking column (e.g., p-value, t-statistics, variable importance). The second column will be used as the ranking metric.
#' @param eQTM An eQTM object containing eQTM data.
#' @param databases A character vector of pathway databases. Supporting: "Reactome", "KEGG", and "GO".
#' @param k_grid A numeric vector of top k CpGs to select. If NULL, inferred automatically.
#' @param stat_grid A numeric vector of statistics thresholds. If NULL, inferred automatically.
#' @param distance_grid A numeric vector of distance thresholds. If NULL, inferred automatically.
#' @param fixed_prune Integer or NULL. Minimum number of votes to retain a pathway. If NULL, will use cuberoot(N) where N is the number of enrichment runs.
#' @param grid_size Integer. Number of values in each grid when auto-generating. Default is 5.
#' @param min_genes_per_hit Minimum number of genes (`Count`) a pathway must include to be considered.
#' @param workers Optional integer. Number of parallel workers. If NULL, use 2 logical cores.
#' @param readable Logical. whether to convert Entrez IDs to gene symbols in enrichment results.
#' @param verbose Logical. whether to print progress messages.
#'
#' @return A named list of data.frames, each containing enrichment results for one database.
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
#'   worker = NULL, # If not specified, will use 2 cores by default
#'   verbose = TRUE
#' )
#'
#' @export
#'
pathway_vote <- function(ewas_data, eQTM,
                         databases = c("Reactome"),
                         k_grid = NULL,
                         stat_grid = NULL,
                         distance_grid = NULL,
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
  if (all(is.na(getData(eQTM)$entrez))) stop("Entrez IDs are required for pathway analysis")
  if (ncol(ewas_data) < 2) stop("ewas_data must have at least two columns: 'cpg' and a ranking column")

  # Automatically use the second column as ranking
  rank_column <- colnames(ewas_data)[2]
  rank_values <- ewas_data[[rank_column]]

  if (!is.numeric(rank_values)) {
    stop("The second column of ewas_data must be numeric for ranking purposes.")
  }

  is_pval_like <- all(rank_values >= 0 & rank_values <= 1, na.rm = TRUE)
  use_abs <- !is_pval_like
  rank_decreasing <- !is_pval_like

  if (verbose) {
    message(sprintf("Auto-selected ranking: %s (decreasing = %s, absolute = %s)",
                    rank_column, rank_decreasing, use_abs))
  }

  if (is.null(k_grid)) {
    k_grid <- auto_generate_k_grid_inflection(ewas_data, rank_column, rank_decreasing, grid_size, verbose)
  }

  if (is.null(stat_grid)) {
    stat_vals <- abs(getData(eQTM)$statistics)
    stat_grid <- round(quantile(stat_vals, probs = seq(0.1, 0.9, length.out = grid_size), na.rm = TRUE), 2)
    if (verbose) message("Auto-selected stat_grid: ", paste(stat_grid, collapse = ", "))
  }

  if (is.null(distance_grid)) {
    dist_vals <- getData(eQTM)$distance
    distance_grid <- round(quantile(dist_vals, probs = seq(0.1, 1, length.out = grid_size), na.rm = TRUE), -3)
    if (verbose) message("Auto-selected distance_grid: ", paste(distance_grid, collapse = ", "))
  }

  ranking_values <- if (use_abs) abs(rank_values) else rank_values
  ewas_data <- ewas_data[order(ranking_values, decreasing = rank_decreasing), ]

  if (verbose) message("Generating gene lists for all k values...")
  all_gene_sets <- list()
  valid_combination_count <- 0

  for (k in k_grid) {
    if (verbose) message(sprintf("Processing top %d CpGs...", k))
    selected_cpgs <- head(ewas_data$cpg, k)
    eQTM_subset <- new("eQTM", data = getData(eQTM)[getData(eQTM)$cpg %in% selected_cpgs, ],
                       metadata = getMetadata(eQTM))

    # Step 1: Generate all candidate gene lists
    raw_results <- generate_gene_lists_grid(eQTM_subset, stat_grid, distance_grid, verbose = verbose)

    # Step 2: Apply entropy + stability-based pruning
    entropy_filtered_lists <- select_gene_lists_entropy_auto(
      gene_lists = raw_result$gene_lists,
      verbose = verbose
    )

    # Step 3: Match back to their corresponding stat/distance params
    kept_indices <- which(vapply(raw_result$gene_lists, function(x) any(sapply(entropy_filtered_lists, function(y) identical(x, y))), logical(1)))

    for (i in kept_indices) {
      gene_list_i <- raw_result$gene_lists[[i]]
      all_gene_sets[[length(all_gene_sets) + 1]] <- list(
        gene_list = gene_list_i,
        param = raw_result$params[[i]],
        k = k
      )
    }
  }

  if (verbose) message(sprintf("Gene filtering completed. %d valid combinations retained.", valid_combination_count))
  if (verbose) message("Running enrichment analysis...")

  gene_lists <- lapply(all_gene_sets, function(x) x$gene_list)

  enrich_results <- furrr::future_map(
    gene_lists,
    function(glist) {
      tryCatch({
        run_enrichment(
          gene_list = glist,
          databases = databases,
          readable = readable,
          verbose = FALSE
        )
      }, error = function(e) {
        warning("Enrichment failed: ", conditionMessage(e))
        NULL
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
