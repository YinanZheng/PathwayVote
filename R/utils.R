gene_filter <- function(eQTM, stat_threshold, distance) {
  if (!inherits(eQTM, "eQTM")) {
    stop("eQTM must be an eQTM object")
  }
  eQTM_data <- getData(eQTM)

  eQTM_filtered <- eQTM_data[abs(eQTM_data$statistics) >= stat_threshold & eQTM_data$distance <= distance, ]

  return(list(entrez = unique(na.omit(eQTM_filtered$entrez)), p_values = eQTM_filtered$p_value))
}

filter_gene_lists <- function(eQTM, stat_grid, distance_grid, overlap_threshold = 0.5, verbose = FALSE) {
  if (!inherits(eQTM, "eQTM")) stop("Input must be an eQTM object.")
  data <- getData(eQTM)

  gene_lists <- list()
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

      gene_list_key <- paste(sort(gene_list), collapse = "_")
      if (any(vapply(gene_lists, function(g) identical(sort(g), sort(gene_list)), logical(1)))) next

      if (verbose) {
        message(sprintf("  Accepted: %d genes, %d CpGs (stat > %.2f, dist < %d)",
                        length(gene_list), length(map), r, d))
      }

      gene_lists[[length(gene_lists) + 1]] <- gene_list
      params[[length(params) + 1]] <- c(stat = r, d = d)
    }
  }

  return(list(
    gene_lists = gene_lists,
    params = params
  ))
}

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

prune_pathways_by_vote <- function(enrich_results,
                                   min_vote_support = 2,
                                   min_genes = 2,
                                   verbose = TRUE) {
  all_hits <- unlist(lapply(enrich_results, function(x) {
    unlist(lapply(x, function(df) {
      if (is.null(df) || !"ID" %in% colnames(df) || !"Count" %in% colnames(df)) return(character(0))
      df$ID[df$Count >= min_genes]
    }))
  }))

  total_before <- length(unique(all_hits))
  freq_table <- table(all_hits)
  keep_ids <- names(freq_table)[freq_table >= min_vote_support]
  total_after <- length(keep_ids)

  if (verbose) {
    message(sprintf("Pruned pathways: %d retained out of %d total (min votes = %d, min genes = %d)",
                    total_after, total_before, min_vote_support, min_genes))
  }

  pruned_results <- lapply(enrich_results, function(x) {
    lapply(x, function(df) {
      if (is.null(df) || !"ID" %in% colnames(df) || !"Count" %in% colnames(df)) return(df)
      df[df$ID %in% keep_ids & df$Count >= min_genes, , drop = FALSE]
    })
  })

  return(pruned_results)
}
