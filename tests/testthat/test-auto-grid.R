test_that("auto_generate_k_grid_inflection returns reasonable k_grid", {
  ewas_data <- data.frame(
    cpg = paste0("cg", sprintf("%06d", 1:1000)),
    p_value = runif(1000, 1e-10, 1),
    statistics = rnorm(1000)
  )

  k_grid <- auto_generate_k_grid_inflection(
    ewas_data = ewas_data,
    rank_column = "p_value",
    rank_decreasing = FALSE,
    grid_size = 5,
    verbose = FALSE
  )

  expect_type(k_grid, "double")
  expect_true(length(k_grid) <= 5)
  expect_true(min(k_grid) >= 10)
  expect_true(max(k_grid) <= nrow(ewas_data))
})

test_that("auto_overlap_threshold computes quantile-based threshold", {
  set.seed(1)
  gene_lists <- list(
    sample(letters, 10),
    sample(letters, 12),
    sample(letters, 8),
    sample(letters, 9)
  )

  thres <- auto_overlap_threshold(gene_lists, quantile_level = 0.75, verbose = FALSE)

  expect_type(thres, "double")
  expect_true(thres >= 0 && thres <= 1)
})
