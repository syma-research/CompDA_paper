# requires the MaAsLin2 package
test_DA <-
  function(x, y, covariates = NULL,
           epsilon,
           debug_path = NULL) {

    x <- t(x)
    rownames(x) <- paste0("Feature", seq(1, nrow(x)))

    metadata <- data.frame(y = y)
    if(!is.null(covariates)) {
      covariates <- as.data.frame(covariates)
      names_covariates <- paste0("metadata", seq(1, ncol(covariates)))
      colnames(covariates) <- names_covariates

      metadata <- cbind(metadata, covariates)
    }
    colnames(x) <-
      rownames(metadata) <-
      paste0("Sample", seq(1, ncol(x)))

    msg <- capture.output(
      fit_Maaslin2 <- Maaslin2::Maaslin2(
        input_data = x,
        input_metadata = metadata,
        output = debug_path,
        min_abundance = 0,
        min_prevalence = 0,
        normalization = "TSS",
        transform = "LOG",
        standardize = FALSE,
        max_significance = 1,
        plot_scatter = FALSE,
        plot_heatmap = FALSE
      ))


    unlink(debug_path, recursive = TRUE)

    results <- fit_Maaslin2$results %>%
      dplyr::filter(metadata == "y")
    p_return <- rep(NA_real_, nrow(x))
    names(p_return) <- rownames(x)
    p_return[results$feature] <- results$pval
    return(p_return)
  }

# requires the MicrobiomeStat package
test_linda <-
  function(x, y) {

    x <- t(x)
    rownames(x) <- paste0("Feature", seq(1, nrow(x)))
    ind_allzero <- !apply(x == 0, 2, all)

    msg <- capture.output(
      suppressWarnings(
        fit_linda <-
          MicrobiomeStat::linda(
            feature.dat = x[!apply(x == 0, 1, all), ind_allzero, drop = FALSE],
            meta.dat = data.frame(y = y[ind_allzero]),
            formula = "~y",
            feature.dat.type = "proportion",
            is.winsor = FALSE, # set this otherwise would remove non-zero values
            verbose = FALSE)$output$y)
    )

    p_return <- rep(NA_real_, nrow(x))
    names(p_return) <- rownames(x)
    p_return[rownames(fit_linda)] <- fit_linda$pvalue

    return(p_return)
  }


test_lasso <-
  function(x, y, covariates = NULL,
           family = "gaussian",
           epsilon,
           debug_path = NULL) {

    x_full <-
      x <-
      log(x + epsilon)

    if(!is.null(covariates))
      x_full <- as.matrix(cbind(x_full, covariates))

    fit_y <-
      glmnet::cv.glmnet(x = x_full,
                        y = y,
                        family = family,
                        alpha = 1)
    beta <- glmnet:::coef.cv.glmnet(
      fit_y,
      s = "lambda.min")[seq(2, ncol(x) + 1), 1]

    # non-zero betas are considered "significant" p-values
    return((beta == 0) * 1)
  }
