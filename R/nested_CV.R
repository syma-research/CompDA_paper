nested_cv_pred <- function(x, y, K = 10) {
  tb_nest <- tibble::tibble(
    k = seq(1, K)
  ) %>%
    dplyr::arrange(k) %>%
    dplyr::mutate(folds = caret::createFolds(y = y, k = K)) %>%
    dplyr::ungroup()

  predictions <- seq(1, nrow(tb_nest)) %>%
    future.apply::future_lapply(
      function(i_run) {
        i_fold <- tb_nest$folds[[i_run]]
        x_train <- x[-i_fold, , drop = FALSE]
        y_train <- y[-i_fold]

        fitted_train <-
          glmnet::cv.glmnet(x = x_train,
                            y = y_train,
                            family = "binomial",
                            alpha = 1,
                            nfolds = 10)
        yhat <- glmnet:::predict.cv.glmnet(
          fitted_train,
          newx = x[i_fold, , drop = FALSE],
          s = "lambda.min")[, 1]
        yhat <- (exp(yhat) / (1 + exp(yhat)))
        return(yhat)
      },
      future.seed = TRUE
    )

  n <- length(y)
  y_pred <- y
  for(i_run in seq(1, nrow(tb_nest))) {
    y_pred[tb_nest$folds[[i_run]]] <- predictions[[i_run]]
  }

  return(y_pred)
}
