logit <- function(x) {
  log(x) - log(1 - x)
}

expit <- function(x) {
  exp(x) / (1 + exp(x))
}

subset_features <- function(x, p) {
  if(p > ncol(x))
    stop("p is greater than the number of features in x!")

  features_select <- x %>%
    apply(2, mean) %>%
    {order(-.)} %>%
    {.[seq(1, p)]}
  x <- x[, features_select]

  ## renormalize
  if(any(apply(x == 0, 1, all)))
    stop("Subsetting unsuccessful!")
  x <- t(apply(x, 1, function(x) x / sum(x)))

  return(x)
}

add_pseudo <- function(x) {
  if(all(x > 0))
    return(x)

  x <- x + min(setdiff(x, 0)) / 2
  x <- t(apply(x, 1, function(x) x / sum(x)))

  return(x)
}

subset_pseudo <- function(x, subset = TRUE) {
  feature_subset <- !apply(x == 0, 2, all)
  x_old <- x

  x <- add_pseudo(x_old[, feature_subset])

  if(!subset) {
    x_old[, feature_subset] <- x
    x <- x_old
  }

  return(list(x = x,
              feature_subset = feature_subset))
}

tss_withzero <- function(x) {
  if(any(x < 0))
    stop("Data shouldn't be TSS'ed!")
  if(all(x == 0))
    return(x)
  return(x / sum(x))
}
