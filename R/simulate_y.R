simulate_y <- function(x, effect_size, n_signal, epsilon,
                       family = "gaussian",
                       seed = 0) {
  if(n_signal != round(n_signal))
    stop("n_signal must be an integer!")

  set.seed(seed)
  ind_TP <- seq(1, ncol(x)) %in%
    sample(seq(1, ncol(x)), size = n_signal, replace = FALSE)
  effects <- runif(min = -2, max = 2, n = n_signal) * effect_size

  xTP_renorm <- x[, ind_TP, drop = FALSE] %>%
    apply(1, tss_withzero) %>%
    t()

  x_TP <- log(xTP_renorm + epsilon)
  x_TP <-
    apply(x_TP, 2, function(x) x - mean(x))

  link_ey <- (x_TP %*% effects)[, 1]
  if(family == "gaussian") {
    y <- link_ey + rnorm(n = nrow(x), mean = 0, sd = 1)
    snr <- var(link_ey)
  }
  else if (family == "binomial") {
    ey <- exp(link_ey) / (1 + exp(link_ey))
    y <- rbinom(n = nrow(x), size = 1, prob = ey)
    snr <- var(ey) / (var(y) - var(ey))
  }

  return(list(ind_TP = ind_TP,
              effects = effects,
              x_TP = x_TP,
              link_ey = link_ey,
              y = y,
              snr = snr,
              seed = seed))
}

# the association between y and x is squred
simulate_y_sqr <- function(x, effect_size, n_signal, epsilon,
                           family = "gaussian",
                           seed = 0) {
  if(n_signal != round(n_signal))
    stop("n_signal must be an integer!")

  set.seed(seed)
  ind_TP <- seq(1, ncol(x)) %in%
    sample(seq(1, ncol(x)), size = n_signal, replace = FALSE)
  effects <- runif(min = -2, max = 2, n = n_signal) * effect_size

  xTP_renorm <- x[, ind_TP, drop = FALSE] %>%
    apply(1, tss_withzero) %>%
    t()
  x_TP <- log(xTP_renorm + epsilon)
  x_TP <-
    apply(x_TP, 2, function(x) x - mean(x))
  x_TP_sqr <-
    apply(x_TP^2, 2, function(x) x - mean(x))

  link_ey <- (x_TP_sqr %*% effects)[, 1]
  if(family == "gaussian") {
    link_ey <- link_ey * sqrt(snr) / sd(linke_ey)
    y <- link_ey + rnorm(n = nrow(x), mean = 0, sd = 1)
  }
  else if (family == "binomial") {
    y <- rbinom(n = nrow(x), size = 1, prob = exp(link_ey) / (1 + exp(link_ey)))
  }
  snr <- var(link_ey) / 1

  return(list(ind_TP = ind_TP,
              effects = effects,
              x_TP = x_TP,
              link_ey = link_ey,
              y = y,
              snr = snr,
              seed = seed))
}

# when there are additional confounder effects
simulate_y_confounded <-
  function(x, z, effect_size, effect_confounded,
           n_signal, epsilon,
           family = "gaussian",
           seed = 0) {
    if(n_signal != round(n_signal))
      stop("n_signal must be an integer!")

    set.seed(seed)
    ind_TP <- seq(1, ncol(x)) %in%
      sample(seq(1, ncol(x)), size = n_signal, replace = FALSE)
    effects <- runif(min = -2, max = 2, n = n_signal) * effect_size
    effects_confounded <- runif(min = -2, max = 2, n = ncol(z))

    xTP_renorm <- x[, ind_TP, drop = FALSE] %>%
      apply(1, tss_withzero) %>%
      t()
    x_TP <- log(xTP_renorm + epsilon)
    x_TP <-apply(x_TP, 2, function(x) x - mean(x))
    z <- apply(z, 2, function(x) x - mean(x))

    link_ey_x <- (x_TP %*% effects)[, 1]
    link_ey_z <- (z %*% effects_confounded)[, 1]
    link_ey_z <- link_ey_z / sd(link_ey_z) * sd(link_ey_x) * effect_confounded
    link_ey <- link_ey_x + link_ey_z
    if(family == "gaussian") {
      link_ey <- link_ey * sqrt(snr) / sd(linke_ey)
      y <- link_ey + rnorm(n = nrow(x), mean = 0, sd = 1)
    }
    else if (family == "binomial") {
      y <- rbinom(n = nrow(x), size = 1, prob = exp(link_ey) / (1 + exp(link_ey)))
    }
    snr <- var(link_ey) / 1

    return(list(ind_TP = ind_TP,
                effects = effects,
                effects_confounded = effects_confounded,
                x_TP = x_TP,
                link_ey = link_ey,
                y = y,
                snr = snr,
                seed = seed))
  }
