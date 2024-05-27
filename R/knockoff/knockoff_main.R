knockoff_R <- function(x, y, covariates = NULL, epsilon, timelimit = 60 * 10,
                       dir_debug) {
  n <- nrow(x)
  n0 <- round(nrow(x) * 0.4)
  p <- ncol(x)
  n1 <- n - n0
  nScreen <- round(n0 / log(n0)) * 2
  FDR <- 0.05

  sampleCol <- sample(1:n, size = n0, replace = F)
  ind_nonzero <- !apply(x[sampleCol, ] == 0, 2, all)
  relativeZ <- log(x + epsilon)[, ind_nonzero, drop = FALSE]
  relativeZ0 = relativeZ[sampleCol,]
  Y0 = y[sampleCol]
  W0 <- (x + epsilon)[sampleCol, ind_nonzero, drop = FALSE]
  W1 <- (x + epsilon)[-sampleCol, ind_nonzero, drop = FALSE]
  Y1 = y[-sampleCol]

  mioFit = screeningStep(relativeZ0, Y0, nScreen, n0 = n0, n1 = n1, timelimit,
                         dir_warning = paste0(dir_debug, "_screen_warning"))
  mioSel = mioFit[[1]]
  save(mioSel, file = paste0(dir_debug, "_selection.RData"))

  mSelect = selectionStep(Y0,Y1,W0,W1,mioSel,FDR,
                          dir_test = paste0(dir_debug, "_select_test_"))

  q1 <- q2 <- q3 <- rep(1, p)
  q1[ind_nonzero][mSelect[[1]]] <- 0
  q2[ind_nonzero][mSelect[[2]]] <- 0
  q3[ind_nonzero][mioFit[[1]]] <- 0

  return(cbind(q1, q2, q3))
}
