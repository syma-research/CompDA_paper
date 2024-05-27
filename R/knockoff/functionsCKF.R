#### screeningStep - runs the screening with MIO
#### Inputs:
#### X - feature matrix
#### Y - response vector
#### nScreen - screening set size
#### n0 - size of screening sample
#### n1 - size of hold out sample
#### timelim - maximum runtime for MIO

screeningStep = function(X0,Y0,nScreen, n0, n1, timelim, dir_warning){
  ##### Load in MIO functions (change path)
  
  library(knockoff)
  library(glmnet)
  library(MethylCapSig)
  library(GUniFrac)
  library(energy)
  library(dirmult)
  library(gurobi)
  
  n = nrow(X0)
  ##### Run MIO to screen
  MIOFit <- suppressWarnings(
    withCallingHandlers(
      bsConst(X0,Y0, form = 4, k = nScreen, 
              intercept = F, time.limit = timelim),
      warning = function(w) {
        write(w$message, dir_warning)
      }))
  # MIOFit = suppressWarnings(bsConst(X0,Y0, form = 4, k = nScreen, intercept = F, time.limit = timelim))
  selected = which(MIOFit$beta != 0)
  
  #### return response
  res = list()
  res[[1]] = selected
  res[[2]] = X0
  res[[3]] = Y0
  res[[4]] = MIOFit
  return(res)
}


sumToZeroLm = function(X,Y){
  library(limSolve)
  A = X
  B = Y
  E = rep(1, ncol(X))
  FRes = 0
  fit = lsei(A,B,E=E, F=FRes)
  return(fit)
}


#### selectionStep - runs the selection procedure with lambda max stat
#### X,Y,X0,Y0,X1,Y1 - datasets submitted
#### screened - screening set
#### type - kn+ = 1, kn = 0
#### FDR - target FDR
selectionStep = function(Y0,Y1,W0,W1,screened, FDR, 
                         dir_test){
  library(knockoff)
  
  # Renormalize composotion
  W0Sel = W0[,screened]
  W1Sel = W1[,screened]
  
  W0tot = rowSums(W0Sel)
  W1tot = rowSums(W1Sel)
  Z0n = W0Sel
  Z1n = W1Sel
  
  for(i in 1:nrow(W0Sel)){
    for(j in 1:ncol(W0Sel)){
      Z0n[i,j] = log(W0Sel[i,j]/W0tot[i])
    }
  }
  
  for(i in 1:nrow(W1Sel)){
    for(j in 1:ncol(W1Sel)){
      Z1n[i,j] = log(W1Sel[i,j]/W1tot[i])
    }
  }
  X0New = scale(Z0n, center = T, scale = F)
  X1New = scale(Z1n, center = T, scale = F)
  
  ### create knockoff matrix with recycling
  X0Sel = normc(X0New, center = F)
  X1Sel = normc(X1New, center = F)
  
  write("this is a test", paste0(dir_test, "1"))
  # create X1 tilde 
  X1SelK = knockoff::create.fixed(X1Sel)$Xk
  write("this is a test", paste0(dir_test, "2"))
  
  # the knockoff matrix is X0 on top of X1 tilde
  Xk = rbind(X0Sel, X1SelK)
  
  ### run lasso
  XSel = rbind(X0Sel,X1Sel)
  responseVector = c(Y0,Y1)
  W = stat.glmnet_lambdasmax(XSel,Xk,responseVector)
  
  ### select using KF
  t = knockoff.threshold(W, FDR, 0)
  selected = which(W>=t)
  
  tp = knockoff.threshold(W, FDR, 1)
  selectedp = which(W>=tp)
  
  ### return response
  res = list()
  res[[1]] = selectedp
  res[[2]] = selected
  res[[3]] = c(t,tp)
  
  return(res)
}

########### Utility functions
normc = function(X,center=T) {
  X.centered = scale(X, center=center, scale=F)
  X.scaled = scale(X.centered, center=F, scale=sqrt(colSums(X.centered^2)))
  X.scaled[,] # No attributes
}
listToMat = function(l){
  nL = length(l)
  nC = length(l[[1]])
  m = matrix(0, nrow = nL, ncol = nC)
  for(i in 1:nL){
    m[i,] = as.numeric(l[[i]])
  }
  
  return(m)
}
###########

########### lambda signed max statistic
lambdaFirstIn = function(coef, lambda){
  coef = t(coef)
  nL = length(lambda)
  p = ncol(coef)/2
  
  vals = apply(coef, 1,nonZero)
  firstIn = lambda[vals]
  return(firstIn)
}
nonZero = function(x){
  w = which(x!=0)
  if(length(w) > 0){
    l = min(w)
  }
  else{
    l = length(x)
  }
  return(l)
}
lambdaSignedMax = function(Z){
  p = length(Z)/2
  orig = Z[1:p]
  kn = Z[(p+1):(2*p)]
  W = pmax(orig, kn)*sign(orig-kn)
}
###########

### LS diff stat 
ZSLSDiff = function(beta){
  p = length(beta)/2
  orig = 1:p
  kn = (p+1):(2*p)
  
  W = abs(beta[orig]) - abs(beta[kn])
  return(W)
}

sanityCheck = function(X0,Y0,X1,Y1, X, Y, relativeZ, relativeZ0, relativeZ1){
  
  Z0dupe = sum(duplicated(t(relativeZ0)))
  Z1dupe = sum(duplicated(t(relativeZ1)))
  
  Zavg = mean(rowSums(exp(relativeZ)))
  Z0avg = mean(rowSums(exp(relativeZ0)))
  Z1avg = mean(rowSums(exp(relativeZ1)))
  
  sanity = F
  Zrows = (Zavg+Z0avg+Z1avg)/3
  dupeCheck = Z0dupe + Z1dupe
  rowCheck = Zrows == 1
  dupeCheckSet = dupeCheck == 0 
  if(rowCheck & dupeCheckSet == T){
    sanity = T
  }
  return(sanity)
}