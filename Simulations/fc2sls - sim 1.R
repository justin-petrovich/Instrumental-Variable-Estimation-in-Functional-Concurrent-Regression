library(fcr)
library(MASS)
library(tidyverse)
library(refund)
library(itsadug)
library(fdapace)
library(face)

setwd("C:/Users/Owner/OneDrive - Saint Vincent College/Research/Time-Varying Elasticities")
setwd("./FCR-IV/Simulation Results")
source("fc2sls.R")
source("fc2sls-methods.R")
source("sampling scheme.R")
source("mse.R")
# source("FCReg_cv.R")
# source("llsmooth_mean.R")
# source("llsmooth_cov.R")


#-------------------------#
## Simulation Parameters ##
#-------------------------#
N <- 100 # sample size
K <- 100 # grid length
ts <- seq(0,1,length.out = K) # evenly-spaced grid from 0 to 1 of length K
scheme = 1
sparsity <- if(scheme==1){
  2:5
}else if(scheme==2){
  2
}else{
  stop("Scheme must be either 1 or 2.")
}
sims <- 100 # number of simulations


#---------------#
## Model Setup ##
#---------------#
beta0 <- rep(1,K)
alpha_f <- function(x){rep(1,length(x))}
alpha <- alpha_f(ts)
beta1_f <- function(x){rep(1,length(x))}
beta1 <- beta1_f(ts)

delta0 <- rep(0,K)
delta1_f <- function(x){rep(0.5,length(x))}
delta1 <- delta1_f(ts)
delta2_f <- function(x){rep(1,length(x))}
delta2 <- delta2_f(ts)

theta0_f <- function(x){rep(1,length(x))}
theta0 <- theta0_f(ts)
theta1_f <- function(x){rep(0.5,length(x))}
theta1 <- theta1_f(ts)


sig.eps <- 1
sig.lam <- 0.5


##--------- X1 -------##
mu.x1 <- rep(0,K)
Cx1_f <- function(t,s,sig2=10,rho=0.5){ # Matern covariance function with nu = 5/2
  d <- abs(outer(t,s,"-"))
  sig2*(1+sqrt(5)*d/rho + 5*d^2/(3*rho^2))*exp(-sqrt(5)*d/rho)}
Cx1 <- Cx1_f(ts,ts)

##-------- Z --------##
muz <- rep(0,K) # mean function of Z
Cz_f <- function(t,s,kappa = 1,V = 2){# Rational quadratic covariance function
  d <- abs(outer(t,s,"-"))
  # tmp <- exp(-(d/V))
  (1 + (d/(sqrt(2*V)*kappa))^2)^(-V)
}
Cz <- 0.5*Cz_f(ts,ts)
# Cz <- Cx1

##------- W --------##
mu.w <- delta0 + delta1*(theta0 + theta1*mu.x1) + delta2*muz


#--------------------#
## Simulation Setup ##
#--------------------#
names <- c("naive.fcr","naive.fcreg","fc2sls.pace","fc2sls.face",
           "fc2sls.fcr","fc2sls.lls","fc2sls.FCReg")
msedf <- data.frame(matrix(NA,sims,7,dimnames = list(NULL,names)))
betalist <- vector(mode = "list",length = length(names))
names(betalist) <- names
betas <- vector(mode = "list",length = 3)
names(betas) <- paste0("s",1:3)

#--------------------------------------------#
## Simulate observed data: 1) Random Design ##
#--------------------------------------------#
set.seed(123)

for(s in 1:sims){
  ##-------- Epsilon* ------##
  eps.star <- matrix(rnorm(N*K,mean = 0,sd = sig.eps),N,K)
  eps.star <- t(apply(eps.star,1,scale))
  
  ##-------- lam ------------##
  lam <- matrix(rnorm(N*K,mean = 0,sd = sig.lam),N,K)
  lam <- t(apply(lam,1,scale))
  
  ##-------- X1 --------##
  X1 <- mvrnorm(N,mu = mu.x1,Sigma = Cx1)
  
  ##-------- Z --------##
  Z <- mvrnorm(N,mu = muz,Sigma = Cz)
  
  ##-------- W* --------##
  W.star <- t(theta0 + theta1*t(X1)) + lam
  
  ##-------- W --------##
  W <- t(delta0 + delta1*t(W.star) + delta2*t(Z))
  
  ##-------- Y --------##
  Y <- t(beta0 + alpha*t(W) + beta1*t(X1)) + eps.star
  
  ##-------- Epsilon --------##
  eps <- t(beta1*t(X1)) + eps.star
  
  ## Check conditions for IV
  covwz <- covwe <- covze <- covzx1 <- covzes <- rep(NA,K)
  for(i in 1:K) covwz[i] <- mean((W[,i] - mu.w[i])*(Z[,i] - muz[i]))
  for(i in 1:K) covwe[i] <- mean((W[,i] - mu.w)*(eps[,i]))
  for(i in 1:K) covze[i] <- mean((Z[,i] - muz[i])*(eps[,i]))
  for(i in 1:K) covzx1[i] <- mean((Z[,i] - muz[i])*(X1[,i] - mu.x1[i]))
  for(i in 1:K) covzes[i] <- mean((Z[,i] - muz[i])*(eps.star[,i]))
  
  data_full <- data.frame(subj = rep(1:N,each = K),
                          argvals = rep(ts,N),
                          Y = c(t(Y)),
                          W = c(t(W)),
                          X1 = c(t(X1)),
                          Z = c(t(Z)))
  data_obs <- sampscheme(scheme = scheme,sparsity = sparsity,data_full)
  ts_obs <- sort(unique(data_obs$argvals))
  
  Lt <- split(data_obs$argvals,data_obs$subj)
  LY <- split(data_obs$Y,data_obs$subj)
  LW <- split(data_obs$W,data_obs$subj)
  LX1 <- split(data_obs$X1,data_obs$subj)
  LZ <- list(LZ1 = split(data_obs$Z,data_obs$subj))
  
  obs_list_naive <- list(
    W = list(Lt = Lt, Ly = LW),
    Y = list(Lt = Lt, Ly = LY)
  )
  
  
  #--------------------#
  ## Naive Estimation ##
  #--------------------#
  nphi <- floor((nrow(data_obs) - 21)/N)
  fcr_naive <- fcr(Y ~ s(argvals,k = -1,bs = "ps") +
                     s(argvals,by = W,k = -1,bs = "ps"),
                   argvals = "argvals",subj = "subj",
                   nPhi = nphi,
                   use_bam = T,         # use especially for larger N
                   # discrete = T,      # may also help speed up computation?
                   data = data_obs)
  fcr_naive$alpha <- get_modelterm(fcr_naive$fit,
                                   select = 1,
                                   n.grid = K,
                                   print.summary = F)$fit
  msedf$naive.fcr[s] <- sum((alpha - fcr_naive$alpha)^2)/K
  
  # cvdat <- FCReg_cv(dat = data_obs[,-c(5,6)],ts = ts,kFolds = 10)
  # bw <- cvdat$h
  bw <- 0.25
  FCReg_naive <- FCReg(vars = obs_list_naive,
                       userBwMu = bw,
                       userBwCov = bw,
                       outGrid = ts,
                       measurementError = F)
  
  msedf$naive.fcreg[s] <- sum((alpha - FCReg_naive$beta)^2)/K
  
  
  #---------------------#
  ## FC2SLS Estimation ##
  #---------------------#
  
  # ## Use fpca.sc to estimate full curves
  # ydata = data_obs[,1:3]; names(ydata) <- c(".id",".index",".value")
  # fit.y <- fpca.sc(ydata = ydata)
  # ### can't set argvals, so the sample must cover the entire grid
  # ### or predictions will only cover whatever section of the grid was observed
  
  
  fc2sls_imp_face <- fc2sls(LY,Lt,LW,LZ,grid.out = ts,fcmethod = "imp-face")
  msedf$fc2sls.face[s] <- mse(ts,ts_obs,alpha,fc2sls_imp_face$beta1)
  
  fc2sls_imp_pace <- fc2sls(LY,Lt,LW,LZ,grid.out = ts,fcmethod = "imp-pace")
  msedf$fc2sls.pace[s] <- mse(ts,ts_obs,alpha,fc2sls_imp_pace$beta1)
  
  # fc2sls_fcr <- fc2sls(LY,Lt,LW,LZ,grid.out = ts,fcmethod = "fcr")
  # msedf$fc2sls.fcr[s] <- mse(ts,ts_obs,alpha,fc2sls_fcr$alpha)
  
  # fc2sls_lls <- fc2sls(LY,Lt,LW,LZ,grid.out = ts,fcmethod = "loclin")
  # msedf$fc2sls.lls[s] <- mse(ts,ts_obs,alpha,fc2sls_lls$alpha)
  
  # fc2sls_FCReg <- fc2sls(LY,Lt,LW,LZ,h1 = 0.25,h2 = 0.25,grid.out = ts,fcmethod = "FCReg")
  # msedf$fc2sls.FCReg[s] <- mse(ts,ts_obs,alpha,fc2sls_FCReg$alpha)
  
  
  
  b1vals <- c(alpha,fcr_naive$alpha,FCReg_naive$beta[1,],
              fc2sls_imp_pace$beta1,
              # fc2sls_lls$alpha,
              # fc2sls_FCReg$alpha,
              # fc2sls_fcr$alpha
              fc2sls_imp_face$beta1
  )
  ylim <- c(min(b1vals),max(b1vals))
  plot(ts,alpha,type = 'l',lwd = 2,col = 'blue',ylim = ylim)
  lines(ts,fcr_naive$alpha,col = 'red')
  lines(ts,FCReg_naive$beta[1,],col = 'red',lty = 2)
  lines(ts,fc2sls_imp_face$beta1,col = 'black')
  lines(ts,fc2sls_imp_pace$beta1,col = 'black',lwd = 2)
  # lines(ts,fc2sls_fcr$alpha,lty = 2)
  # if(length(fc2sls_lls$alpha)!=length(alpha)){
  #   lines(ts_obs,fc2sls_lls$alpha,lty = 3)
  # }else{
  #   lines(ts,fc2sls_lls$alpha,lty = 3)
  # }
  # lines(ts,fc2sls_FCReg$alpha,lty = 4)
  legend("topright",
         legend = c("truth",names[-c(5:7)]),
         col = c("blue","red","red",rep("black",2)),
         lty = c(1,1,2,1,1),
         lwd = c(2,1,1,2,1))
  
  # Store betas from first 3 simulations to plot
  if(s <= 3){
    betalist[["naive.fcr"]] <- fcr_naive$alpha
    betalist[["naive.fcreg"]] <- c(FCReg_naive$beta)
    betalist[["fc2sls.pace"]] <- c(fc2sls_imp_pace$beta1)
    betalist[["fc2sls.face"]] <- c(fc2sls_imp_face$beta1)
    # betalist[["fc2sls.fcr"]] <- fc2sls_fcr$alpha
    # betalist[["fc2sls.lls"]] <- fc2sls_lls$alpha
    # betalist[["fc2sls.FCReg"]] <- fc2sls_FCReg$alpha
    betas[[s]] <- betalist
  }
  
  
  if(s %% 10 == 0){
    s
  }
}

colMeans(msedf)
apply(msedf,2,median)

results <- list(msedf = msedf,betas = betas,N = N,K = K,ts = ts,
                sims = sims,scheme = scheme,sparsity = sparsity,
                beta0 = beta0,alpha = alpha,beta1 = beta1,
                delta0 = delta0,delta1 = delta1, delta2 = delta2,
                theta0 = theta0, theta1 = theta1,
                sig.eps = sig.eps, sig.lam = sig.lam,
                Cx1 = Cx1, Cz = Cz)
save(results,file = "mse-N100-s1scheme1.RData")





# ### Leave One Out CV
# hl <- 20
# h_vec <- 10^(seq(from = 0,to = -4,length.out = hl))
# # cv_vec <- rep(NA,length = hl)
# cv <- matrix(0,ncol = hl,nrow = N)
# 
# for(hi in 1:hl){
#   h <- h_vec[hi]
# 
#   ## select curve i to leave out
#   for(i in 1:N){
#     dat_cv <- data_obs[-which(data_obs$subj==i),]
#     # muhat_cv <- llsmooth_mean(K = K,Z.tilde = dat_cv$Z,t.tilde = dat_cv$argvals,h = h)
#     muhat_cv <- GetMeanCurve(Ly = LZ,Lt = Lt,optns = list(userBwMu = h,
#                                                           dataType = "Sparse"))
#     muhat_cv <- ConvertSupport(fromGrid = muhat_cv$workGrid,toGrid = ts,mu = muhat_cv$mu)
#     
#     # CV contribution
#     rows <- which(data_obs$subj==i)
#     Kj <- length(rows)
#     tj <- data_obs$argvals[rows]
#     muhat_tj <- muhat_cv[match(tj,ts)]
#     cv[i,hi] <- sum((data_obs$Z[rows] - muhat_tj)^2)
#   }
# }
# 
# h <- h_vec[which.min(colMeans(cv))]
