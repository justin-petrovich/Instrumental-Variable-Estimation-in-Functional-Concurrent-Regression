library(fcr)
library(MASS)
library(tidyverse)
library(refund)
library(itsadug)
library(fdapace)
library(face)

setwd("C:/Users/Owner/OneDrive - Saint Vincent College/Research/Time-Varying Elasticities/FCR-IV/Simulation Results")
source("fc2sls.R")
source("fc2sls-methods.R")
source("sampling scheme.R")
source("mse.R")


#-------------------------#
## Simulation Parameters ##
#-------------------------#
N <- 200 # sample size
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
beta1_f <- function(x){1 + x*1*sin(2*pi*x)}
beta1 <- beta1_f(ts)
beta2_f <- function(x){5 + x*5*cos(2*pi*x)}
beta2 <- beta2_f(ts)

delta0 <- rep(0,K)
delta1_f <- function(x){1 + 0.5*sin(2*pi*x)}
delta1 <- delta1_f(ts)
delta2_f <- function(x){1 + sin(2*pi*x)}
delta2 <- delta2_f(ts)
delta3_f <- function(x){1 + cos(2*pi*x)}
delta3 <- delta3_f(ts)

theta0_f <- function(x){1 + sin(2*pi*x)}
theta0 <- theta0_f(ts)
theta1_f <- function(x){1 + 0.5*cos(2*pi*x)}
theta1 <- theta1_f(ts)


sig.eps <- 1
sig.nu <- 0.5

#--------------#
## Data Setup ##
#--------------#

##--------- X2 -------##
mu.x2 <- rep(0,K)
Cx2_f <- function(t,s,sig2=1,rho=0.5){ # Matern covariance function with nu = 5/2
  d <- abs(outer(t,s,"-"))
  sig2*(1+sqrt(5)*d/rho + 5*d^2/(3*rho^2))*exp(-sqrt(5)*d/rho)}
Cx2 <- Cx2_f(ts,ts)

##-------- Z1 --------##
mu.z1 <- rep(0,K) # mean function of Z
Cz1 <- Cx2

##-------- Z2 --------##
# mu.z2 <- (1/mean(ts)^4)*(mean(ts)^2-(ts - mean(ts))^2)
# Cz_f <- function(t,s,kappa = 1,V = 2){# Rational quadratic covariance function
#   d <- abs(outer(t,s,"-"))
#   # tmp <- exp(-(d/V))
#   (1 + (d/(sqrt(2*V)*kappa))^2)^(-V)
# }
# Cz2 <- Cz_f(ts,ts)

# ##--------X1----------##
# mu.x1 <- delta0 + delta1*(theta0 + theta1*mu.x2) + delta2*muz

#--------------------#
## Simulation Setup ##
#--------------------#
names <- c("naive.fcr","naive.fcreg","naive.pace","fc2sls.pace","fc2sls.face",
           "fc2sls.fcr","fc2sls.lls","fc2sls.FCReg")
msedf <- data.frame(matrix(NA,sims,8,dimnames = list(NULL,names)))
betalist <- vector(mode = "list",length = length(names))
names(betalist) <- names
betas <- vector(mode = "list",length = 3)
names(betas) <- paste0("s",1:3)

set.seed(123)

#--------------------------------------------#
## Simulate observed data: 1) Random Design ##
#--------------------------------------------#
for(s in 1:sims){
  ##-------- Epsilon* ------##
  eps.star <- matrix(rnorm(N*K,mean = 0,sd = sig.eps),N,K)
  eps.star <- t(apply(eps.star,1,scale))
  
  ##-------- Nu ------------##
  nu <- matrix(rnorm(N*K,mean = 0,sd = sig.nu),N,K)
  nu <- t(apply(nu,1,scale))
  
  ##-------- Z1 --------##
  Z1 <- mvrnorm(N,mu = mu.z1,Sigma = Cz1)
  
  ##-------- Z2 --------##
  # Z2 <- mvrnorm(N,mu = mu.z2,Sigma = Cz2)
  Z2 <- Z1^2
  
  ##--------- X2 -------##
  X2 <- mvrnorm(N,mu = mu.x2,Sigma = Cx2)
  
  ##-------- X* --------##
  X.star <- t(theta0 + theta1*t(X2)) + nu
  
  ##-------- X1 --------##
  X1 <- t(delta0 + delta1*t(X.star) + delta2*t(Z1) + delta3*t(Z2))
  
  ##-------- Y --------##
  Y <- t(beta0 + beta1*t(X1) + beta2*t(X2)) + eps.star
  
  data_full <- data.frame(subj = rep(1:N,each = K),
                          argvals = rep(ts,N),
                          Y = c(t(Y)),
                          X1 = c(t(X1)),
                          X2 = c(t(X2)),
                          Z1 = c(t(Z1)),
                          Z2 = c(t(Z2)))
  data_obs <- sampscheme(scheme = scheme,sparsity = sparsity,data_full)
  ts_obs <- sort(unique(data_obs$argvals))
  
  Lt <- split(data_obs$argvals,data_obs$subj)
  LY <- split(data_obs$Y,data_obs$subj)
  LX1 <- split(data_obs$X1,data_obs$subj)
  LX2 <- split(data_obs$X2,data_obs$subj)
  LZ1 <- split(data_obs$Z1,data_obs$subj)
  LZ2 <- split(data_obs$Z2,data_obs$subj)
  # Lt.d <- split(data_full$argvals,data_full$subj)
  # Ly.d <- split(data_full$Y,data_full$subj)
  # LX1.d <- split(data_full$X1,data_full$subj)
  # LX2.d <- split(data_full$X2,data_full$subj)
  # LZ.d <- split(data_full$Z,data_full$subj)
  
  obs_list_naive <- list(
    X1 = list(Lt = Lt, Ly = LX1),
    Y = list(Lt = Lt, Ly = LY)
  )
  
  Zvars <- list(LZ1 = LZ1, LZ2 = LZ2)
  
  #--------------------#
  ## Naive Estimation ##
  #--------------------#
  nphi <- floor((nrow(data_obs) - 21)/N)
  fcr_naive <- fcr(Y ~ s(argvals,k = -1,bs = "ps") +
                     s(argvals,by = X1,k = -1,bs = "ps"),
                   argvals = "argvals",subj = "subj",
                   nPhi = nphi,
                   use_bam = T,         # use especially for larger N
                   # discrete = T,      # may also help speed up computation?
                   data = data_obs)
  fcr_naive$beta1 <- get_modelterm(fcr_naive$fit,
                                   select = 1,
                                   n.grid = K,
                                   print.summary = F)$fit
  msedf$naive.fcr[s] <- sum((beta1 - fcr_naive$beta1)^2)/K
  
  # cvdat <- FCReg_cv(dat = data_obs[,-c(5,6)],ts = ts,kFolds = 10)
  # bw <- cvdat$h
  bw <- 0.25
  FCReg_naive <- FCReg(vars = obs_list_naive,
                       userBwMu = bw,
                       userBwCov = bw,
                       outGrid = ts,
                       measurementError = F)
  
  msedf$naive.fcreg[s] <- sum((beta1 - FCReg_naive$beta)^2)/K
  
  
  # Ypace <- FPCA(Ly = LY,Lt = Lt,optns = list(methodBwCov = "GCV",
  #                                            methodBwMu = "GCV",
  #                                            dataType = "Sparse",
  #                                            error = TRUE,
  #                                            usergrid = FALSE,
  #                                            nRegGrid = K,
  #                                            methodRho = "trunc"))
  # 
  # X1pace <- FPCA(Ly = LX1,Lt = Lt,optns = list(methodBwCov = "GCV",
  #                                              methodBwMu = "GCV",
  #                                              dataType = "Sparse",
  #                                              error = TRUE,
  #                                              usergrid = FALSE,
  #                                              nRegGrid = K,
  #                                              methodRho = "trunc"))
  # 
  # Ydense <- fitted(Ypace)
  # X1dense <- fitted(X1pace)
  # pacedat <- list(Y = Ydense,X1 = X1dense)
  # naive.pace <- pffr(Y ~ X1,data = pacedat,yind = ts)
  # naive.pace$beta1 <- coef(naive.pace)$smterms[[2]][['value']]
  # msedf$naive.pace[s] <- sum((beta1 - naive.pace$beta1)^2)/K
  
  #---------------------#
  ## FC2SLS Estimation ##
  #---------------------#
  
  # ## Use fpca.sc to estimate full curves
  # ydata = data_obs[,1:3]; names(ydata) <- c(".id",".index",".value")
  # fit.y <- fpca.sc(ydata = ydata)
  # ### can't set argvals, so the sample must cover the entire grid
  # ### or predictions will only cover whatever section of the grid was observed
  
  
  fc2sls_imp_face <- fc2sls(LY,Lt,LX1,Zvars,Lt,grid.out = ts,fcmethod = "imp-face")
  msedf$fc2sls.face[s] <- mse(ts,ts_obs,beta1,fc2sls_imp_face$beta1)
  
  fc2sls_imp_pace <- fc2sls(LY,Lt,LX1,Zvars,Lt,grid.out = ts,fcmethod = "imp-pace")
  msedf$fc2sls.pace[s] <- mse(ts,ts_obs,beta1,fc2sls_imp_pace$beta1)
  
  # fc2sls_fcr <- fc2sls(LY,Lt,LX1,Zvars,grid.out = ts,fcmethod = "fcr")
  # msedf$fc2sls.fcr[s] <- mse(ts,ts_obs,beta1,fc2sls_fcr$beta1)
  
  # fc2sls_lls <- fc2sls(LY,LX1,LZ,Lt,grid.out = ts,fcmethod = "loclin")
  # msedf$fc2sls.lls[s] <- mse(ts,ts_obs,beta1,fc2sls_lls$beta1)
  
  # fc2sls_FCReg <- fc2sls(LY,Lt,LX1,Zvars,h1 = 0.25,h2 = 0.25,grid.out = ts,fcmethod = "FCReg")
  # msedf$fc2sls.FCReg[s] <- mse(ts,ts_obs,beta1,fc2sls_FCReg$beta1)
  
  b1vals <- c(beta1,
              fcr_naive$beta1,
              FCReg_naive$beta[1,],
              fc2sls_imp_pace$beta1,
              # fc2sls_lls$beta1,
              # fc2sls_FCReg$beta1,
              # fc2sls_fcr$beta1,
              fc2sls_imp_face$beta1
  )
  ylim <- c(min(b1vals),max(b1vals))
  plot(ts,beta1,type = 'l',lwd = 2,col = 'blue',ylim = ylim)
  lines(ts,fcr_naive$beta1,col = 'red')
  lines(ts,FCReg_naive$beta[1,],col = 'red',lty = 2)
  # lines(ts,naive.pace$beta1,col = 'red',lty = 3)
  lines(ts,fc2sls_imp_face$beta1,col = 'black')
  lines(ts,fc2sls_imp_pace$beta1,col = 'black',lwd = 2)
  # lines(ts,fc2sls_fcr$beta1,lty = 2)
  # if(length(fc2sls_lls$beta1)!=length(beta1)){
  #   lines(ts_obs,fc2sls_lls$beta1,lty = 3)
  # }else{
  #   lines(ts,fc2sls_lls$beta1,lty = 3)
  # }
  # lines(ts,fc2sls_FCReg$beta1,lty = 4)
  legend("topright",
         legend = c("truth",names[-c(1,6,7,8)]),
         col = c("blue","red","red",rep("black",2)),
         lty = c(1,2,3,1,1),
         lwd = c(2,1,1,2,1))
  
  # Store betas from first 3 simulations to plot
  if(s <= 3){
    betalist[["naive.fcr"]] <- fcr_naive$beta1
    betalist[["naive.fcreg"]] <- c(FCReg_naive$beta)
    # betalist[["naive.pace"]] <- c(naive.pace$beta1)
    betalist[["fc2sls.pace"]] <- c(fc2sls_imp_pace$beta1)
    betalist[["fc2sls.face"]] <- c(fc2sls_imp_face$beta1)
    # betalist[["fc2sls.fcr"]] <- fc2sls_fcr$beta1
    # betalist[["fc2sls.lls"]] <- fc2sls_lls$beta1
    # betalist[["fc2sls.FCReg"]] <- fc2sls_FCReg$beta1
    betas[[s]] <- betalist
  }

  if(s %% 10 ==0){
    cat(s,"\n")
  }

}

colMeans(msedf)
results <- list(msedf = msedf,betas = betas,N = N,K = K,ts = ts,
                sims = sims,scheme = scheme,sparsity = sparsity,
                beta0 = beta0,beta1 = beta1,beta2 = beta2,
                delta0 = delta0,delta1 = delta1, delta2 = delta2,
                theta0 = theta0, theta1 = theta1,
                sig.eps = sig.eps, sig.nu = sig.nu,
                Cx2 = Cx2, Cz1 = Cz1)
save(results,file = "mse-N200-s3scheme1.RData")
