library(face)
library(fdapace)
library(refund)
library(itsadug)

imp_face <- function(dat,zcols,N,grid.out){
  #### !! Need to make sure that dat is in the correct form
  #### Add some checks/warnings here
  
  K <- length(grid.out)
  
  ## use FACE to estimate full curves
  ydata <- dat[,1:3]; names(ydata)[3] <- "y"
  x1data <- dat[,c(1,2,4)]; names(x1data)[3] <- "y"

  newdata <- data.frame(subj = rep(1:N,each = K),
                        argvals = rep(grid.out,times = N),
                        y = rep(NA,N*K))
  fit.y <- face.sparse(data = ydata,newdata = rbind(ydata,newdata),
                       argvals.new = grid.out,knots = 7)
  fit.x <- face.sparse(data = x1data,newdata = rbind(x1data,newdata),
                       argvals.new = grid.out,knots = 7)
  
  Ydense <- matrix(fit.y$y.pred[-c(1:nrow(ydata))],nrow = N,ncol = K,byrow = T)
  X1dense <- matrix(fit.x$y.pred[-c(1:nrow(x1data))],nrow = N,ncol = K,byrow = T)
  
  data1 <- vector(mode = "list",length = 1 + length(zcols))
  data1[[1]] <- X1dense
  znames <- paste0("Z",1:length(zcols))
  names(data1) <- c("X1",znames)
  
  for(i in 1:length(zcols)){
    zdata <- dat[,c(1,2,zcols[i])]; names(zdata)[3] <- "y" 
    fit.z <- face.sparse(data = zdata,newdata = rbind(zdata,newdata),
                         argvals.new = grid.out,knots = 7)
    Zdense <- matrix(fit.z$y.pred[-c(1:nrow(zdata))],nrow = N,ncol = K,byrow = T)
    data1[[i+1]] <- Zdense
  }
  
  ## Use functional concurrent models for densely sampled data
  # Stage 1
  s1formula <- as.formula(paste("X1 ~ ",paste(znames, collapse = "+")))
  stage1 <- pffr(formula = s1formula,data = data1,yind = grid.out)
  X1hat <- fitted(stage1)
  
  # Stage 2
  data2 <- list(Y = Ydense,X1hat = X1hat)
  stage2 <- pffr(Y ~ X1hat,data = data2,yind = grid.out)
  beta1hat_2sls <- coef(stage2)$smterms[[2]][['value']]
  beta0 <- NA
  
  out <- list(beta0 = beta0,
              beta1 = beta1hat_2sls)
  return(out)
}





imp_pace <- function(Lt,LX,Zvars,LY,grid.out){
  
  K <- length(grid.out)
  
  ## Use PACE to estimate full curves
  Ypace <- FPCA(Ly = LY,Lt = Lt,optns = list(methodBwCov = "GCV",
                                             methodBwMu = "GCV",
                                             dataType = "Sparse",
                                             error = TRUE,
                                             usergrid = FALSE,
                                             nRegGrid = K,
                                             methodRho = "trunc"))
  
  X1pace <- FPCA(Ly = LX,Lt = Lt,optns = list(methodBwCov = "GCV",
                                              methodBwMu = "GCV",
                                              dataType = "Sparse",
                                              error = TRUE,
                                              usergrid = FALSE,
                                              nRegGrid = K,
                                              methodRho = "trunc"))
  
  Ydense <- fitted(Ypace)
  X1dense <- fitted(X1pace)
  
  numz <- length(Zvars)
  if(numz >= length(Lt)){
    stop("Dimension of Zvars is too high")
  }else if(numz > 100){
    warning("Dimension of Zvars is large. Make sure that Zvars is a list of instrumental variables
            such that each variable is itself a list containing the observed values.")
  }
  data1 <- vector(mode = "list",length = 1 + numz)
  data1[[1]] <- X1dense
  znames <- paste0("Z",1:numz)
  names(data1) <- c("X1",znames)
  
  for(i in 1:numz){
    Zpace <- FPCA(Ly = Zvars[[i]],
                  Lt = Lt,
                  optns = list(methodBwCov = "GCV",
                               methodBwMu = "GCV",
                               dataType = "Sparse",
                               error = TRUE,
                               usergrid = FALSE,
                               nRegGrid = K,
                               methodRho = "trunc"))
    Zdense <- fitted(Zpace)
    data1[[i + 1]] <- Zdense
  }
  
  
  ## Use functional concurrent models for densely sampled data
  # Stage 1
  s1formula <- as.formula(paste("X1 ~ ",paste(znames, collapse = "+")))
  stage1 <- pffr(formula = s1formula,data = data1,yind = grid.out)
  X1hat <- fitted(stage1)
  
  # Stage 2
  data2 <- list(Y = Ydense,X1hat = X1hat)
  stage2 <- pffr(Y ~ X1hat,data = data2,yind = grid.out)
  beta1hat_2sls <- coef(stage2)$smterms[[2]][['value']]
  beta0 <- NA
  
  
  out <- list(beta0 = beta0,
              beta1 = beta1hat_2sls)
  return(out)
}





FCReg_2s <- function(Lt,LX,Zvars,LY,h1,h2,grid.out,N){
  ## Stage 1:
  numz <- length(Zvars)
  vars <- vector(mode = "list",length = 1 + numz)
  vars[[1]] <- list(Ly = LX,Lt = Lt)
  for(i in 2:(numz+1)){
    vars[[i]] <- list(Ly = Zvars[[i-1]],Lt = Lt)
  }
  znames <- paste0("Z",1:numz)
  names(vars) <- c("X1",znames)
  
  FCmod1 <- FCReg(vars = vars, h1, h1, outGrid = grid.out)
  
  LXhat <- vector(mode = "list",length = N)
  for(i in 1:N){
    gid <- match(Lt[[i]],grid.out)
    LXhat[[i]] <- FCmod1$beta0[gid] + 
      FCmod1$beta[1,gid]*vars[[1]][["Ly"]][[i]]
  }
  
  ## Stage 2:
  vars <- list(Xhat = list(Ly = LXhat,Lt = Lt),
               Y = list(Ly = LY,Lt = Lt))
  FCmod2 <- FCReg(vars = vars, h2, h2, outGrid = grid.out)
  
  beta0star <- FCmod2$beta0
  beta_2sls <- FCmod2$beta[1,]
  
  out <- list(beta0 = beta0star,
              beta1 = beta_2sls)
  return(out)
}





fcr_2s <- function(dat,zcols,N,grid.out){
  nphi <- floor((nrow(data_obs) - 21)/N)
  K <- length(grid.out)
  
  znum <- length(zcols)
  ## Stage 1:
  zsmooths <- paste("s(argvals, by = LZ",1:znum,", k = -1, bs = 'ps')",sep = "",collapse = " + ")
  s1formula <- as.formula(paste("X ~ s(argvals,k = -1,bs = 'ps') + ",zsmooths))
  fcrmod1 <- fcr(formula = s1formula,
                 use_bam = T,
                 argvals = "argvals",subj = "subj",data = dat,nPhi = nphi)
  m1pred <- predict(fcrmod1,newdata = dat)
  dat$Xhat <- m1pred$insample_predictions
  
  ##Stage 2:
  fcrmod2 <- fcr(Y ~ s(argvals,k = -1,bs = "ps") +
                   s(argvals,by = Xhat,k = -1,bs = "ps"),
                 use_bam = T,
                 argvals = "argvals",subj = "subj",data = dat,nPhi = nphi)

  beta0star <- NA
  beta_2sls <- get_modelterm(fcrmod2$fit,
                             select = 1,
                             n.grid = K,
                             print.summary = F)$fit
  
  out <- list(beta0 = beta0star,
              beta1 = beta_2sls)
  return(out)
}





loclin_2s <- function(Lt,LX,LZ,LY,grid.obs,N){
  ## STAGE 1:
  ## 1) Estimate mean functions
  # muhatZ <- llsmooth_mean(K = K,Z.tilde = dat$Z,t.tilde = dat$argvals,h = h)
  # muhatX <- llsmooth_mean(K = K,Z.tilde = dat$X,t.tilde = dat$argvals,h = h)
  muhatZ <- GetMeanCurve(Ly = LZ,Lt = Lt,optns = list(methodBwMu = "GCV",
                                                      dataType = "Sparse",
                                                      kernel = "gauss"))
  muhatZ <- ConvertSupport(fromGrid = muhatZ$workGrid,toGrid = grid.obs,mu = muhatZ$mu)
  muhatX <- GetMeanCurve(Ly = LX,Lt = Lt,optns = list(methodBwMu = "GCV",
                                                      dataType = "Sparse",
                                                      kernel = "gauss"))
  muhatX <- ConvertSupport(fromGrid = muhatX$workGrid,toGrid = grid.obs,mu = muhatX$mu)
  
  ## 2) and 3) Estimate smooth covariance/cross-covariance functions
  
  # Gzz_smooth <- llsmooth_cov(N = N,K = K,LX = LZ,LZ = LZ,Lt = Lt,h = h,
  #                            muhatX = muhatZ,muhatZ = muhatZ,ts = grid.obs)
  # Gxz_smooth <- llsmooth_cov(N = N,K = K,LX = LX,LZ = LZ,Lt = Lt,h = h,
  #                            muhatX = muhatX,muhatZ = muhatZ,ts = grid.obs)
  # fpca_z <- FPCA(Ly = LZ,Lt = Lt,optns = list(methodBwMu = "GCV",
  #                                             methodBwCov = "GCV",
  #                                             usergrid = TRUE,
  #                                             dataType = "Sparse",
  #                                             methodRho = "trunc"))
  # Gzz_smooth <- fpca_z$smoothedCov
  Gzz <- GetCrCovYX(Ly1 = LZ,Lt1 = Lt,Ymu1 = muhatZ,
                    Ly2 = LZ,Lt2 = Lt,Ymu2 = muhatZ,
                    rmDiag = TRUE)
  Gzz_smooth <- ConvertSupport(Cov = Gzz$smoothedCC,fromGrid = Gzz$smoothGrid[,1],toGrid = grid.obs)
  Gxz <- GetCrCovYX(Ly1 = LX,Lt1 = Lt,Ymu1 = muhatX,
                    Ly2 = LZ,Lt2 = Lt,Ymu2 = muhatZ,
                    rmDiag = TRUE)
  Gxz_smooth <- ConvertSupport(Cov = Gxz$smoothedCC,fromGrid = Gxz$smoothGrid[,1],toGrid = grid.obs)
  
  
  ## 4) Estimate regression coefficients
  beta1hat <- diag(Gxz_smooth)/diag(Gzz_smooth)
  beta0hat <- muhatX - beta1hat*muhatZ
  LXhat <- vector(mode = "list",length = N)
  for(i in 1:N){
    id <- match(Lt[[i]],grid.obs)
    LXhat[[i]] <- beta0hat[id] + beta1hat[id]*LZ[[i]]
  }
  # dat$Xhat <- unlist(LXhat)
  # 
  ## STAGE 2:
  ## 1) Estimate mean functions
  # muhatY <- llsmooth_mean(K = K,Z.tilde = dat$Y,t.tilde = dat$argvals,h = h)
  # muhatXhat <- llsmooth_mean(K = K,Z.tilde = dat$Xhat,t.tilde = dat$argvals,h = h)
  muhatY <- GetMeanCurve(Ly = LY,Lt = Lt,optns = list(methodBwMu = "GCV",
                                                      dataType = "Sparse",
                                                      kernel = "gauss"))
  muhatY <- ConvertSupport(fromGrid = muhatY$workGrid,toGrid = grid.obs,mu = muhatY$mu)
  muhatXhat <- GetMeanCurve(Ly = LXhat,Lt = Lt,optns = list(methodBwMu = "GCV",
                                                            dataType = "Sparse",
                                                            kernel = "gauss"))
  muhatXhat <- ConvertSupport(fromGrid = muhatXhat$workGrid,toGrid = grid.obs,mu = muhatXhat$mu)
  
  ## 2) and 3) Estimate smooth covariance/cross-covariance functions
  # Gxx_smooth <- llsmooth_cov(N = N,K = K,LX = LXhat,LZ = LXhat,Lt = Lt,h = h,
  #                            muhatX = muhatXhat,muhatZ = muhatXhat,ts = grid.obs)
  # Gxy_smooth <- llsmooth_cov(N = N,K = K,LX = LXhat,LZ = LY,Lt = Lt,h = h,
  #                            muhatX = muhatXhat,muhatZ = muhatY,ts = grid.obs)
  # fpca_Xhat <- FPCA(Ly = LXhat,Lt = Lt,optns = list(methodBwMu = "GCV",
  #                                                   methodBwCov = "GCV",
  #                                                   usergrid = TRUE,
  #                                                   dataType = "Sparse",
  #                                                   methodRho = "trunc"))
  # Gxx_smooth <- fpca_Xhat$smoothedCov
  Gxx <- GetCrCovYX(Ly1 = LXhat,Lt1 = Lt,Ymu1 = muhatXhat,
                    Ly2 = LXhat,Lt2 = Lt,Ymu2 = muhatXhat,
                    rmDiag = TRUE)
  Gxx_smooth <- ConvertSupport(Cov = Gzz$smoothedCC,fromGrid = Gzz$smoothGrid[,1],toGrid = grid.obs)
  Gxy <- GetCrCovYX(Ly1 = LXhat,Lt1 = Lt,Ymu1 = muhatXhat,
                    Ly2 = LY,Lt2 = Lt,Ymu2 = muhatY,
                    rmDiag = TRUE)
  Gxy_smooth <- ConvertSupport(Cov = Gxy$smoothedCC,fromGrid = Gxy$smoothGrid[,1],toGrid = grid.obs)
  
  
  ## 4) Estimate regression coefficients
  beta_2sls <- diag(Gxy_smooth)/diag(Gxx_smooth)
  beta0star <- muhatY - beta_2sls*muhatXhat
  
  out <- list(beta0 = beta0star,
              beta1 = beta_2sls)
  return(out)
}








