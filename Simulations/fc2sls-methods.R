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
