#-------------#
## FC - 2SLS ##
#-------------#

# fcmethod is one of "loclin", "FCReg", "fcr", "imp-face", or "imp-pace".
# "loclin" uses local linear smoothing for both stages and is coded by hand
# "FCReg" and "fcr" use the corresponding functions at each stage of the regression.
# "imp-face" imputes the functions using face, then uses pffr on the dense curves
# "imp-pace" imputes the functions using pace, then uses pffr on the dense curves

fc2sls <- function(LY,Lt,LX,Zvars,h1,h2,grid.out,fcmethod = "imp-face"){

  grid.obs <- sort(unique(unlist(Lt)))
  N <- length(LX)
  
  method <- fcmethod %in% c("loclin","FCReg","fcr","imp-face","imp-pace")
  if(!method)
    stop("Must choose a method from loclin, FCReg, fcr, imp-face, or imp-pace.")
  
  # vars.l <- sapply(vars,length)
  # if(length(unique(vars.l)!=1))
  #   stop("Functions are assumed to be observed concurrently, so each function must have the same
  #        number of observed time points.")
  
  if(fcmethod=="FCReg"){
    est <- FCReg_2s(Lt,LX,Zvars,LY,h1,h2,grid.out,N)
    
  }else if(fcmethod=="imp-pace"){
    est <- imp_pace(Lt,LX,Zvars,LY,grid.out)
    
  }else if(fcmethod=="loclin"){
    est <- loclin_2s(Lt,LX,LZ,LY,grid.obs,N)
    
  }else{
    
    dat <- data.frame(
      subj = rep(1:N,times = unlist(lapply(Lt,length))),
      argvals = unlist(Lt),
      Y = unlist(LY),
      X = unlist(LX)
    )

    Zdat <- data.frame(sapply(Zvars,unlist))
    dat <- cbind(dat,Zdat)
    zcols <- match(names(Zdat),colnames(dat))

    if(fcmethod=="fcr"){
      est <- fcr_2s(dat,zcols,N,grid.out)
      
    }else if(fcmethod=="imp-face"){  
      est <- imp_face(dat,zcols,N,grid.out)
      
    }
  }
  
  return(est)
}