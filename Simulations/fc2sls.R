#-------------#
## FC - 2SLS ##
#-------------#

# fcmethod is one of "imp-face" or "imp-pace".
# "imp-face" imputes the functions using face, then uses pffr on the dense curves
# "imp-pace" imputes the functions using pace, then uses pffr on the dense curves

fc2sls <- function(LY,Lt,LX,Zvars,h1,h2,grid.out,fcmethod = "imp-face"){

  grid.obs <- sort(unique(unlist(Lt)))
  N <- length(LX)
  
  method <- fcmethod %in% c("imp-face","imp-pace")
  if(!method)
    stop("Must choose a method from imp-face or imp-pace.")
  
  # vars.l <- sapply(vars,length)
  # if(length(unique(vars.l)!=1))
  #   stop("Functions are assumed to be observed concurrently, so each function must have the same
  #        number of observed time points.")
  
  if(fcmethod=="imp-pace"){
    est <- imp_pace(Lt,LX,Zvars,LY,grid.out)
    
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

    est <- imp_face(dat,zcols,N,grid.out)
      
  }
  
  return(est)
}