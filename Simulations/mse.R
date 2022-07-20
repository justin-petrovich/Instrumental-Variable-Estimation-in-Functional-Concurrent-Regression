mse <- function(grid,grid.obs,beta,betahat){
  if(length(beta)!=length(betahat)){
    notobs <- setdiff(grid,grid.obs)
    id <- match(notobs,grid)
    beta <- beta[-id]
  }
  K <- length(beta)
  sum((beta - betahat)^2)/K
}
