Simu_Multi_Norm <- function(x_len, sd = 1, pho = 0.5) {
  V <- matrix(data = NA,
              nrow = x_len,
              ncol = x_len)
  for (i in 1:x_len) {
    for (j in 1:x_len) {
      V[i, j] <- pho ^ abs(i - j)
    }
  }
  V <- (sd ^ 2) * V
  return(V)
}

#' @title Simulation model 1
#'
#' @description similar to the simulation model in SCAD (Jianqing Fan and Runze Li,2001)
#' @import MASS
#' @param n sample size.
#' @param sigma the standard deviation of the error used in the simulation model.
#' @param beta the parameters needed to be estimated which is assumed in the simulation model.
#'
#' @importFrom stats rnorm
#' @return the data for similation 
#' @examples
#' \dontrun{
#' data<-simulation1(n=40,sigma=1,beta=c(3, 1.5, 2, 0, 0, 0, 0, 0))
#' }
#' @export
simulation1 <- function(n=40,sigma=1,beta = c(3, 1.5, 2, 0, 0, 0, 0, 0)){
  d <- length(beta)
  x <- mvrnorm(n, mu = rep(0, d), Simu_Multi_Norm(x_len = d, sd  = 1, pho = 0.5))
  epsilon  <-  rnorm(n, 0, sigma)
  y  <-  x %*% beta + epsilon
  return(list(y=y,x=x))
}