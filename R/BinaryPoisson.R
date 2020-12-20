#' @title A function to fit data with a mixed distribution of Poisson distribution using EM algorithm.
#' @description A function to use a mixed distribution of Poisson distribution with mean mu1 and Poisson distribution with mean mu1 to fit data, which is known as the number of occurrences of a random event in unit time.And a is the percentage of the first distribution.
#' @param x An event (such as death) is arranged in increasing order of the number of occurrences per day within a year, and there is no number of days to complete the maximum number of occurrences with 0.
#' @return  an estimate of BinaryPoisson parameter (a,mu1,mu2)
#' @examples
#' \dontrun{
#' x=c(162,267,271,185,111,61,27,8,3,1)
#' BinaryPoisson(x)
#' }
#' @export
BinaryPoisson<-function(x=c(162,267,271,185,111,61,27,8,3,1)){
  n <- length(x)
  d <- 0:(n-1)
  a <- 0.4 ; mu1 <- 1 ; mu2 <-3;
  t <- function(a,m1,m2){
    b <- x
    for(i in 1:n){
      b[i] <- a*exp(-m1)*m1^(i-1)/(a*exp(-m1)*m1^(i-1)+(1-a)*exp(-m2)*m2^(i-1))
    }
    b
  }
  for(j in 1:(100*n)){
    y <- x*t(a,mu1,mu2)
    z <- x-y
    a <- sum(y)/sum(x)
    mu1 <- d%*%y/sum(y)
    mu2 <- d%*%z/sum(z)
  }
  return(list(a=a,mu1=mu1,mu2=mu2))
}
