# Define the penalty derivative
p_lambda_d<-function(theta,a,lambda){
  l=length(theta)
  p1=rep(0,l)
  for(i in 1:l){
    if((abs(theta[i])>lambda)){
      if((a*lambda) >= abs(theta[i])){
        p1[i]=(a*lambda - theta[i])/((a - 1)*lambda)
      }
      else{
        p1[i]=0
      }
    }
    else{
      p1[i]=lambda
    }
  }
  return(p1)
}

Sigmaf=function(beta,a,lambda) {
  l=length(beta)
  Sigmabeta=matrix(0,l,l)
  diag(Sigmabeta)=p_lambda_d(abs(beta),a,lambda)/abs(beta)
  return (Sigmabeta)
}

scad=function(x,y,a,lambda) {
  n = length(y)
  d = dim(x)[2]
  beta0=solve(t(x)%*%x)%*%t(x)%*%y
  betahat=beta0
  #remove 0 in betahat
  {
    if(min(abs(betahat))<1e-5){
      r=which(abs(betahat)<1e-5)
      z=x[,-r]
      betahat=betahat[-r]
    }
    else {
      z=x
      betahat=betahat
    }
  }
  iter=0
  while((max(abs(t(z)%*%z%*%betahat-t(z)%*%y+n*p_lambda_d(abs(betahat),a,lambda)*sign(betahat)))>0.01)&&(iter<=1000)){
    betahat=solve(t(z)%*%z+n*Sigmaf(betahat,a,lambda))%*%t(z)%*%y
    if(min(abs(betahat))<1e-5){
      r=which(abs(betahat)<1e-5)
      z=z[,-r]
      betahat=betahat[-r]
    }
    iter=iter+1
  }
  betahate=rep(0,d)
  for (i in 1:d){
    for (j in 1:ncol(z)){
      if (length(which(x[,i]==z[,j]))==n) betahate[i]=betahat[j]
    }
  }
  px=z%*%solve(t(z)%*%z+n*Sigmaf(betahat,a,lambda))%*%t(z)
  gcv=sum((y-x%*%betahate)^2)/(n*(1-sum(diag(px))/n)^2)
  result=list(betahat=betahate,iter=iter,gcv=gcv)
  return(result)
}


#' @title A function to realize SCAD.(Jianqing Fan and Runze Li,2001)
#' @description Use scad method to solve parameter dimensionality reduction estimation in linear model.
#' @param x This is an n*d matrix which is also known as design matrix.
#' @param y This is an n*1 vector which is also known as response matrix
#' @return an n*1 vector of the estimated parameter beta
#' @import MASS
#' @examples
#' \dontrun{
#' data <- simulation1(n=40,sigma=1,beta=c(3, 1.5, 2, 0, 0, 0, 0, 0))
#' scadest(data$x, data$y)
#' }
#' @export
scadest <- function(x, y) {
  as = seq(3.5, 4, 0.1)
  lambdas = seq(0, 1, 0.1)
  gcv1 = matrix(0, nrow = length(as), ncol = length(lambdas))
  
  for (i in 1:length(as)) {
    for (j in 1:length(lambdas))
      gcv1[i, j] = scad(x, y, as[i], lambdas[j])$gcv
  }
  b = which(gcv1 == gcv1[which.min(gcv1)], arr.ind = T)
  agcv = as[b[1]]
  lambdagcv = lambdas[b[2]]
  beta1 = as.matrix(scad(x, y, agcv, lambdagcv)$betahat)
  print(list(beta1))
}






  
  

  
  