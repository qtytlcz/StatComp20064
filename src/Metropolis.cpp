#include <Rcpp.h>
using namespace Rcpp;
//' @title A Metropolis sampler using Rcpp
//' @description A Metropolis sampler using Rcpp
//' @param N the number of samples
//' @param x0 the initial value
//' @param sigma the standard deviation of the normal distribution
//' @return a random sample of size \code{N} 
//' @useDynLib StatComp20064
//' @examples
//' \dontrun{
//' cpp.rw1=Metropolis(sigma=0.05,x0=25,N=2000)
//' plot(1:2000,cpp.rw1[,1],type='l',ylab="x",xlab='iteration',main='sd=0.05(Cpp)')
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix Metropolis(double sigma=0.05,double x0=25,int N=2000) {
  NumericMatrix mat(N, 2);
  mat(0,0)=x0;mat(0,1)=0;int k=0;
  for (int i=1; i<N;i++) {
    double z = runif(1,0,1)[0];
    double y = rnorm(1,mat(i-1, 0),sigma)[0];
    if (z < (exp(-abs(y))/exp(-abs(mat(i-1, 0))))) mat(i, 0) = y ;
    else { mat(i, 0) = mat(i-1, 0);++k;}
    mat(i, 1) = k;
  }
  return(mat);
}

 






