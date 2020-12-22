## ----echo=FALSE---------------------------------------------------------------
knitr::kable(summary(iris))

## ----echo=FALSE, warning=FALSE,fig.width=6------------------------------------
iris <- datasets::iris
iris2 <- iris[, -5]
species_labels <- iris[, 5]
library(colorspace) # get nice colors
species_col <- rev(rainbow_hcl(3))[as.numeric(species_labels)]
# Plot a SPLOM:
pairs(
  iris2,
  col = species_col,
  lower.panel = NULL,
  cex.labels = 2,
  pch = 19,
  cex = 1.2
)

## -----------------------------------------------------------------------------
n <- 10000
u <- runif(n)
x <- 2 / ((1 - u) ^ (1 / 2))
#After many attempts, discard some values with low probability of occurrence, and select appropriate coordinates for drawing
x <- x[x <= 10]
y <- seq(1, 10, 0.5)
hist(
x,
prob = TRUE,
breaks = y,
main = "the density histogram of the sample",
xlab = "random variates"
)
y <- seq(2, 10, 0.1)
lines(y, 8 / (y ^ 3), col = "red", lwd = 2)


## -----------------------------------------------------------------------------
n <- 10000
u1 <- runif(n, min = -1, max = 1)
u2 <- runif(n, min = -1, max = 1)
u3 <- runif(n, min = -1, max = 1)
u4 <- u3
u4 <- u2[(abs(u3) > abs(u2)) & (abs(u3) > abs(u1))]
hist(
  u4,
  prob = TRUE,
  breaks = seq(-1, 1, 0.1),
  main = "the density histogram of the sample",
  xlab = "random variates"
)
y <- seq(-1, 1, 0.025)
lines(y, 0.75 * (1 - y ^ 2), col = "red", lwd = 2)

## -----------------------------------------------------------------------------
n <- 1000
r <- 4
beta <- 2
lambda <- rgamma(n, r, beta)
x <- rexp(n, lambda)
#After many attempts, discard some values with low probability of occurrence, and select appropriate coordinates for drawing
x <- x[x <= 6]
hist(
  x,
  prob = TRUE,
  breaks = seq(0, 6, 0.5),
  main = "the density histogram of the sample",
  xlab = "random variates"
)
y <- seq(0, 6, 0.025)
lines(y, 64 / (2 + y) ^ 5, col = "red", lwd = 2)

## -----------------------------------------------------------------------------
Monte_carlo <- function(n) {
  set.seed(12306)
  # generate n random number from distribution U[0,pi/3]
  r <- runif(n, min = 0, max = pi / 3)
  # use the sample mean to estimate the population mean
  MC_value <- pi / 3 * mean(sin(r))
  # calcullate the true value use function: integrate
  true_value <- 0 + integrate(sin, 0, pi / 3)$value
  return(c(MC_value, true_value))
}
# simulation
#Use x to represent our result
x <- Monte_carlo(1000000)
names(x) <- c("MC-estimate", "true-value")
knitr::kable(x)

## -----------------------------------------------------------------------------
set.seed(12306)
MC.Phi <- function(n, R = 10000, antithetic = TRUE) {
t <- NULL 
for (i in 1:n) {
  u <- runif(R / 2)
  if (!antithetic)
    v <- runif(R / 2)
  else
    v <- 1 - u
  u <- c(u, v)
  t <- c(t, mean(exp(u)))
}
  t
}
# simulation
n <- 100
MC1 <- MC.Phi(n, anti = FALSE)
MC2 <- MC.Phi(n, anti = TRUE)
x <- c(mean(MC1), mean(MC2), exp(1) - 1)
names(x) <-
  c("simple Monte Carlo method",
    "antithetic variate approach",
    "true-value")
knitr::kable(x)

## -----------------------------------------------------------------------------
print((var(MC1) - var(MC2))/var(MC1))

## -----------------------------------------------------------------------------
set.seed(12306)
x <- seq(1,5,length.out = 20)
    g <- exp(-x^2/2)*x^2/sqrt(2*pi)
    f1 <- exp(-x+1)
    f2 <- 1/2*(x-1)^2 *exp(-x+1)
#figure (a)
plot(f1~x,type = "l",col=2,lty=2,main="figure (a)")
lines(g~x,col=1,lty=1)
lines(f2~x,col=3,lty=3)
legend("topright", legend =c("g", "f1", "f2"),
           lty = 1:3, lwd = 2, inset = 0.02,col=1:3)
#figure (b)
plot(g/f2~x,type = "l",col=2,main="figure (b)")
lines(g/f1~x,col=1)
legend("topright", legend =c("g/f1", "g/f2"),
           lty = 1:2, lwd = 2, inset = 0.02,col=1:2)


m <- 10000
  theta.hat <- se <- numeric(2)
  g <- function(x) exp(-x^2/2)*x^2/sqrt(2*pi) * (x > 1)
x <- rexp(m, rate= 1)+1 #using f1
  fg <- g(x)/exp(-x+1)
  theta.hat[1] <- mean(fg)
  se[1] <- sd(fg)
x <- rgamma(m, shape=3, rate = 1)+1 #using f2
  fg <- g(x)/(1/2*(x-1)^2 *exp(-x+1))
  theta.hat[2] <- mean(fg)
  se[2] <- sd(fg)
  res <- rbind(theta=round(theta.hat,3), se=round(se,3))
  colnames(res) <- paste0('f',1:2)
  knitr::kable(res,align='c')

## -----------------------------------------------------------------------------
F<-function(x) return((1-exp(-x))/(1-exp(-1))) # distribution function
F_inverse<-function(x) return(-log(1-(1-exp(-1))*x)) # the inverse function of distribution function
G<-function(x) return((1-exp(-1))/(5*(1+x^2))) # the function of the random variable, we want to calculate its expectation
f<-function(x) return(5*exp(-x)/(1-exp(-1))) # density function of random variable in each subinterval

quant<-F_inverse(seq(0,1,by=1/5)) # interval endpoints of the 5 subintervals
g<-(quant[-1]-quant[-6])^-1 # the density function of uniform distridution in each subinterval
theta_hat<-numeric(100)
for(k in 1:100){
theta<-numeric(5)
for(i in 1:5){
  # use acceptance-rejection method to generate the random number
 
  random_vector<-numeric(0)
  
  c<-f(quant[i])/g[i]
  while (length(random_vector)<2000) {
    Y<-runif(1,min=quant[i],max=quant[i+1])
    U<-runif(1)
    if(U<(f(Y)/(c*g[i]))) random_vector<-c(random_vector,Y)
  }
  
  theta[i]<-mean(G(random_vector))
  
}


theta_hat[k]<-sum(theta)
}
list(theta_hat=theta_hat,sd=sd(theta_hat))

## -----------------------------------------------------------------------------
cpt<-function(m,n=20,alpha=0.05){
  # m: the number of random experiment
  # n: the number of random numbers generated in each experiment
  # alpha: the significance level.
  # generate matrix of random numbers of dim m*n
  chisq<-matrix(rchisq(m*n,df=2),nrow = m,ncol = n)
  # calculate the symmetric t confidence interval according to each row of random numbers
  interval_cal<-function(x) return(c(mean(x)-var(x)/sqrt(n)*qt(1-alpha/2,df=n-1),
                                 mean(x)-var(x)/sqrt(n)*qt(alpha/2,df=n-1)))
  interval_matix<-t(apply(chisq,1,interval_cal))
  # calculate the empirical confidence level
  return(1/m*sum(2> interval_matix[,1]&2<interval_matix[,2]))
  
}
cpt(m=1000)

## -----------------------------------------------------------------------------
upper_bound<-replicate(1000,expr={
  n<-20
  alpha<-0.05
  x <- rchisq(n, df=2)
  (n-1) * var(x) / qchisq(alpha, df = n-1)
})
cpt_variance<-mean(upper_bound>4)

cat('the empirical confidence level for variance is:', cpt_variance)

## ----warning=FALSE------------------------------------------------------------
set.seed(12306)
library(knitr)
sig <-0.05 # significance level
m <-1000 # times of simulations
n <-500 # number of replications in each simulation
cv <- qnorm(1-sig/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))

# sk is the fuction used to compute the sample skewness coeff.

sk <- function(x) {
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}

alpha<- seq(0.5,10,by=0.5) # parameter of symmetric beta distribution
power<-numeric(length(alpha))

for(i in 1:length(alpha)){
   reject_i<-replicate(m,expr={
   rand_num<-rbeta(n,alpha[i],alpha[i])  
   skew<-sk(rand_num) 
   as.integer(abs(skew)>=cv) 
     
   })
   
  power[i]<-mean(reject_i)
}
plot(alpha ,power,type='l', xlab='alpha',ylab='power', main='power versus alpha for Beta(alpha,alpha)')


## ----results='asis'-----------------------------------------------------------
knitr::kable (rbind(alpha,power),format = 'html',row.names = T,digits = 2)

## -----------------------------------------------------------------------------
# the degree of freedom of the T distribution
v<-1:100

m <-1000 # times of simulations
n <-500 # number of replications in each simulation
cv <- qnorm(1-sig/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3)))) # the critical value

# store the power 
power_t<-numeric(100)

for(j in v){

reject.t<- vector('numeric',length=m)

for(i in 1:m){
  rand_num<-rt(n,df=j)
  skew<-sk(rand_num)
  reject.t[i]<-as.integer(abs(skew)>=cv)
}
power_t[j]<-mean(reject.t)
}

plot(v,power_t,xlab='df',ylab='power',type='l',main=' power versus v for t(v)')

abline(h=0.1)

## -----------------------------------------------------------------------------
set.seed(12306)
a <- 0.055
n <- c(10,30,50,100,200,500)
mu1 <- mu2 <- 0
sigma1 <- 1
sigma2 <- 1.5
m <- 1e4
result <- matrix(0, length(n), 2)

countFtest <- function(x, y){
  X <- x - mean(x)
  Y <- y - mean(y)
  x1 <- sum(X > max(Y)) + sum(X < min(Y))
  y1 <- sum(Y > max(X)) + sum(Y < min(X))
return (as.integer(max(c(x1, y1)) > 5))
}

for (i in 1:length(n)){ 
  ni <- n[i]
  tests <- replicate(m, expr={
  x <- rnorm(ni, mu1, sigma1)
  y <- rnorm(ni, mu2, sigma2)
  Fp <- var.test(x, y)$p.value
  Ftest <- as.integer(Fp <= a)
c(countFtest(x, y), Ftest)
})
result[i, ] <- rowMeans(tests)
}
data.frame(n=n, CF=result[, 1], Fp=result[, 2])

## ----eval=FALSE---------------------------------------------------------------
#  library(MASS)
#  #computes the sample skewness sttistic
#  mskr <- function(X){
#  n <- nrow(X)
#  xbar <- colMeans(X)
#  sigma.hat <- cov(X) * (n - 1) / n
#  b <- sum(((t(t(X) - xbar))%*%solve(sigma.hat)%*%(t(X) - xbar))^3) / n^2
#  return (b)
#  }
#  
#  
#  alpha=0.05 # the significance level
#  d=2 # the dimension of multivariate normal distribution
#  sigma=diag(10,d) # the covariance of multivariate normal distribution
#  n=c(10,20,30,50,100,500) # sample sizes
#  cv=qchisq(1-alpha,d*(d+1)*(d+2)/6) # critical value for the skewness test
#  p.reject=numeric(length(n)) # store sim. results
#  m=1000
#  for (i in 1:length(n)) {
#    sktests=numeric(m) # test decisions
#    for (j in 1:m) {
#      x=mvrnorm(n[i],rep(0,d),sigma)
#      sktests[j]=as.integer(abs(n[i]*mskr(x)/6) >= cv )
#    }
#    p.reject[i]=mean(sktests) #proportion rejected
#  }
#  data.frame(n=n, p.reject=p.reject)

## ----eval=FALSE---------------------------------------------------------------
#  # repeat Example 6.10 for  multivariate skewness test.
#  alpha=0.1 # the significance level
#  n=30 # the size of sample
#  m=1000 # the number of replicates
#  d=2 # the dimension of multivariate normal distribution
#  sigma1=diag(d)
#  sigma2=100*diag(d)
#  epsilon=c(seq(0,0.1,0.025),seq(0.15,0.4,0.05),seq(0.55,1,0.15))
#  N=length(epsilon)
#  pwr=numeric(N)
#  cv=qchisq(1-alpha,d*(d+1)*(d+2)/6) # critical value for the skewness test
#  for(i in 1:N){
#  e=epsilon[i]
#  sktests <- numeric(m)
#  for (j in 1:m) {
#  x=matrix(0,n,d)
#  for(k in 1:n) {if(runif(1)<=1-e) x[k,]=mvrnorm(1,rep(0,2),sigma1)
#  else x[k,]=mvrnorm(1,rep(0,2),sigma2) }
#  sktests[j] <- as.integer(n*abs(mskr(x))/6>= cv)
#  }
#  pwr[i]=mean(sktests)
#  }
#  #plot power vs epsilon
#  plot(epsilon, pwr, type = "b",xlab = bquote(epsilon),ylim = c(0,1))
#  abline(h=0.1, lty = 3)
#  se=sqrt(pwr*(1-pwr)/m) #add standard errors
#  lines(epsilon, pwr+se, lty = 3)
#  lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
set.seed(12306)
LSAT <- c(576,635,558,578,666,580,555,661,651,605,653,575,545,572,594)
GPA <- c(339,330,281,303,344,307,300,343,336,313,312,274,276,288,296)
x <- cbind(LSAT,GPA)
n <- 15
b.cor <- function(x,i) cor(x[i,1],x[i,2])
theta.hat <- b.cor(x,1:n)
theta.jack <- numeric(n)
for(i in 1:n)theta.jack[i] <- b.cor(x,(1:n)[-i])
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <-sqrt((n-1) *mean((theta.jack - mean(theta.jack))^2))
list(bias.jack=bias.jack, se.jack=se.jack)

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12306)
#  library(boot)
#  x <- c(3,5,7,18,43,85,91,98,100,130,230,487)
#  qqnorm(x);qqline(x)
#  boot.mean <- function(x,i) mean(x[i])
#  de <- boot(data=x,statistic=boot.mean,R=2000)
#  ci <- boot.ci(de,type=c("norm","basic","perc","bca"))
#  ci

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12306)
#  library(bootstrap)
#  n <- 88
#  data <- scor
#  lambda <- eigen(var(data))$values
#  theta.hat <- lambda[1]/sum(lambda)
#  
#  theta.jack <- numeric(n)
#  m <- matrix(0,5,5)
#  for (i in 1:n){
#    lambda <- eigen(var(data[-i,]))$values
#    theta.jack[i] <- lambda[1]/sum(lambda)
#  }
#  
#  
#  bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
#  se.jack <-sqrt((n-1) *mean((theta.jack - mean(theta.jack))^2))
#  list(bias.jack=bias.jack, se.jack=se.jack)

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12306)
#  
#  # calculate maximum number of extreme points for pair x,y
#  maxout <- function(x, y) {
#    X <- x - mean(x)
#    Y <- y - mean(y)
#    outx <- sum(X > max(Y)) + sum(X < min(Y))
#    outy <- sum(Y > max(X)) + sum(Y < min(X))
#    return(max(c(outx, outy)))
#  }
#  
#  # the statistics passed to boot
#  stat<-function(z,ix,n){
#    x<-z[ix][1:n]
#    y<-z[ix][-(1:n)]
#    maxout(x,y)
#  }
#  
#  # this function is used to calculate p value
#  permu_count5<-function(n1,n2,mu=0,sd1,sd2){
#    x<-rnorm(n1,mu,sd1)
#    y<-rnorm(n2,mu,sd2)
#    z<-c(x,y)
#    R=999
#    boot_obj<-boot(z,statistic = stat,R=R,sim='permutation',n=n1)
#    count<-c( boot_obj$t0, boot_obj$t)
#    p.value<-mean(count>=count[1])
#    return(p.value)
#  }
#  n<-1000
#  p_value<-numeric(n)
#  
#  # calculate the empirical type I error rate
#  for(i in 1:n) p_value[i]<-permu_count5(n1=20,n2=30,sd1=1,sd2=1)
#  cat('the empirical type I error rate is:',mean(p_value<0.05),'\n')
#  
#  # calculate the power
#  for(i in 1:n) p_value[i]<-permu_count5(n1=20,n2=30,sd1=1,sd2=2)
#  cat('the empirical power is:',mean(p_value<0.05))

## ----include=FALSE, eval=FALSE------------------------------------------------
#  library(boot)
#  library(boot)
#  library(energy)
#  library(Ball)
#  library(RANN)
#  library(ggplot2)
#  m <- 100 #permutation samples
#  p <- 2 # dimension of data
#  n1 <- n2 <- 50 #the sample size of x and y
#  R<-999 #boot parameter
#  k<-3 #boot parameter
#  n <- n1 + n2
#  N = c(n1,n2)
#  
#  # the function of NN method
#  Tn <- function(z, ix, sizes,k){
#    n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#    if(is.vector(z)) z <- data.frame(z,0);
#    z <- z[ix, ];
#    NN <- nn2(data=z, k=k+1)
#    block1 <- NN$nn.idx[1:n1,-1]
#    block2 <- NN$nn.idx[(n1+1):n,-1]
#    i1 <- sum(block1 < n1 + .5)
#    i2 <- sum(block2 > n1+.5)
#    (i1 + i2) / (k * n)
#  }
#  
#  eqdist.nn <- function(z,sizes,k){
#    boot.obj <- boot(data=z,statistic=Tn,R=R,sim = "permutation",
#                   sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(ts>=ts[1])
#    list(statistic=ts[1],p.value=p.value)
#  }
#  p.values <- matrix(NA,m,3)

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12306)
#  sd <- 1.5
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*p),ncol=p)
#    y <- matrix(rnorm(n2*p,sd=sd),ncol=p)
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value#NN method
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#energy methods
#    p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value# ball method
#  }
#  alpha <- 0.05;
#  (pow <- colMeans(p.values<alpha))
#  
#  power <- data.frame(methods = c('NN','energy','Ball'),pow)
#  ggplot(power,aes(methods,pow))+#plot
#    geom_col(fill = 'palegreen3')+
#    coord_flip()

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12306)
#  mu <- 0.5
#  sd <- 1.5
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*p),ncol=p)
#    y <- matrix(rnorm(n2*p,mean=mu,sd=sd),ncol=p)
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value#NN method
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#energy methods
#    p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value# ball method
#  }
#  alpha <- 0.05;
#  pow <- colMeans(p.values<alpha)
#  pow
#  
#  power <- data.frame(methods = c('NN','energy','Ball'),pow)
#  ggplot(power,aes(methods,pow))+#plot
#    geom_col(fill = 'palegreen3')+
#    coord_flip()

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12306)
#  mu <- 0.5
#  sd <- 2
#  for(i in 1:m){
#    x <- matrix(rt(n1*p,df=1),ncol=p)
#    y1 = rnorm(n2*p);  y2 = rnorm(n2*p,mean=mu,sd=sd)
#    w = rbinom(n, 1, .5) # 50:50 random choice
#    y <- matrix(w*y1 + (1-w)*y2,ncol=p)# normal mixture
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value#NN method
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#energy methods
#    p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value# ball method
#  }
#  alpha <- 0.05;
#  pow <- colMeans(p.values<alpha)
#  pow
#  
#  power <- data.frame(methods = c('NN','energy','Ball'),pow)
#  ggplot(power,aes(methods,pow))+
#    geom_col(fill = 'palegreen3')+
#    coord_flip()

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12306)
#  for(i in 1:m){
#    x <- matrix(rt(n1*p,df=1),ncol=p)
#    y <- matrix(rnorm(n2*p,sd=1.5),ncol=p)
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value#NN method
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#energy methods
#    p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value# ball method
#  }
#  alpha <- 0.05;
#  pow <- colMeans(p.values<alpha)
#  pow
#  
#  power <- data.frame(methods = c('NN','energy','Ball'),pow)
#  ggplot(power,aes(methods,pow))+
#    geom_col(fill = 'palegreen3')+
#    coord_flip()

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12306)
#  mu <- 0.5
#  sd <- 2
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*p),ncol=p)
#    y1 = rnorm(n2*p);  y2 = rnorm(n2*p,mean=mu,sd=sd)
#    w = rbinom(n, 1, .5) # 50:50 random choice
#    y <- matrix(w*y1 + (1-w)*y2,ncol=p)# normal mixture
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value#NN method
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#energy methods
#    p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value# ball method
#  }
#  alpha <- 0.05;
#  pow <- colMeans(p.values<alpha)
#  pow
#  
#  power <- data.frame(methods = c('NN','energy','Ball'),pow)
#  ggplot(power,aes(methods,pow))+#plot
#    geom_col(fill = 'palegreen3')+
#    coord_flip()

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12306)
#  N = c(n1,n2*10)
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*p),ncol=p);
#    y <- cbind(rnorm(n2*10),rnorm(n2*10,mean = 0.5));
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value#NN method
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#energy methods
#    p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value# ball method
#  }
#  alpha <- 0.05;
#  pow <- colMeans(p.values<alpha)
#  pow
#  
#  power <- data.frame(methods = c('NN','energy','Ball'),pow)
#  ggplot(power,aes(methods,pow))+
#    geom_col(fill = 'palegreen3')+
#    coord_flip()

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12306)
#  # pdf of standard Laplace distribution
#  laplace <- function(x) return(0.5*exp(-abs(x)))
#  
#  rw.Metropolis <- function(sigma, x0, N) {
#  # N is the number of iterations
#  x <- numeric(N)
#  x[1] <- x0 # x0 is the initial value
#  u <- runif(N) # u determines whether accept Y as x(t+1) or not
#  k <- 0 # k denotes the times of rejection
#  
#  for (i in 2:N) {
#    # the candidate is from x[i-1] plus a normal increment ~ N(0,sigma)
#    y <- rnorm(1, x[i - 1], sigma)
#    if (u[i] <= (laplace(y) / laplace(x[i - 1])))
#    x[i] <- y
#    else {
#    x[i] <- x[i - 1]
#    k <- k + 1
#    }
#  }
#  return(list(x = x, k = k))
#  }
#  
#  sigma <- c(.05, 0.5, 2, 16);N <- 2000;x0 <- 25
#  rw1 <- rw.Metropolis( sigma[1], x0, N)
#  rw2 <- rw.Metropolis( sigma[2], x0, N)
#  rw3 <- rw.Metropolis( sigma[3], x0, N)
#  rw4 <- rw.Metropolis( sigma[4], x0, N)
#  
#  plot(1:2000,rw1$x,type='l',ylab="x",xlab='iteration',main='sd=0.05')
#  abline(h=c(-3*sqrt(2),3*sqrt(2)))
#  plot(1:2000,rw2$x,type='l',ylab="x",xlab='iteration',main='sd=0.5')
#  abline(h=c(-3*sqrt(2),3*sqrt(2)))
#  plot(1:2000,rw3$x,type='l',ylab="x",xlab='iteration',main='sd=2')
#  abline(h=c(-3*sqrt(2),3*sqrt(2)))
#  plot(1:2000,rw4$x,type='l',ylab="x",xlab='iteration',main='sd=16')
#  abline(h=c(-3*sqrt(2),3*sqrt(2)))

## ----eval=FALSE---------------------------------------------------------------
#  accept_rate<-c(1-rw1$k/(N-1), 1-rw2$k/(N-1), 1-rw3$k/(N-1), 1-rw4$k/(N-1))
#  names(accept_rate)<-c('sd=0.05','sd=0.5','sd=2','sd=16')
#  print(accept_rate)

## ----eval=FALSE---------------------------------------------------------------
#  rw3 <- rw.Metropolis(2, 0, 10000)
#  y <- rw3$x[1001:10000]
#  ##Real density curve and histogram of random number generation
#  hist(y, breaks=seq(-10,10, by=0.3),
#         freq=FALSE,main='MCMC of  standard Laplace')
#  curve(laplace, from=-10, to=10, add=TRUE,col="red", lwd=3)
#  ##qqplot
#  a=ppoints(100)
#  QR=c(log(2*a[a<=0.5]),-log(2*(1-a[a>0.5]))) #quantiles of Laplace
#  Q=quantile(y, a)
#  qqplot(QR, Q, main="qqplot",
#         xlab="Laplace Quantiles", ylab="Sample Quantiles")
#  lines(c(min(y)-1,max(y)+1),c(min(y)-1,max(y)+1),col='blue')

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12306)
#  # the function for computing the mean of each chain
#  Gelman.Rubin = function(psi) {
#    # psi[i,j] is the statistic psi(X[i,1:j])
#    # for chain in i-th row of X
#    psi = as.matrix(psi)
#    n = ncol(psi)
#    k = nrow(psi)
#    psi.means = rowMeans(psi) # row means
#    B = n * var(psi.means) # between variance est.
#    psi.w = apply(psi, 1, "var") # within variances
#    W = mean(psi.w) # within est.
#    v.hat = W * (n - 1) / n + (B / n) # upper variance est.
#    r.hat = v.hat / W # G-R statistic
#    return(r.hat)
#  }
#  
#  # implement a random walk Metropolis sampler
#  rw.Metropolis = function(sigma, x0, N) {
#    x = numeric(N)
#    x[1] = x0
#    u = runif(N)
#    for (i in 2:N) {
#      y = rnorm(1, x[i - 1], sigma)
#      if (u[i] <= (laplace(y) / laplace(x[i - 1])))
#        x[i] = y # accept y
#      else
#        x[i] = x[i - 1]
#    }
#    return(x)
#  }
#  
#  # set parameters
#  k = 4 # number of chains
#  sigma = 2 # variances
#  n = 10000 # length of each chain
#  b = 1000 # burn-in length
#  x0 <- c(-30,-10,  10, 30) # initialize x0
#  
#  # generate the chains
#  X = matrix(0, nrow = k, ncol = n)
#  for (i in 1:k)
#  X[i, ] = rw.Metropolis(sigma, x0[i], n)
#  
#  # compute diagnostic statistics
#  psi = t(apply(X, 1, cumsum))
#  for (i in 1:nrow(psi))
#  psi[i, ] = psi[i, ] / (1:ncol(psi))
#  print(Gelman.Rubin(psi))
#  
#  #plot psi for the four chains
#  for (i in 1:k)
#  plot(psi[i, (b + 1):n],
#  type = "l",
#  xlab = i,
#  ylab = bquote(psi))
#  
#  #plot the sequence of R-hat statistics
#  rhat = rep(0, n)
#  for (j in (b + 1):n)
#  rhat[j] = Gelman.Rubin(psi[, 1:j])
#  plot(rhat[(b + 1):n],
#  type = "l",
#  xlab = "",
#  ylab = "R")
#  abline(h = 1.2, lty = 2)

## ----eval=FALSE---------------------------------------------------------------
#  
#  Sk_1 <- function(a, k) {
#    q <- sqrt(a ^ 2 * (k - 1) / (k - a ^ 2))
#    return (1 - pt(q, df = k - 1))
#  }
#  Sk <- function(a, k) {
#  q <- sqrt(a ^ 2 * k / (k + 1 - a ^ 2))
#  return (1 - pt(q, df = k))
#  }
#  difSK <- function(x, k) {
#  Sk_1(x, k) - Sk(x, k)
#  }
#  kset <- c(4:25, 100, 500, 1000)
#  out <- 1:length(kset)
#  for (i in 1:length(kset)) {
#  out[i] <- uniroot(
#  difSK,
#  lower = 0 + 1e-5,
#  upper = sqrt(kset[i]) - 1e-5,
#  k = kset[i]
#  )$root
#  }
#  out

## ----eval=FALSE---------------------------------------------------------------
#  kset[abs(out-sqrt(kset)) < sqrt(kset)*0.01]

## ----eval=FALSE---------------------------------------------------------------
#  n <- 1:length(kset)
#  Kwrongnum <- n[abs(out-sqrt(kset)) < sqrt(kset)*0.01]
#  
#  #Example : k=23
#  k=23
#  xx <- seq(0.01,sqrt(k)-1e-5,length=1000)
#  y <- difSK(xx,k)
#  plot(xx,y,type="l",col="red")
#  abline(h=0, lty=1)
#  
#  #Example : k=1000
#  k=1000
#  xx <- seq(0.01,sqrt(k)-1e-5,length=1000)
#  y <- difSK(xx,k)
#  plot(xx,y,type="l",col="red")
#  abline(h=0, lty=1)
#  
#  #change upper to 3
#  
#  for (i in Kwrongnum) {
#    out[i] <- uniroot(difSK,
#    lower = 0 + 1e-5,
#    upper = 3,
#    k = kset[i])$root
#  }
#  names(out) <- kset
#  
#  out
#  

## ----warning=FALSE, eval=FALSE------------------------------------------------
#  library(nloptr)
#  # Mle function
#  eval_f0 <- function(x,x1,n.A=444,n.B=132,nOO=361,nAB=63) {
#    #x[1] mean p , x1[1] mean p0
#    #x[2] mean q , x1[2] mean q0
#    r1<-1-sum(x1)
#    nAA<-n.A*x1[1]^2/(x1[1]^2+2*x1[1]*r1)
#    nBB<-n.B*x1[2]^2/(x1[2]^2+2*x1[2]*r1)
#    r<-1-sum(x)
#    return(-2*nAA*log(x[1])-2*nBB*log(x[2])-2*nOO*log(r)-
#             (n.A-nAA)*log(2*x[1]*r)-(n.B-nBB)*log(2*x[2]*r)-nAB*log(2*x[1]*x[2]))
#  }
#  
#  
#  # constraint function
#  eval_g0 <- function(x,x1,n.A=444,n.B=132,nOO=361,nAB=63) return(sum(x)-0.999999)
#  
#  
#  opts <- list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1.0e-8)
#  mle<-NULL
#  r<-matrix(0,1,2)
#  r<-rbind(r,c(0.2,0.35))# the beginning value of p0 and q0
#  j<-2
#  while (sum(abs(r[j,]-r[j-1,]))>1e-8) {
#  res <- nloptr( x0=c(0.3,0.25),
#                 eval_f=eval_f0,
#                 lb = c(0,0), ub = c(1,1),
#                 eval_g_ineq = eval_g0,
#                 opts = opts, x1=r[j,],n.A=444,n.B=132,nOO=361,nAB=63)
#  j<-j+1
#  r<-rbind(r,res$solution)
#  mle<-c(mle,eval_f0(x=r[j,],x1=r[j-1,]))
#  }
#  r  #the result of EM algorithm
#  list(mle=-mle) #the max likelihood values
#  

## ----eval=FALSE---------------------------------------------------------------
#  formulas<-list(
#    mpg~disp,
#    mpg~I(1/disp),
#    mpg~disp+wt,
#    mpg~I(1/disp)+wt
#  )
#  
#  # lapply
#  fit1<-lapply(1:4,function(i) lm(formula = formulas[[i]],data=mtcars))
#  
#  # forloop
#  
#  fit2<-vector('list',length = 4)
#  for(i in seq_along(formulas)) fit2[[i]]<-lm(formula = formulas[[i]],data=mtcars)
#  
#  fit1
#  fit2

## ----eval=FALSE---------------------------------------------------------------
#  trials <- replicate(100,t.test(rpois(10, 10), rpois(7, 10)),simplify = FALSE)
#  
#  # sapply with anonymous function
#  pvalue<-sapply(trials,function(test) test$p.value)
#  
#  # sapply without anonymous function
#  pvalue1<-sapply(trials,'[[','p.value')
#  
#  pvalue
#  pvalue1

## ----eval=FALSE---------------------------------------------------------------
#  testlist <- list(mtcars, cars)
#  lapply(testlist, function(x) vapply(x, mean, numeric(1)))

## ----eval=FALSE---------------------------------------------------------------
#  lmapply <- function(X, f, f.value, simplify = FALSE) {
#    out <- Map(function(x)
#    vapply(x, f, f.value), X)
#    if (simplify == TRUE) {
#    return(simplify2array(out))
#    }
#    unlist(out, recursive = FALSE)
#  }
#  lmapply(testlist, mean, numeric(1))
#  

## ----warning=FALSE, eval=FALSE------------------------------------------------
#  library(Rcpp)
#  ## 1.  R random number generater
#  
#  # pdf of standard Laplace distribution
#  laplace<-function(x) return(1/2*exp(-abs(x)))
#  
#  rw.Metropolis <- function(sigma, x0, N) {
#  
#  # N is the number of iterations
#  x <- numeric(N)
#  # x0 is the initial value
#  x[1] <- x0
#  
#  # u determines whether accept Y as x(t+1) or not
#  u <- runif(N)
#  
#  # k denotes the times of rejection
#  k <- 0
#  
#  for (i in 2:N) {
#    # the candidate is from x[i-1] plus a normal increment ~ N(0,sigma)
#    y <- rnorm(1, x[i - 1], sigma)
#    if (u[i] <= (laplace(y) / laplace(x[i - 1])))
#    x[i] <- y
#    else {
#    x[i] <- x[i - 1]
#    k <- k + 1
#    }
#  }
#  return(list(x = x, k = k))
#  }
#  
#  
#  ## 2. C++ random number generater: function(Metropolis)
#  library(StatComp20064)
#  
#  
#  sigma <- c(.05, .5, 2, 16);N <- 2000;x0 <- 25
#  
#  rw1 <- rw.Metropolis( sigma[1], x0, N)
#  rw2 <- rw.Metropolis( sigma[2], x0, N)
#  rw3 <- rw.Metropolis( sigma[3], x0, N)
#  rw4 <- rw.Metropolis( sigma[4], x0, N)
#  
#  cpp.rw1<-Metropolis( sigma[1], x0, N)
#  cpp.rw2<-Metropolis( sigma[2], x0, N)
#  cpp.rw3<-Metropolis( sigma[3], x0, N)
#  cpp.rw4<-Metropolis( sigma[4], x0, N)
#  
#  
#  plot(1:2000,rw1$x,type='l',ylab="x",xlab='iteration',main='sd=0.05(R)')
#  
#  plot(1:2000,rw2$x,type='l',ylab="x",xlab='iteration',main='sd=0.5(R)')
#  abline(h=c(-3*sqrt(2),3*sqrt(2)))
#  plot(1:2000,rw3$x,type='l',ylab="x",xlab='iteration',main='sd=2(R)')
#  abline(h=c(-3*sqrt(2),3*sqrt(2)))
#  plot(1:2000,rw4$x,type='l',ylab="x",xlab='iteration',main='sd=16(R)')
#  abline(h=c(-3*sqrt(2),3*sqrt(2)))
#  
#  
#  plot(1:2000,cpp.rw1[,1],type='l',ylab="x",xlab='iteration',main='sd=0.05(Cpp)')
#  plot(1:2000,cpp.rw2[,1],type='l',ylab="x",xlab='iteration',main='sd=0.5(Cpp)')
#  abline(h=c(-3*sqrt(2),3*sqrt(2)))
#  plot(1:2000,cpp.rw3[,1],type='l',ylab="x",xlab='iteration',main='sd=2(Cpp)')
#  abline(h=c(-3*sqrt(2),3*sqrt(2)))
#  plot(1:2000,cpp.rw4[,1],type='l',ylab="x",xlab='iteration',main='sd=16(Cpp)')
#  abline(h=c(-3*sqrt(2),3*sqrt(2)))

## ----eval=FALSE---------------------------------------------------------------
#  
#  qqplot(rw1$x[500:2000],cpp.rw1[500:2000,1],xlab='R',ylab='cpp',main='sd=0.05')
#  qqplot(rw2$x[500:2000],cpp.rw2[500:2000,1],xlab='R',ylab='cpp',main='sd=0.5')
#  qqplot(rw3$x[500:2000],cpp.rw3[500:2000,1],xlab='R',ylab='cpp',main='sd=2')
#  qqplot(rw4$x[500:2000],cpp.rw4[500:2000,1],xlab='R',ylab='cpp',main='sd=16')

## ----eval=FALSE---------------------------------------------------------------
#  library(microbenchmark)
#  n <- 2000
#  ts1 <-
#  microbenchmark(R = rw.Metropolis(0.05, 25, n), cpp = Metropolis(0.05, 25, n))
#  
#  ts2 <-
#  microbenchmark(R = rw.Metropolis(0.5, 25, n), cpp = Metropolis(0.5, 25, n))
#  
#  ts3 <-
#  microbenchmark(R = rw.Metropolis(2, 25, n), cpp = Metropolis(2, 25, n))
#  
#  ts4 <-
#  microbenchmark(R = rw.Metropolis(16, 25, n), cpp = Metropolis(16, 25, n))
#  
#  summary(ts1)[, c(1, 3, 5, 6)]
#  summary(ts2)[, c(1, 3, 5, 6)]
#  summary(ts3)[, c(1, 3, 5, 6)]
#  summary(ts4)[, c(1, 3, 5, 6)]

