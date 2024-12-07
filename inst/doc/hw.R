## -----------------------------------------------------------------------------
set.seed(12345)
n<-1000
sigma<-1 #set the value of sigma
x<-rnorm(n,sd=sigma)
y<-rnorm(n,sd=sigma)
r<-sqrt(x^2+y^2)
h<-rweibull(n,shape=2)
par(mfrow=c(1,2))
hist(r) #algorithm to generate random samples
hist(h) #Rayleigh distribution

## -----------------------------------------------------------------------------
sigma<-2 #set the value of sigma
x<-rnorm(n,sd=sigma)
y<-rnorm(n,sd=sigma)
r<-sqrt(x^2+y^2)
h<-rweibull(n,shape=2,scale=2)
par(mfrow=c(1,2))
hist(r) #algorithm to generate random samples
hist(h) #Rayleigh distribution

## -----------------------------------------------------------------------------
sigma<-3 #set the value of sigma
x<-rnorm(n,sd=sigma)
y<-rnorm(n,sd=sigma)
r<-sqrt(x^2+y^2)
h<-rweibull(n,shape=2,scale=3)
par(mfrow=c(1,2))
hist(r) #algorithm to generate random samples
hist(h) #Rayleigh distribution

## -----------------------------------------------------------------------------
set.seed(12345)
n=1000
p1=0.75
p2=1-p1
x=rnorm(n)
y=rnorm(n,mean=3)
r=sample(c(0,1),n,replace=TRUE,prob=c(p2,p1))
z=r*x+(1-r)*y
hist(z)

## -----------------------------------------------------------------------------
p1=0.5
p2=1-p1
r=sample(c(0,1),n,replace=TRUE,prob=c(p2,p1))
z=r*x+(1-r)*y
hist(z)

## -----------------------------------------------------------------------------
p1=0.25
p2=1-p1
r=sample(c(0,1),n,replace=TRUE,prob=c(p2,p1))
z=r*x+(1-r)*y
hist(z)

## -----------------------------------------------------------------------------
r=1
beta=1
n=1000
t=10
nt=rpois(n,t)
xt=numeric(n)
i=j=0
for(i in 1:n){
  for(j in 1:nt[i]){
    xt[i]=xt[i]+rgamma(1,r,beta)
  }
}
hist(xt)  
m1=mean(xt)
v1=var(xt)
cat("the estimation of mean is ",m1,",the variance is ",v1,".the theoretical value of mean is 10,the variance is 20.")

## -----------------------------------------------------------------------------
r=2
beta=2
n=1000
t=10
nt=rpois(n,t)
xt=numeric(n)
i=j=0
for(i in 1:n){
  for(j in 1:nt[i]){
    xt[i]=xt[i]+rgamma(1,r,beta)
  }
}
hist(xt)  
m2=mean(xt)
v2=var(xt)
cat("the estimation of mean is ",m2,",the variance is ",v2,".the theoretical value of mean is 10,the variance is 15.")

## -----------------------------------------------------------------------------
r=4
beta=2
n=1000
t=10
nt=rpois(n,t)
xt=numeric(n)
i=j=0
for(i in 1:n){
  for(j in 1:nt[i]){
    xt[i]=xt[i]+rgamma(1,r,beta)
  }
}
hist(xt)  
m3=mean(xt)
v3=var(xt)
cat("the estimation of mean is ",m3,",the variance is ",v3,".the theoretical value of mean is 20,the variance is 50.")

## -----------------------------------------------------------------------------
MC_beta_cdf<-function(x,n){
  t<-runif(n,min=0,max=x)
  theta_hat<-mean(30*t^2*(1-t)^2)*x
  return(theta_hat)
}
cb<-numeric(0)
for(i in 1:9){
  cb[i]<-MC_beta_cdf(i/10,10000)
}
par(mfrow=c(1,2))
plot(cb)
xb1<-pbeta(seq(0,1,by=0.1),shape1=3,shape2=3)
plot(xb1)

## -----------------------------------------------------------------------------
MC.Phi <- function(x,sigma, R = 10000, antithetic = FALSE) {
  u <- runif(R/2)
  if (antithetic) v <- 1 - u else v <- runif(R/2)
  u <- c(u, v)
  g <- x*(u*x)*exp(-(u*x)^2/(2*sigma^2)) # x*u ~ N(0,x)
  cdf <- mean(g) / (sigma^2)
  cdf
}
m <- 1000
MC1 <- MC2 <- numeric(m)
x <- 1
for (i in 1:m) {
  MC1[i] <- MC.Phi(x,1, R = 1000, antithetic = FALSE)
  MC2[i] <- MC.Phi(x,1, R = 1000, antithetic = TRUE)
}
round(c(sd(MC1),sd(MC2),sd(MC2)/sd(MC1)),5)

## -----------------------------------------------------------------------------
n=10000
x1<-rnorm(n)
g<-function(x){
  x^2*exp(-x^2/2)/sqrt(2*pi)*(x>=1)
}
fg1<-g(x1)/dnorm(x1)
theta1<-mean(fg1)
sd1<-sd(fg1)

## -----------------------------------------------------------------------------
x2<-rexp(n)
fg2<-g(x2)/exp(-x2)
theta2<-mean(fg2)
sd2<-sd(fg2)

theta<-c(theta1,theta2)
sd<-c(sd1,sd2)
res<-rbind(theta=round(theta,3),sd=round(sd,3))
res

## -----------------------------------------------------------------------------
fsa<-function(x){  #fast sorting algorithm function
  if(length(x)<=1) {
    return(x)
  }
  p<-x[1]
  l<-x[x<p]
  r<-x[x>p]
  return(c(fsa(l),p,c(fsa(r))))
}
x1<-sample(1e4)
x2<-sample(2e4)
x3<-sample(4e4)
x4<-sample(6e4)
x5<-sample(8e4)

t1<-t2<-t3<-t4<-t5<-numeric(100)
for(i in 1:100){
  t1[i]<-system.time({fsa(x1)})[1]
  t2[i]<-system.time({fsa(x2)})[1]
  t3[i]<-system.time({fsa(x3)})[1]
  t4[i]<-system.time({fsa(x4)})[1]
  t5[i]<-system.time({fsa(x5)})[1]
}
t<-c(mean(t1),mean(t2),mean(t3),mean(t4),mean(t5))

## -----------------------------------------------------------------------------
n<-c(1e4,2e4,4e4,6e4,8e4)
x0<-n*log(n)
model<-lm(t~x0)
r1<-model$coefficients[1]
r2<-model$coefficients[2]

## -----------------------------------------------------------------------------
plot(n,t)
lines(n,r1+r2*n*log(n))

## -----------------------------------------------------------------------------
skewness<-function(x){
  #compute the sample skewness statistics
  xbar<-mean(x)
  m1<-mean((x-xbar)^3)
  m2<-mean((x-xbar)^2)
  m3<-m1/m2^1.5
  return(m3)
}

eq<-function(x,q){
  #estimate the quantiles
  t<-quantile(x,q)
  return(t)
}

n=1e4
m=1e2
sk<-numeric(n)
for(i in 1:n){
  sk[i]<-skewness(rnorm(m))
}

## -----------------------------------------------------------------------------
q<-eq(sk,c(0.025,0.05,0.95,0.975))
varq<-function(q){
  #compute the standard error
  f<-exp(-q^2/2)/sqrt(2*pi)
  r<-q*(1-q)/n/f^2
  return(r)
}
varq(q)

## -----------------------------------------------------------------------------
#compare the estimated quantiles with the quantiles of the large sample approximation sqrt(b_1) approx N(0,6/n)
com<-rnorm(n,sd=sqrt(6/m))
qcom<-eq(com,c(0.025,0.05,0.95,0.975))
result1<-function(q,qcom){
  #result reporting
  print(q)
  print(qcom)
}
result1(q,qcom)

## -----------------------------------------------------------------------------
library(MASS)
set.seed(12345)
data<-mvrnorm(10000,mu=c(0,0),Sigma=matrix(c(1,0,0,1),2,2))
cor.test(data[,1],data[,2],method=c("pearson"))
cor.test(data[,1],data[,2],method=c("spearman"))
cor.test(data[,1],data[,2],method=c("kendall"))

## -----------------------------------------------------------------------------
x<-rexp(1000)
y<-rnorm(1000)
cor.test(x,y,method="pearson")
cor.test(x,y,method="spearman")
cor.test(x,y,method="kendall")

## -----------------------------------------------------------------------------
fwerbon<-fwerbh<-fdrbon<-fdrbh<-tprbon<-tprbh<-numeric(10000)
for(j in 1:10000){
  p0<-runif(950)
  p1<-rbeta(50,0.1,1)
  p<-sort(c(p0,p1))
  p.adj1<-p.adjust(p,method='bonferroni')
  p.adj2<-p.adjust(p,method='BH')
  fwer<-fdr<-tpr<-numeric(2)
  for(i in 1:1000){
    if(p.adj1[i]<0.1/1000) {fwer[1]<-fwer[1]+1}
    if(p.adj2[i]<0.1/1000) {fwer[2]<-fwer[2]+1}
    if(p.adj1[i]<0.1*i/1000) {fdr[1]<-fdr[1]+1}
    if(p.adj2[i]<0.1*i/1000) {fdr[2]<-fdr[2]+1}
  }
  fwerbon[j]<-fwer[1]/1000
  fwerbh[j]<-fwer[2]/1000
  fdrbon[j]<-fdr[1]/1000
  fdrbh[j]<-fdr[2]/1000
  tprbon[j]<-1-fdrbon[j]
  tprbh[j]<-1-fdrbh[j]
}

## -----------------------------------------------------------------------------

r<-matrix(c(mean(fwerbon),mean(fwerbh),mean(fdrbon),mean(fdrbh),mean(tprbon),mean(tprbh)),byrow=TRUE,nrow=3)
colnames(r)<-c('bonferroni','BH')
rownames(r)<-c('FWER','FDR','TPR')
round(r,4)

## -----------------------------------------------------------------------------
# compute the MLE of lambda_hat
library(boot)
f<-function(lambda){
  return(length(x)*log(lambda)-lambda*sum(x))
}
x<-c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
optimize(f,interval=c(0,2),maximum = TRUE)

## -----------------------------------------------------------------------------
#Bootstrap method
B=1e4
set.seed(12345)
lambdastar<-numeric(B)
lambda1<-1/mean(x)
for(i in 1:B){
  xstar<-sample(x,replace=TRUE)
  lambdastar[i]<-1/mean(xstar)
}
round(c(bias=mean(lambdastar)-lambda1,se.boot=sd(lambdastar)),5)

## -----------------------------------------------------------------------------
#basic method
fb<-function(x,i){return(mean(x[i]))}
bx<-boot(data=x,statistic=fb,R=999)
ci<-boot.ci(bx,type=c("norm","basic","perc","bca"))
ci

## -----------------------------------------------------------------------------
library(bootstrap)
data1<-scor
dc<-cov(data1)
f1<-function(dc){
  #compute the theta
  l<-eigen(dc,only.values = TRUE)$values
  theta<-l[1]/sum(l)
  return(theta)
}
theta<-f1(dc)

## -----------------------------------------------------------------------------
n<-dim(data1)[1]
thetahat<-numeric(n)
f2<-function(data1){
  #jackknife
  for(i in 1:n){
    data2<-data1[-i]
    c<-cov(data2)
    thetahat[i]<-f1(c)
  }
  return(thetahat)
}
thetahat<-f2(data1)
bias.jack<-(n-1)*(mean(thetahat)-theta)
se.jack<-sqrt((n-1)*mean((thetahat-theta)^2))


f3<-function(theta,bias.jack,se.jack){
  #show the result
  round(c(original=theta,bias.jack=bias.jack,se.jack=se.jack),3)
}
f3(theta,bias.jack,se.jack)


## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
  e1[k] <- magnetic[k] - yhat1
  
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
  J2$coef[3] * chemical[k]^2
  e2[k] <- magnetic[k] - yhat2
  
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
  yhat3 <- exp(logyhat3)
  e3[k] <- magnetic[k] - yhat3
  
  J4 <- lm(y ~ x + I(x^2) + I(x^3))
  yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] +
  J4$coef[3] * chemical[k]^2 + J4$coef[4]*chemical[k]^3
  e4[k] <- magnetic[k] - yhat4
}
c(mean(e1^2),mean(e2^2),mean(e3^2),mean(e4^2))

## -----------------------------------------------------------------------------
summary(J1)
summary(J2)
summary(J3)
summary(J4)

## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)

## -----------------------------------------------------------------------------
library(twosamples) #call the function cvm_test
R <- 999 #number of replicates
z <- c(x, y) #pooled sample
K <- 1:26
options(warn = -1)
D0 <- cvm_test(x, y)[2]

cvm<-function(){
  D<-numeric(R)
  for (i in 1:R) {
    #generate indices k for the first sample
    k <- sample(K, size = 14, replace = FALSE)
    x1 <- z[k]
    y1 <- z[-k] 
    D[i] <- cvm_test(x1, y1)[2]
  }
  return(D)
}
D<-cvm()
cp<-function(D0,D){
  p <- mean(c(D0, D) >= D0)
  return(p)
}
p<-cp(D0,D)
options(warn = 0)
p

## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
attach(chickwts)
x <- sort(as.vector(weight[feed == "sunflower"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)

## -----------------------------------------------------------------------------
R<-999
z<-c(x,y)
K<-1:24
options(warn=-1)
D0<-cor(x,y,method="spearman")

spearman<-function(){
  D<-numeric(R)
  for (i in 1:R) {
    k <- sample(K, size = 12, replace = FALSE)
    x1 <- z[k]
    y1 <- z[-k] 
    D[i] <- cor(x1,y1,method="spearman")
  }
  return(D)
}
D<-spearman()

sp<-function(D0,D){
  p <- mean(c(D0, D) >= D0)
  return(p)
}
p<-sp(D0,D)
options(warn = 0)

r<-function(p){
  if(p<=0.05){print("reject the null hypothesis")}
  else {print("accept the null hypothesis")}
}
r(p)

## -----------------------------------------------------------------------------
cp<-cor.test(x,y)[3]
r(cp)

## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
fc<-function(x,eta,theta){
  return(1/(theta*pi*(1+((x-eta)/theta)^2)))
}

MH<-function(){
  m<-10000
  eta<-0
  theta<-1
  x<-numeric(m)
  x[1]<-rnorm(1,1)
  u<-runif(m)
  
  for(i in 2:m){
    xt<-x[i-1]
    y<-rnorm(1,xt)
    num<-fc(y,eta,theta)*dnorm(xt,y)
    den<-fc(xt,eta,theta)*dnorm(y,xt)
    if(u[i] <= num/den){
      x[i]<-y
    }
    else{
      x[i]<-xt
    }
  }
  return(x)
}

## -----------------------------------------------------------------------------
q<-c(10,20,30,40,50,60,70,80,90)/100
x<-MH()[1001:10000]
nc<-rcauchy(9000)
m1<-quantile(x,q,na.rm=TRUE)
m2<-quantile(nc,q)

r<-function(m1,m2){
  print("the deciles of the generated observations is ")
  print(m1)
  print("the deciles of the standard Cauchy distribution is ")
  print(m2)
}
r(m1,m2)

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}


## -----------------------------------------------------------------------------
n=9000

#generate the chains
for(k in 2:10){
  X <- matrix(0, nrow=k, ncol=n)
  for (i in 1:k){
    X[i, ] <- MH()[1001:10000]
  }
  
#Compute the diagnostic statistics
  psi <- t(apply(X, 1, cumsum))
  for (i in 1:nrow(psi)){
    psi[i,] <- psi[i,] / (1:ncol(psi))
  }
  if(Gelman.Rubin(psi)<1.2) break
}
print(c(Gelman.Rubin(psi),k))

## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
N=5000
burn=1000
a=b=1
n=100
Gibbs<-function(){
  X<-matrix(0,N,2)
  
  y<-0.5
  x<-50
  X[1,]<-c(x,y)
  
  for(i in 2:N){
    x2<-X[i-1,2]
    X[i,1]<-rbinom(1,100,x2)
    x1<-X[i-1,1]
    X[i,2]<-rbeta(1,x1+a,n-x1+b)
  }
  b<-burn+1
  x<-X[b:N,] 
  return(x)
}

fcom<-function(x,n){
  return(factorial(x)*factorial(n-x)/factorial(n))
}

fd<-function(x,y){
  #density f(x,y)
  return(fcom(x,n)*(y^(x+a-1))*((1-y)^(n-x+b-1)))
}

x<-Gibbs()
result<-fd(x[,1],x[,2])

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}

## -----------------------------------------------------------------------------
n=4000

#generate the chains
for(k in 2:10){
  X1 <- matrix(0, nrow=k, ncol=n)
  X2 <- matrix(0, nrow=k, ncol=n)
  for (i in 1:k){
    X1[i, ] <- Gibbs()[,1]
    X2[i, ] <- Gibbs()[,2]
  }
  
#Compute the diagnostic statistics
  psi1 <- t(apply(X1, 1, cumsum))
  for (i in 1:nrow(psi1)){
    psi1[i,] <- psi1[i,] / (1:ncol(psi1))
  }
  psi2 <- t(apply(X2, 1, cumsum))
  for (i in 1:nrow(psi2)){
    psi2[i,] <- psi2[i,] / (1:ncol(psi2))
  }
  if(Gelman.Rubin(psi1)<1.2&Gelman.Rubin(psi2)<1.2) break
}
print(c(Gelman.Rubin(psi1),Gelman.Rubin(psi2)<1.2,k))

## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------

an<-function(a){
  #Euclidean norm
  return(sqrt(sum(a^2)))
}

lf<-function(n){
  #log factorial function
  if(n==0){
    return(0)
  }else{
    return(log(n)+lf(n-1))
  }
}

kt<-function(a,k,d){
  #k_th term
  lterm<-log(an(a))+lgamma((d+1)/2)+lgamma(k+3/2)-lf(k)-k*log(2)-log(2*k+1)-log(2*k+2)-lgamma(k+d/2+1)
  return((-1)^k*exp(lterm))
}


## -----------------------------------------------------------------------------
fsum<-function(a,k,d){
  f<-0
  for(i in 0:k){
    f<-f+kt(a,i,d)
  }
  return(f)
}

## -----------------------------------------------------------------------------
a<-c(1,2)
k<-2000
d<-2
fsum(a,k,d)

## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------

ck<-function(a,k){
  return(sqrt((a^2)*k/(k+1-a^2)))
}

f<-function(a){
  pt(ck(a,k),df=k)-pt(ck(a,k-1),df=k-1)
}

## -----------------------------------------------------------------------------
fa<-function(k){
  res<-uniroot(f,c(0.1,sqrt(k)-0.1))
  return(res$root)
}

a<-numeric(22)
for(k in 4:25){
  a[k-3]<-fa(k)
}
names(a)<-4:25
print(a)

## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
T0<-c(0.54,0.48,0.33,0.43,1.00,1.00,0.91,1.00,0.21,0.85)

MLE<-function(T0){
  return(length(T0)/sum(T0))
}

lambda0<-lambda<-MLE(T0)
that<-1+1/lambda
T0[T0==1.00]<-that


N=10000
tol<-.Machine$double.eps^0.5
for(i in 1:N){
  t<-lambda
  lambda<-MLE(T0)
  that<-1+1/lambda
  T0[T0==t]<-that
  if(abs(lambda-t)/t<tol) break
}

round(c(MLE=lambda0,EM=lambda),4)

## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
library(boot)
A1<-rbind(c(2,1,1),c(1,-1,3))
b1<-c(1,3)
a<-c(4,2,9)
simplex(a=a,A1=A1,b1=b1,maxi=FALSE)


## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
formulas <- list(
    mpg ~ disp,
    mpg ~ I(1 / disp),
    mpg ~ disp + wt,
    mpg ~ I(1 / disp) + wt
)

mod1<-lapply(formulas,function(x){lm(x,data=mtcars)})

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})

formulas<-list(y~x)
mod2<-list()  #save the result
for(i in seq_along(bootstraps)){
  y<-bootstraps[[i]]$mpg
  x<-bootstraps[[i]]$disp
  mod2<-append(mod2,lapply(formulas,lm))
}

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

rsq1<-numeric(4)
rsq2<-numeric(10)

for(i in 1:4){
  rsq1[i]<-rsq(mod1[[i]])
}

for(i in 1:10){
  rsq2[i]<-rsq(mod2[[i]])
}

round(c(rsq1=rsq1),3)
round(c(rsq2=rsq2),3)

## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)

pv1<-sapply(trials,function(x) x$p.value)

pv2<-numeric(100)
for(i in 1:100){
  pv2[i]<-trials[[i]]$p.value
}
round(pv1,3)
round(pv2,3)

## -----------------------------------------------------------------------------
lapply1<-function(x,fun,fun.value){
  rv<-vapply(x,FUN = fun,FUN.VALUE = fun.value)
  Map(as.matrix,rv)
}

x<-data.frame(cbind(x1=c(1,2,3,4),x2=c(5,6,7,8)))
result<-lapply1(x,mean,numeric(1))
result

## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
fastchisq.test<-function(x,y){
  if(is.numeric(x)&&is.numeric(y)){
    if(any(is.na(x)|is.na(y))) print("the vectors have missing values")
    else return(chisq.test(x,y))
  }
  else print("not numeric vectors")
}

#example
x<-c(1,2,2,1)
y<-c(1,1,1,2)
fastchisq.test(x,y)

x<-c(1,2,2,1)
y<-c(1,2,2,NA)
fastchisq.test(x,y)

x<-c(1,2,2,1)
y<-c(1,2,2,"a")
fastchisq.test(x,y)

## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------

ftable<-function(x,y){
  if(!is.integer(x)||!is.integer(y)||any(is.na(x)|is.na(y))){
    print("error")
  }else{
    table(x,y)
  }
}

#example
x<-as.integer(c(1,2,3))
y<-as.integer(c(3,1,2))
ftable(x,y)

x<-c(0.1,2,3)
y<-c(3,1,0.2)
ftable(x,y)

## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
library(Rcpp)
dir_cpp<-'D:/a/'
sourceCpp(paste0(dir_cpp,"Gibbs.cpp"))
N <- 5000
burn <- 1000
a <- 1
b <- 1
n <- 100


result <- Gibbs(N, burn, n, a, b)
result<-result[1001:5000,]

fcom<-function(x,n){
  return(factorial(x)*factorial(n-x)/factorial(n))
}

fd<-function(x,y){
  #density f(x,y)
  return(fcom(x,n)*(y^(x+a-1))*((1-y)^(n-x+b-1)))
}

fcpp<-fd(result[,1],result[,2])

## -----------------------------------------------------------------------------
Gibbsr<-function(){
  X<-matrix(0,N,2)
  
  y<-0.5
  x<-50
  X[1,]<-c(x,y)
  
  for(i in 2:N){
    x2<-X[i-1,2]
    X[i,1]<-rbinom(1,100,x2)
    x1<-X[i-1,1]
    X[i,2]<-rbeta(1,x1+a,n-x1+b)
  }
  b<-burn+1
  x<-X[b:N,] 
  return(x)
}

x<-Gibbsr()
fr<-fd(x[,1],x[,2])

## -----------------------------------------------------------------------------
qqplot(log(fcpp),log(fr))

## -----------------------------------------------------------------------------
library(microbenchmark)
ts<-microbenchmark(timecpp=Gibbs(N, burn, n, a, b),timer=Gibbsr())
summary(ts)[,c(1,3,5,6)]

