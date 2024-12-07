## ----eval=FALSE---------------------------------------------------------------
#  NR<-function(x,theta,method){
#    if(method==1){
#      return(theta+sum(2*(x-theta)/(1+(x-theta)^2))/sum(2*(1-(x-theta)^2)/(1+(x-theta)^2)^2))
#    }
#    if(method==2){
#      t<-theta+sum(2*(x-theta)/(1+(x-theta)^2))/sum(2*(1-(x-theta)^2)/(1+(x-theta)^2)^2)
#      k<-1
#      theta<-theta
#      while(abs(theta-t)>1e-3){
#        theta<-t
#        t<-theta+sum(2*(x-theta)/(1+(x-theta)^2))/sum(2*(1-(x-theta)^2)/(1+(x-theta)^2)^2)
#        k<-k+1
#      }
#      return(t)
#    }
#    if(method!=1 && method !=2){
#      print("please input the correct method")
#    }
#  }

## ----eval=FALSE---------------------------------------------------------------
#  double NC(const NumericVector &x, double theta, int method) {
#    if (method == 1) {
#      return (theta + sum(2.0 * (x - theta) / (1.0 + pow(x - theta, 2))) / sum(2.0 * (1.0 - pow(x - theta, 2)) / pow(1.0 + pow(x - theta, 2), 2)));
#    } else if (method == 2) {
#      double t = theta + sum(2.0 * (x - theta) / (1.0 + pow(x - theta, 2))) / sum(2.0 * (1.0 - pow(x - theta, 2)) / pow(1.0 + pow(x - theta, 2), 2));
#      int k = 1;
#      while (std::abs(theta - t) > 1e-3) {
#        theta = t;
#       t = theta + sum(2.0 * (x - theta) / (1.0 + pow(x - theta, 2))) / sum(2.0 * (1.0 - pow(x - theta, 2)) / pow(1.0 + pow(x - theta, 2), 2));
#       k++;
#       if (k > 1000) {
#         Rcpp::Rcout << "error" << std::endl;
#         break;
#       }
#     }
#     return t;
#    }else{
#      Rcpp::Rcout << "please input the correct method" << std::endl;
#      return 0;
#   }
#  }

## ----eval=FALSE---------------------------------------------------------------
#  FSR<-function(x,theta,method){
#    if(method==1){
#      return(theta+2*mean(2*(x-theta)/(1+(x-theta)^2)))
#    }
#    if(method==2){
#      t<-theta+2*mean(2*(x-theta)/(1+(x-theta)^2))
#      k<-1
#      theta<-theta
#      while(abs(theta-t)>1e-3){
#        theta<-t
#        t<-theta+2*mean(2*(x-theta)/(1+(x-theta)^2))
#        k<-k+1
#      }
#      return(t)
#    }
#    if(method!=1 && method !=2){
#      print("please input the correct method")
#    }
#  }

## ----eval=FALSE---------------------------------------------------------------
#  double FSC(const NumericVector &x, double theta, int method) {
#    if (method == 1) {
#      return (theta + 2.0 * mean(2.0 * (x - theta) / (1.0 + pow(x - theta, 2))));
#    } else if (method == 2) {
#      double t = theta + 2.0 * mean(2.0 * (x - theta) / (1.0 +pow(x - theta, 2)));
#      int k = 1;
#      while (std::abs(theta - t) > 1e-3) {
#        theta = t;
#        t = theta + 2.0 * mean(2.0 * (x - theta) / (1.0 + pow(x - theta, 2)));
#        k++;
#        if (k > 1000) {
#          Rcpp::Rcout << "error" << std::endl;
#          break;
#        }
#      }
#      return t;
#    }else{
#      Rcpp::Rcout << "please input the correct method" << std::endl;
#      return 0;
#    }
#  }

## -----------------------------------------------------------------------------
library(SA24204136)
set.seed(1234)
n<-1000
theta<-1
n1<-n2<-f1<-f2<-numeric(n)
for(i in 1:n){
  x<-rcauchy(n,location=1)
  n1[i]<-NR(x,theta,1)
  n2[i]<-NR(x,theta,2)
  f1[i]<-FSR(x,theta,1)
  f2[i]<-FSR(x,theta,2)
}

sdn1<-sqrt(sum((n1-mean(n1))^2)/(1000-1))
sdn2<-sqrt(sum((n2-mean(n2))^2)/(1000-1))
sdf1<-sqrt(sum((f1-mean(f1))^2)/(1000-1))
sdf2<-sqrt(sum((f2-mean(f2))^2)/(1000-1))

round(c(OnestepN_theta=mean(n1),OnestepN_sd=sdn1,MultistepN_theta=mean(n2),MultistepN_sd=sdn2,OnestepF_theta=mean(f1),OnestepF_sd=sdf1,MultistepF_theta=mean(f2),MultistepF_sd=sdf2),5)

## -----------------------------------------------------------------------------
n1<-n2<-f1<-f2<-numeric(n)
for(i in 1:n){
  x<-rcauchy(n,location=1)
  n1[i]<-NC(x,theta,1)
  n2[i]<-NC(x,theta,2)
  f1[i]<-FSC(x,theta,1)
  f2[i]<-FSC(x,theta,2)
}

sdn1<-sqrt(sum((n1-mean(n1))^2)/(1000-1))
sdn2<-sqrt(sum((n2-mean(n2))^2)/(1000-1))
sdf1<-sqrt(sum((f1-mean(f1))^2)/(1000-1))
sdf2<-sqrt(sum((f2-mean(f2))^2)/(1000-1))

round(c(OnestepN_theta=mean(n1),OnestepN_sd=sdn1,MultistepN_theta=mean(n2),MultistepN_sd=sdn2,OnestepF_theta=mean(f1),OnestepF_sd=sdf1,MultistepF_theta=mean(f2),MultistepF_sd=sdf2),5)

