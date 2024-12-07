
#' @importFrom Rcpp evalCpp
#' @importFrom stats rcauchy
#' @useDynLib SA24204136
NULL

#' @title Newton-Raphson algorithm for cauchy(theta) using R.
#' @description Using three input to get the result of Newton-Raphson algorithm of cauchy distribution.
#' @param x the data of cauchy distribution
#' @param theta the initial value 
#' @param method method=1 for one step iteration and method=2 for multi-step iteration
#' @return the estimation of param theta
#' @examples
#' \dontrun{
#' x<-rcauchy(200,location=1)
#' theta<-median(x)
#' thetahat<-NR(x,theta,1)
#' print(thetahat)
#' }
#' @export
NR<-function(x,theta,method){
  if(method==1){
    return(theta+sum(2*(x-theta)/(1+(x-theta)^2))/sum(2*(1-(x-theta)^2)/(1+(x-theta)^2)^2))
  }
  if(method==2){
    t<-theta+sum(2*(x-theta)/(1+(x-theta)^2))/sum(2*(1-(x-theta)^2)/(1+(x-theta)^2)^2)
    k<-1
    theta<-theta
    while(abs(theta-t)>1e-3){
      theta<-t
      t<-theta+sum(2*(x-theta)/(1+(x-theta)^2))/sum(2*(1-(x-theta)^2)/(1+(x-theta)^2)^2)
      k<-k+1
      if(k>1000){
        print("error")
        break
      }
    }
    return(t)
  }
  if(method!=1 && method !=2){
    print("please input the correct method")
  }
}

#' @title Fisher scoring algorithm for cauchy(theta) using R.
#' @description Using three input to get the result of Fisher scoring algorithm of cauchy distribution.
#' @param x the data of cauchy distribution
#' @param theta the initial value 
#' @param method method=1 for one step iteration and method=2 for multi-step iteration
#' @return the estimation of param theta
#' @examples
#' \dontrun{
#' x<-rcauchy(200,location=1)
#' theta<-median(x)
#' thetahat<-FSR(x,theta,1)
#' print(thetahat)
#' }
#' @export
FSR<-function(x,theta,method){
  if(method==1){
    return(theta+2*mean(2*(x-theta)/(1+(x-theta)^2)))
  }
  if(method==2){
    t<-theta+2*mean(2*(x-theta)/(1+(x-theta)^2))
    k<-1
    theta<-theta
    while(abs(theta-t)>1e-3){
      theta<-t
      t<-theta+2*mean(2*(x-theta)/(1+(x-theta)^2))
      k<-k+1
      if(k>1000){
        print("error")
        break
      }
    }
    return(t)
  }
  if(method!=1 && method !=2){
    print("please input the correct method")
  }
}


