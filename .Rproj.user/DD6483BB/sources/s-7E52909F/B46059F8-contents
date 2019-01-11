#' @title GiniEstimation
#' @description Estimate by simulation the mean, median and deciles of hat{G} if X is continuous random variable
#' @param cdf for the supposed distribution(continuous)
#' @param mu for theoritical expectation
#' @param n for size of random numbers
#' @param m for repeating times
#' @param a for the better xlim for histogram(u can choose a better one when u observed its distribution in histogram with random a b )
#' @param b for the better xlim for histogram(u can choose a better one when u observed its distribution in histogram with random a b )
#' @return mean median and deciles ,a density histograms of the replicates\code{n}
#' @examples
#' \dontrun{
#' GiniEstimation(rlnorm,exp(0.5),1e3,1e2,0.45,0.65)#The standard lognormal
#' }
#' @export
GiniEstimation <-  function(cdf,mu,n,m,a,b)
{
  Gcdf<-numeric(m)
  set.seed(1234)
  v <- numeric(n)

  for(i in 1:m)
  {
    x<-sort(cdf(n))     #sort the x

    for(j in 1:n)
    {
      v[j] <- (2*j-n-1)*x[j]
    }
    s <- sum(v)
    Gcdf[i]<-(1/(n^2*mu))*s   #caculate the G for every i
  }
  mean<-mean(Gcdf)
  median<-quantile(Gcdf,0.5)
  deciles<-quantile(Gcdf,probs=seq(0,1,0.1))
  print(mean)
  print(median)
  print(deciles)
  hist(Gcdf,freq=F,xlim=c(a,b),main="Gini ratio of supposed distribution")  #
}
NULL
