#' @title a Monte Carlo estimate of the Beta function
#' @description a Monte Carlo estimate of the Beta(a, b) cdf
#' @param n the number of sampling times
#' @param a for Beta function
#' @param b for Beta function
#' @return a value \code{n}
#' @examples
#' \dontrun{
#' MCofBeta(1e4,3,3)
#' c(MCofBeta(1e4,3,3),beta(3,3))
#' #compare the MC estimation and the theoretical value
#' }
#' @export
MCofBeta <- function(n,a,b){
  x <- runif(n, 0, 1)    #generate random x
  theta.hat <- mean(x^(a-1)*(1-x)^(b-1))*1
  theta.hat
}
NULL

