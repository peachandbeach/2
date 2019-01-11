#' @title InteofCauchy
#' @description a function to compute the cdf of the Cauchy distribution.
#' @param eta the location parameter
#' @param theta the scale parameter
#' @param q for distribution function F(q)
#' @return a value \code{n}
#' @examples
#' \dontrun{
#' InteofCauchy(22,0,1)
#' }
#' @export
InteofCauchy <- function(q,eta,theta){
f <- function(x, eta, theta) {
  1/(1+((x-eta)/theta)^2)/theta/pi
}
CDF <- integrate(f,lower=-Inf,upper=q,rel.tol=.Machine$double.eps^0.25,eta=eta,theta=theta)$value
CDF
}
NULL

