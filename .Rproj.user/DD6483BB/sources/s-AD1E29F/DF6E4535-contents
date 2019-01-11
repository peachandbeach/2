#' @title  compute the cdf of the Cauchy distribution
#' @description Use a function to compute the cdf of the Cauchy distribution, theta and eta are scale and location parameter respectively
#' @param theta  scale parameter
#' @param eta  location parameter
#' @param y  vector of quantiles.
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' theta = 3
#' eta = 2
#' y <- seq(-10, 10, length(10))
#' Cau(theta, eta, y)
#' }
#' @export
Cau <- function(theta, eta, y)
{
  .dcauchy <- function(theta, eta, x)
  {
    #compute the density function of cauchy distribution
    if(theta <= 0) return(NA)

    r <- theta * pi * (1 + ((x-eta)/theta)^2)
    return(1 / r)
  }

  n <- length(y)
  res <- numeric(n)
  for(i in 1:n)
  {
    res[i] <- integrate(.dcauchy, lower = -Inf, upper = y[i],
                        rel.tol = .Machine$double.eps^0.25,
                        theta = theta, eta = eta)$value
  }
  return(res)
}

