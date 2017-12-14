#' ZOIP Distribution
#'
#' The rZOIP function defines the random number generating function for the ZOIP distribution.
#'
#'x has ZOIP distribution with shape parameters "\eqn{\mu}", scale "\eqn{\sigma}", proportion of zeros "\eqn{p0}" and proportion of ones "\eqn{ p1} ",
#'has density: \eqn{p0} if \eqn{x = 0}, \eqn{p1} if \eqn{x = 1}, \eqn{(1-p0-p1) f (x; \mu, \sigma)} yes \eqn{0 <x <1}.
#'
#'where \eqn{p0 \ge 0} represents the probability that \eqn{x = 0, p1 \ge 0} represents the probability
#'that \eqn{x = 1, 0 \le p0 + p1 \le 1} and \eqn{f (x; \mu, \sigma)} represents some of the functions of
#'probability density for proportional data, such as the beta distribution with its different parameterizations
#'and the simplex distribution.
#'
#'When family =' R-S 'uses the beta distribution with beta parameterization Rigby and Stasinopoulos (2005) which has a beta distribution function.
#'\eqn{\mu} is the parameter of mean and shape, plus \eqn{\sigma} is the dispersion parameter of the distribution.
#'family =' F-C 'distribution Beta parametrization Ferrari and Cribari-Neto (2004), where \eqn{\sigma = \phi}, \eqn{\phi} is a precision parameter.
#'family =' Original 'beta distribution original parametrization where \eqn{\mu = a}, a parameter of form 1; \eqn{\sigma = b}, b parameter of form 2.
#'family =' Simplex 'simplex distribution. proposed by Barndorff-Nielsen and Jørgensen (1991)
#'
#' @param n number of observations. If length (n)> 1, the length is taken to be the number required.
#' @param mu vector of location parameters.
#' @param sigma vector of scale parameters.
#' @param p0 parameter of proportion of zeros.
#' @param p1 Parameter of proportion of ones.
#' @param family choice of the parameterization or distribution, family = 'R-S' parameterization beta distribution Rigby and Stasinopoulos, 'F-C' distribution Beta parametrization Ferrari and Cribari-Neto, 'Original' Beta distribution classic parameterization, 'Simplex' simplex distribution.
#' @examples
#' library(ZOIP)
#' a1<-rZOIP(n=1000, mu = 0.2, sigma = 0.5, p0 = 0.2, p1 = 0.2,family='R-S')
#' a2<-rZOIP(n=1000, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='F-C')
#' a3<-rZOIP(n=1000, mu = 0.6, sigma = 2.4, p0 = 0.2, p1 = 0.2,family='Original')
#' system.time(a4<-rZOIP(n=10, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='Simplex'))
#'
#' plot(density(a1))
#' plot(density(a2))
#' plot(density(a3))
#' plot(density(a4))
#'
#' a1<-rZOIP(n=1000, mu = 0.2, sigma = 0.5, p0 = 0.2, p1 = 0,family='R-S')
#' a2<-rZOIP(n=1000, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0,family='F-C')
#' a3<-rZOIP(n=1000, mu = 0.6, sigma = 2.4, p0 = 0.2, p1 = 0,family='Original')
#' system.time(a4<-rZOIP(n=10, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0,family='Simplex'))
#'
#' plot(density(a1))
#' plot(density(a2))
#' plot(density(a3))
#' plot(density(a4))
#'
#' a1<-rZOIP(n=1000, mu = 0.2, sigma = 0.5, p0 = 0, p1 = 0.2,family='R-S')
#' a2<-rZOIP(n=1000, mu = 0.2, sigma = 3, p0 = 0, p1 = 0.2,family='F-C')
#' a3<-rZOIP(n=1000, mu = 0.6, sigma = 2.4, p0 = 0, p1 = 0.2,family='Original')
#' system.time(a4<-rZOIP(n=10, mu = 0.2, sigma = 3, p0 = 0, p1 = 0.2,family='Simplex'))
#'
#' plot(density(a1))
#' plot(density(a2))
#' plot(density(a3))
#' plot(density(a4))
#' @export

rZOIP<-function (n, mu = 0.5, sigma = 0.1,p0 = 0.08333333, p1 = 0.08333333,family='R-S')
{
  if (any(family != 'R-S') && any(family != 'F-C') && any(family != 'Original') && any(family != 'Simplex'))
    stop(paste("family must be in R-S, F-C, Original, Simplex", "\n", ""))
  if (any(family != 'Original') && (any(mu <= 0) | any(mu >= 1)))
    stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(family == 'Original') && any(mu <= 0))
    stop(paste("mu is shape1 must greated than 0", "\n", ""))
  if (any(family == 'R-S') && (any(sigma <= 0) | any(sigma >= 1)))
    stop(paste("sigma must be between 0 and 1", "\n", ""))
  if (any(family != 'R-S') && any(sigma <= 0))
    stop(paste("sigma must greated than 0", "\n", ""))
  if (any(p0>=1))
    stop(paste("p0 must be lower than 1", "\n", ""))
  if (any(p1>=1))
    stop(paste("p1 must be lower than 1", "\n", ""))
  if (any(p0+p1 < 0) | any(p0+p1 >= 1))
    stop(paste("p0+p1 must be between 0 and 1", "\n", ""))
  if (any(p0 < 0) | any(p0 > 1))
    stop(paste("p0 must be between 0 and 1", "\n", ""))
  if (any(n <= 0))
    stop(paste("n must be a positive integer", "\n", ""))
  n <- ceiling(n)
  p <- stats::runif(n)
  r <- qZOIP(p, mu = mu, sigma = sigma, p0 = p0, p1 = p1, family=family)
  r
}
