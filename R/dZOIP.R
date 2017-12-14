#' ZOIP Distribution
#'
#'The dZOIP function defines the probability density function of the ZOIP distribution.
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
#' @param x quantiles vector.
#' @param mu vector of location parameters.
#' @param sigma vector of scale parameters.
#' @param p0 parameter of proportion of zeros.
#' @param p1 Parameter of proportion of ones.
#' @param family choice of the parameterization or distribution, family = 'R-S' parameterization beta distribution Rigby and Stasinopoulos, 'F-C' distribution Beta parametrization Ferrari and Cribari-Neto, 'Original' Beta distribution classic parameterization, 'Simplex' simplex distribution.
#' @param log logical; if TRUE, the probabilities of p will be given as log (p).
#' @examples
#' library(ZOIP)
#' dZOIP(x=0.5, mu = 0.2, sigma = 0.5, p0 = 0.2, p1 =0.2,family='R-S',log = FALSE)
#' dZOIP(x=0.5, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='F-C',log = FALSE)
#' dZOIP(x=0.5, mu = 0.6, sigma = 2.4, p0 = 0.2, p1 = 0.2,family='Original',log = FALSE)
#' dZOIP(x=0.5, mu = 0.2, sigma = 3, p0 = 0, p1 = 0,family='Simplex',log = FALSE)
#' dZOIP(x=0.5, mu = 0.2, sigma = 0.5, p0 = 0.2, p1 =0,family='R-S',log = FALSE)
#' dZOIP(x=0.5, mu = 0.2, sigma = 0.5, p0 = 0, p1 =0.2,family='R-S',log = FALSE)
#' dZOIP(x=0.5, mu = 0.2, sigma = 0.5, p0 = 0, p1 =0,family='R-S',log = FALSE)
#' @export


dZOIP<-function (x, mu = 0.5, sigma = 0.1, p0 = 0.08333333, p1 = 0.08333333,family='R-S', log = FALSE) {
  if (any(family != 'R-S') && any(family != 'F-C') && any(family != 'Original') && any(family != 'Simplex'))
    stop(paste("family must be in R-S, F-C, Original, Simplex", "\n", ""))
  if (any(family != 'Original' && family != 'Simplex') && (any(mu <= 0) | any(mu >= 1)))
    stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(family == 'Original') && any(mu <= 0))
    stop(paste("mu is shape1 must higher than 0", "\n", ""))
  if (any(family == 'R-S') && (any(sigma <= 0) | any(sigma >= 1)))
    stop(paste("sigma must be between 0 and 1", "\n", ""))
  if (any(family != 'R-S') && any(sigma <= 0))
    stop(paste("sigma must higher than 0", "\n", ""))
  if (any(p0>=1))
    stop(paste("p0 must be lower than 1", "\n", ""))
  if (any(p1>=1))
    stop(paste("p1 must be lower than 1", "\n", ""))
  if ((any(p0+p1 < 0) | any(p0+p1 >1)))
   stop(paste("p0+p1 must be between 0 and 1", "\n", ""))
  if (any(p0 < 0) | any(p0 > 1))
    stop(paste("p0 must be between 0 and 1", "\n", ""))
  if (any(p1 < 0) | any(p1 > 1))
    stop(paste("p1 must be between 0 and 1", "\n", ""))
  if (any(x < 0) | any(x > 1))
    stop(paste("x must be 0<=x<=1, i.e. 0 to 1 inclusively",
               "\n", ""))
  if ((any(x==0) | any(x==1)) && any(p0==0) && any(p1==0) )
    stop(paste("x must be 0<x<1, desity is not inflated",
               "\n", ""))

  if(family == 'R-S'){
    a <- mu * (1 - sigma^2)/(sigma^2)
    b <- a * (1 - mu)/mu
  }else if(family == 'F-C'){
    a <- mu * sigma
    b <- (1 - mu) * sigma
  }else if(family == 'Original'){
    a <- mu
    b <- sigma
  }
  if(p0>0 || p1>0){
    nu <- p0/(1-p0-p1)
    tau <- p1/(1-p0-p1)
  }
  logfy <- rep(0, length(x))
  if(family != 'Simplex'){logfy <- ifelse((x > 0 & x < 1), stats::dbeta(x, shape1 = a, shape2 = b
                                                                 , ncp = 0, log = TRUE), 0)}

  if(family == 'Simplex'){logfy <- ifelse((x > 0 & x < 1), dsim2(x, m = mu, s = sigma, log = TRUE), 0)}
  logfy <- ifelse((x == 0), log(nu), logfy)
  logfy <- ifelse((x == 1), log(tau), logfy)
  if(p0>0 && p1>0){
    logfy <- logfy - log(1 + nu + tau)
  }else if (p0>0 && p1==0){ logfy <- logfy - log(1 + nu)
  }else if (p0==0 && p1>0){ logfy <- logfy - log(1 + tau)
  }else if (p0==0 && p1==0){ logfy <- logfy}
  if (log == FALSE){
    fy <- exp(logfy)
  }else fy <- logfy
  fy
}

dsim2 <- function (y, m, s, log = FALSE)
{
  #if (any(m < 0) || any(m > 1))
  #  stop("m must contain values between 0 and 1")
  if (any(s <= 0))
    stop("s must be positive")
  tmp <- -((y - m) / (m * (1 - m))) ^ 2 / (2 * y * (1 - y) * s) -
    (log(2 * pi * s) + 3 * (log(y) + log(1 - y))) / 2
  if (!log)
    tmp <- exp(tmp)
  tmp
}

