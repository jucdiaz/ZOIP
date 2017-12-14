#' ZOIP Distribution
#'
#'The pZOIP function defines the cumulative distribution function of the ZOIP distribution.
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
#' @param q quantiles vector.
#' @param mu vector of location parameters.
#' @param sigma vector of scale parameters.
#' @param p0 parameter of proportion of zeros.
#' @param p1 Parameter of proportion of ones.
#' @param family choice of the parameterization or distribution, family = 'R-S' parameterization beta distribution Rigby and Stasinopoulos, 'F-C' distribution Beta parametrization Ferrari and Cribari-Neto, 'Original' Beta distribution classic parameterization, 'Simplex' simplex distribution.
#' @param lower.tail logical; if TRUE (default), probabilities will be P [X <= x], otherwise, P [X> x].
#' @param log.p logical; if TRUE, the probabilities of p will be given as log (p).
#' @examples
#' library(ZOIP)
#' pZOIP(q=0.5, mu = 0.2, sigma = 0.5, p0 = 0.2, p1 = 0.2,family='R-S',log = FALSE)
#' pZOIP(q=0.5, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='F-C',log = FALSE)
#' pZOIP(q=0.5, mu = 0.6, sigma = 2.4, p0 = 0.2, p1 = 0.2,family='Original',log = FALSE)
#' pZOIP(q=0.5, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='Simplex',log = FALSE)
#'
#' pZOIP(q=0.5, mu = 0.2, sigma = 0.5, p0 = 0, p1 = 0.2,family='R-S',log = FALSE)
#' pZOIP(q=0.5, mu = 0.2, sigma = 3, p0 = 0, p1 = 0.2,family='F-C',log = FALSE)
#' pZOIP(q=0.5, mu = 0.6, sigma = 2.4, p0 = 0, p1 = 0.2,family='Original',log = FALSE)
#' pZOIP(q=0.5, mu = 0.2, sigma = 3, p0 = 0, p1 = 0.2,family='Simplex',log = FALSE)
#'
#' pZOIP(q=0.5, mu = 0.2, sigma = 0.5, p0 = 0.2, p1 = 0,family='R-S',log = FALSE)
#' pZOIP(q=0.5, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0,family='F-C',log = FALSE)
#' pZOIP(q=0.5, mu = 0.6, sigma = 2.4, p0 = 0.2, p1 = 0,family='Original',log = FALSE)
#' pZOIP(q=0.5, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0,family='Simplex',log = FALSE)
#'
#' pZOIP(q=0.5, mu = 0.2, sigma = 0.5, p0 = 0, p1 = 0,family='R-S',log = FALSE)
#' pZOIP(q=0.5, mu = 0.2, sigma = 3, p0 = 0, p1 = 0,family='F-C',log = FALSE)
#' pZOIP(q=0.5, mu = 0.6, sigma = 2.4, p0 = 0, p1 = 0,family='Original',log = FALSE)
#' pZOIP(q=0.5, mu = 0.2, sigma = 3, p0 = 0, p1 = 0,family='Simplex',log = FALSE)
#' @export

pZOIP<-function (q, mu = 0.5, sigma = 0.1, p0 = 0.08333333, p1 = 0.08333333,family='R-S',
                 lower.tail = TRUE, log.p = FALSE)
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
  if (any(p1 < 0) | any(p1 > 1))
    stop(paste("p1 must be between 0 and 1", "\n", ""))
  if (any(q < 0) | any(q > 1))
    stop(paste("y must be 0<=y<=1, i.e. 0 to 1 inclusively",
               "\n", ""))
  if ((any(q==0) | any(q==1)) && any(p0==0) && any(p1==0) )
    stop(paste("y must be 0<y<1, desity is not inflated",
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
  }else if(p0==0 && p1==0){
    nu<-0
    tau<-0
  }
  if(family != 'Simplex'){cdf <- ifelse((q > 0 & q < 1), nu + stats::pbeta(q, shape1 = a,
                                                                    shape2 = b, ncp = 0, lower.tail = TRUE, log.p = FALSE),
                                        0)}
  if(family == 'Simplex'){cdf <- ifelse((q > 0 & q < 1), nu + rmutil::psimplex(q, m=mu, s=sigma),0)}
  cdf <- ifelse((q == 0), nu, cdf)
  cdf <- ifelse((q == 1), 1 + nu + tau, cdf)
  if(p0>0 && p1>0){
    cdf <- cdf/(1 + nu + tau)
  }else if (p0>0 && p1==0){cdf <- cdf/(1 + nu)
  }else if (p0==0 && p1>0){cdf <- cdf/(1 + tau)
  }else if (p0==0 && p1==0){ cdf <- cdf}
  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf = 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
