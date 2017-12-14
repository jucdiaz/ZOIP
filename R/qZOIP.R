#' ZOIP Distribution
#'
#' The qZOIP function defines the quantile function of the ZOIP distribution.
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
#' @param p vector of probabilities.
#' @param mu vector of location parameters.
#' @param sigma vector of scale parameters.
#' @param p0 parameter of proportion of zeros.
#' @param p1 Parameter of proportion of ones.
#' @param family choice of the parameterization or distribution, family = 'R-S' parameterization beta distribution Rigby and Stasinopoulos, 'F-C' distribution Beta parametrization Ferrari and Cribari-Neto, 'Original' Beta distribution classic parameterization, 'Simplex' simplex distribution.
#' @param lower.tail logical; if TRUE (default), probabilities will be P [X <= x], otherwise, P [X> x].
#' @param log.p logical; if TRUE, the probabilities of p will be given as log (p).
#' @examples
#' library(ZOIP)
#' qZOIP(p=0.7, mu = 0.2, sigma = 0.5, p0 = 0.2, p1 = 0.2,family='R-S',log = FALSE)
#' qZOIP(p=0.7, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='F-C',log = FALSE)
#' qZOIP(p=0.7, mu = 0.6, sigma = 2.4, p0 = 0.2, p1 = 0.2,family='Original',log = FALSE)
#' qZOIP(p=0.7, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='Simplex',log = FALSE)
#'
#' qZOIP(p=0.7, mu = 0.2, sigma = 0.5, p0 = 0.2, p1 = 0,family='R-S',log = FALSE)
#' qZOIP(p=0.7, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0,family='F-C',log = FALSE)
#' qZOIP(p=0.7, mu = 0.6, sigma = 2.4, p0 = 0.2, p1 = 0,family='Original',log = FALSE)
#' qZOIP(p=0.7, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0,family='Simplex',log = FALSE)
#'
#' qZOIP(p=0.7, mu = 0.2, sigma = 0.5, p0 = 0, p1 = 0.2,family='R-S',log = FALSE)
#' qZOIP(p=0.7, mu = 0.2, sigma = 3, p0 = 0, p1 = 0.2,family='F-C',log = FALSE)
#' qZOIP(p=0.7, mu = 0.6, sigma = 2.4, p0 = 0, p1 = 0.2,family='Original',log = FALSE)
#' qZOIP(p=0.7, mu = 0.2, sigma = 3, p0 = 0, p1 = 0.2,family='Simplex',log = FALSE)
#'
#' qZOIP(p=0.7, mu = 0.2, sigma = 0.5, p0 = 0, p1 = 0,family='R-S',log = FALSE)
#' qZOIP(p=0.7, mu = 0.2, sigma = 3, p0 = 0, p1 = 0,family='F-C',log = FALSE)
#' qZOIP(p=0.7, mu = 0.6, sigma = 2.4, p0 = 0, p1 = 0,family='Original',log = FALSE)
#' qZOIP(p=0.7, mu = 0.2, sigma = 3, p0 = 0, p1 = 0,family='Simplex',log = FALSE)
#' @export

qZOIP<-function (p, mu = 0.5, sigma = 0.1, p0 = 0.08333333, p1 = 0.08333333,family='R-S',
                 lower.tail = TRUE, log.p = FALSE){
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
  if (any(p < 0) | any(p > 1))
    stop(paste("p must be between 0 and 1", "\n", ""))
  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p

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


  if(family != 'Simplex'){
    if(p0>0 && p1>0){
      suppressWarnings(q <- ifelse((p <= (nu/(1 + nu + tau))),
                                   0, stats::qbeta((p - (nu/(1 + nu + tau)))/(1/(1 + nu + tau)),
                                            shape1 = a, shape2 = b, lower.tail = TRUE, log.p = FALSE)))
    }else if(p0>0 && p1==0){
      suppressWarnings(q <- ifelse((p <= (nu/(1 + nu))), 0, stats::qbeta((p -
                                                                     (nu/(1 + nu)))/(1/(1 + nu)), shape1 = a, shape2 = b,
                                                                  lower.tail = TRUE, log.p = FALSE)))
    }else if(p0==0 && p1>0){
      suppressWarnings(q <- ifelse((p >= (1/(1 + tau))), 1, stats::qbeta((p *
                                                                     (1 + tau)), shape1 = a, shape2 = b, lower.tail = TRUE,
                                                                  log.p = FALSE)))
    }else if(p0==0 && p1==0){
      suppressWarnings(q <-stats::qbeta(p,shape1 = a, shape2 = b, lower.tail = TRUE, log.p = FALSE))
    }
  }
  fun_simp_aux1<-function(p){
    ifelse((p <= (nu/(1 + nu + tau))),0,
           ifelse(((p - (nu/(1 + nu + tau)))/(1/(1 + nu + tau))<=0) || ((p - (nu/(1 + nu + tau)))/(1/(1 + nu + tau))>=1)
                  , NA,rmutil::qsimplex((p - (nu/(1 + nu + tau)))/(1/(1 + nu + tau)),m=mu, s=sigma)))
  }
  fun_simp_aux2<-function(p){
    ifelse((p <= (nu/(1 + nu))),0,
           ifelse(((p -(nu/(1 + nu)))/(1/(1 + nu))<=0) | ((p -(nu/(1 + nu)))/(1/(1 + nu))>=1)
                  , NA,rmutil::qsimplex((p -(nu/(1 + nu)))/(1/(1 + nu)),m=mu, s=sigma)))
  }

  fun_simp_aux3<-function(p){
    ifelse((p >= (1/(1 + tau))),1,
           ifelse((p *(1 + tau)<=0) | (p *(1 + tau)>=1)
                  , NA,rmutil::qsimplex(p *(1 + tau),m=mu, s=sigma)))
  }

  if(family == 'Simplex'){
    if(p0>0 && p1>0){ suppressWarnings(q <- apply(as.matrix(p),1,fun_simp_aux1))
    }else if(p0>0 && p1==0){
      suppressWarnings(q <- apply(as.matrix(p),1,fun_simp_aux2))
    }else if(p0==0 && p1>0){
      suppressWarnings(q <- apply(as.matrix(p),1,fun_simp_aux3))
    }else if(p0==0 && p1==0){
      suppressWarnings(q <-rmutil::qsimplex(p,m=mu, s=sigma))
    }
  }
  if(p0>0 && p1>0){
    q <- ifelse((p >= ((1 + nu)/(1 + nu + tau))), 1, q)
  }
  q
}

