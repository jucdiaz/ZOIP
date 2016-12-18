#' Distribucion ZOIP
#'
#' La funcion ZOIP() define las distribuciones para datos proporcionales inflados en ceros
#' y unos mas conocidas, como la distribucion beta y simplex, ademas permite que se tengan
#' distintas parametrizaciones.La funcion ZOIP() posee 4 parametros dos de ellos es la
#' porporcion de ceros y unos que quiere que posea la distribucion para datos proporcionales
#' de su eleccion.
#'
#'Las funciones dZOIP, pZOIP, qZOIP y rZOIP define la densidad, la funcion de distribucion,
#'la funcion cuantil y la funcion generadora de numeros aleatorios para la distribucion ZOIP.
#'
#' @param  x,q vector de cuantiles
#' @param mu vector de parametros de localizacion
#' @param sigma vector de parametros de escala
#' @param p0 parametro de proporcion de ceros
#' @param p1 parametro de proporcion de unos
#' @param family eleccion de la parametrizacion o distribucion deseada
#' @param log,log.p logico; si TRUE, las probabilidades de p estaran dadas como log(p).
#' @param lower.tail logico; si TRUE (default), probabilidades estaran P[X <= x], en otro caso, P[X > x]
#' @param p vector de probabilidades.
#' @param n numero de observaciones. si length(n) > 1, el tamañano se tomara como el numero requerido
#' @export
#' @examples
#' dZOIP(x=0.5, mu = 0.2, sigma = 0.5, p0 = 0.2, p1 = 0.2,family='R-S',log = FALSE)
#' dZOIP(x=0.5, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='F-C',log = FALSE)
#' dZOIP(x=0.5, mu = 0.6, sigma = 2.4, p0 = 0.2, p1 = 0.2,family='Original',log = FALSE)
#' dZOIP(x=0.5, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='Simplex',log = FALSE)
#' pZOIP(q=0.5, mu = 0.2, sigma = 0.5, p0 = 0.2, p1 = 0.2,family='R-S',log = FALSE)
#' pZOIP(q=0.5, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='F-C',log = FALSE)
#' pZOIP(q=0.5, mu = 0.6, sigma = 2.4, p0 = 0.2, p1 = 0.2,family='Original',log = FALSE)
#' pZOIP(q=0.5, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='Simplex',log = FALSE)
#' qZOIP(p=0.7, mu = 0.2, sigma = 0.5, p0 = 0.2, p1 = 0.2,family='R-S',log = FALSE)
#' qZOIP(p=0.7, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='F-C',log = FALSE)
#' qZOIP(p=0.7, mu = 0.6, sigma = 2.4, p0 = 0.2, p1 = 0.2,family='Original',log = FALSE)
#' qZOIP(p=0.7, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='Simplex',log = FALSE)
#' a1<-rZOIP(n=1000, mu = 0.2, sigma = 0.5, p0 = 0.2, p1 = 0.2,family='R-S')
#' a2<-rZOIP(n=1000, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='F-C')
#' a3<-rZOIP(n=1000, mu = 0.6, sigma = 2.4, p0 = 0.2, p1 = 0.2,family='Original')
#' system.time(a4<-rZOIP(n=10, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='Simplex'))
#'
#' plot(density(a1))
#' plot(density(a2))
#' plot(density(a3))
#' plot(density(a4))

dZOIP<-function (x, mu = 0.5, sigma = 0.1, p0 = 0.08333333, p1 = 0.08333333,family='R-S', log = FALSE) {
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
  if (any(p0 <= 0) | any(p0 >= 1))
    stop(paste("p0 must be between 0 and 1", "\n", ""))
  if (any(p1 <= 0) | any(p1 >= 1))
    stop(paste("p1 must be between 0 and 1", "\n", ""))
  if (any(x < 0) | any(x > 1))
    stop(paste("x must be 0<=x<=1, i.e. 0 to 1 inclusively",
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

  nu <- p0/(1-p0-p1)
  tau <- p1/(1-p0-p1)
  logfy <- rep(0, length(x))
  if(family != 'Simplex'){logfy <- ifelse((x > 0 & x < 1), dbeta(x, shape1 = a, shape2 = b
                                                                 , ncp = 0, log = TRUE), 0)}
  if(family == 'Simplex'){logfy <- ifelse((x > 0 & x < 1), log(dsimplex(x, mu=mu, sig=sigma)), 0)}
  logfy <- ifelse((x == 0), log(nu), logfy)
  logfy <- ifelse((x == 1), log(tau), logfy)
  logfy <- logfy - log(1 + nu + tau)
  if (log == FALSE)
    fy <- exp(logfy)
  else fy <- logfy
  fy
}

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
  if (any(p0 <= 0) | any(p0 >= 1))
    stop(paste("p0 must be between 0 and 1", "\n", ""))
  if (any(p1 <= 0) | any(p1 >= 1))
    stop(paste("p1 must be between 0 and 1", "\n", ""))
  if (any(q < 0) | any(q > 1))
    stop(paste("y must be 0<=y<=1, i.e. 0 to 1 inclusively",
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

  nu <- p0/(1-p0-p1)
  tau <- p1/(1-p0-p1)

  if(family != 'Simplex'){cdf <- ifelse((q > 0 & q < 1), nu + pbeta(q, shape1 = a,
                                                                    shape2 = b, ncp = 0, lower.tail = TRUE, log.p = FALSE),
                                        0)}
  if(family == 'Simplex'){cdf <- ifelse((q > 0 & q < 1), nu + psimplex(q, mu=mu, sig=sigma),0)}
  cdf <- ifelse((q == 0), nu, cdf)
  cdf <- ifelse((q == 1), 1 + nu + tau, cdf)
  cdf <- cdf/(1 + nu + tau)
  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf = 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

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
  if (any(p0 <= 0) | any(p0 >= 1))
    stop(paste("p0 must be between 0 and 1", "\n", ""))
  if (any(p1 <= 0) | any(p1 >= 1))
    stop(paste("p1 must be between 0 and 1", "\n", ""))
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

  nu <- p0/(1-p0-p1)
  tau <- p1/(1-p0-p1)
  if(family != 'Simplex'){
    suppressWarnings(q <- ifelse((p <= (nu/(1 + nu + tau))),
                                 0, qbeta((p - (nu/(1 + nu + tau)))/(1/(1 + nu + tau)),
                                          shape1 = a, shape2 = b, lower.tail = TRUE, log.p = FALSE)))
  }
  if(family == 'Simplex'){
    suppressWarnings(q <- ifelse((p <= (nu/(1 + nu + tau))),
                                 0, qsimplex((p - (nu/(1 + nu + tau)))/(1/(1 + nu + tau)),
                                             mu=mu, sig=sigma)))
  }
  q <- ifelse((p >= ((1 + nu)/(1 + nu + tau))), 1, q)
  q
}

rZOIP<-function (n, mu = 0.5, sigma = 0.1,p0 = 0.08333333, p1 = 0.08333333,family='R-S')
{
  require('simplexreg')
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
  if (any(p0 <= 0) | any(p0 >= 1))
    stop(paste("p0 must be between 0 and 1", "\n", ""))
  if (any(p1 <= 0) | any(p1 >= 1))
    stop(paste("p1 must be between 0 and 1", "\n", ""))
  if (any(n <= 0))
    stop(paste("n must be a positive integer", "\n", ""))
  n <- ceiling(n)
  p <- runif(n)
  r <- qZOIP(p, mu = mu, sigma = sigma, p0 = p0, p1 = p1, family=family)
  r
}
