#' Distribucion ZOIP
#'
#' La funcion qZOIP define la funcion cuantil de la distribucion ZOIP.
#'
#'x tiene distribucion ZOIP con parametros de forma "mu", de escala "sigma", de proporcion de ceros "p0" y de proporcion de unos "p1",
#'si tiene densidad: p0 si x=0, p1 si x=1, (1-p0-p1)f(x;mu,sigma)) si 0<x<1.
#'
#'donde p0 ≥ 0 representa la probabilidad que x = 0, p1 ≥ 0 representa la probabilidad
#'de que x = 1, 0 ≤ p0 + p1 ≤ 1 y f(x; μ, sigma) representa algunas de las funciones de
#'densidad de probabilidad para datos proporcionales, como la distribucion beta con sus diferentes parametrizaciones
#'y la distribucion simplex.
#'
#'Cuando family='R-S' se utiliza la distribucion beta con parametrizacion beta Rigby y Stasinopoulos (2008) el cual tiene una funcion de distribucion beta
#'f(x; μ, sigma) = B(μ, sigma)yμ((1−sigma^2)/sigma^2)−1(1 − y)(1−μ)((1−sigma^2)/sigma^2)−1. mu es el parametro de media y forma, ademas sigma es el parametro de dispersion de la distribucion.
#'family='F-C' distribucion Beta parametrizacion Ferrari y Cribari-Neto (2004), donde sigma=phi, phi es un parametro de precision.
#'family='Original' distribucion beta parametrizacion original donde mu=a, a parametro de forma 1; sigma=b, b parametro de forma 2.
#'family='simplex' distribucion simplex. propuesta por Barndorff-Nielsen and Jørgensen (1991)'
#'
#' @usage dZOIP(x, mu = 0.5, sigma = 0.1, p0 = 0.08333333, p1 = 0.08333333,family='R-S', log = FALSE)
#' @usage pZOIP(q, mu = 0.5, sigma = 0.1, p0 = 0.08333333, p1 = 0.08333333,family='R-S',lower.tail = TRUE, log.p = FALSE)
#' @usage qZOIP(p, mu = 0.5, sigma = 0.1, p0 = 0.08333333, p1 = 0.08333333,family='R-S',lower.tail = TRUE, log.p = FALSE)
#' @usage rZOIP(n, mu = 0.5, sigma = 0.1,p0 = 0.08333333, p1 = 0.08333333,family='R-S')
#' @param p vector de probabilidades.
#' @param mu vector de parametros de localizacion.
#' @param sigma vector de parametros de escala.
#' @param p0 parametro de proporcion de ceros.
#' @param p1 parametro de proporcion de unos.
#' @param family eleccion de la parametrizacion o distribucion deseada, family='R-S' parametrizacion distribucion beta Rigby y Stasinopoulos, 'F-C' distribucion Beta parametrizacion Ferrari y Cribari-Neto, Original distribucion beta parametrizacion original, 'Simplex' distribucion simplex.
#' @param log logico; si TRUE, las probabilidades de p estaran dadas como log(p).
#' @export
#' @examples
#' qZOIP(p=0.7, mu = 0.2, sigma = 0.5, p0 = 0.2, p1 = 0.2,family='R-S',log = FALSE)
#' qZOIP(p=0.7, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='F-C',log = FALSE)
#' qZOIP(p=0.7, mu = 0.6, sigma = 2.4, p0 = 0.2, p1 = 0.2,family='Original',log = FALSE)
#' qZOIP(p=0.7, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='Simplex',log = FALSE)


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
