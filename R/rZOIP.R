#' Distribucion ZOIP
#'
#' la funcion rZOIP define la funcion generadora de numeros aleatorios para la distribucion ZOIP.
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
#' pZOIP(q, mu = 0.5, sigma = 0.1, p0 = 0.08333333, p1 = 0.08333333,family='R-S',lower.tail = TRUE, log.p = FALSE)
#' qZOIP(p, mu = 0.5, sigma = 0.1, p0 = 0.08333333, p1 = 0.08333333,family='R-S',lower.tail = TRUE, log.p = FALSE)
#' rZOIP(n, mu = 0.5, sigma = 0.1,p0 = 0.08333333, p1 = 0.08333333,family='R-S')
#' @param n numero de observaciones. si length(n) > 1, el tamañano se tomara como el numero requerido.
#' @param mu vector de parametros de localizacion.
#' @param sigma vector de parametros de escala.
#' @param p0 parametro de proporcion de ceros.
#' @param p1 parametro de proporcion de unos.
#' @param family eleccion de la parametrizacion o distribucion deseada, family='R-S' parametrizacion distribucion beta Rigby y Stasinopoulos, 'F-C' distribucion Beta parametrizacion Ferrari y Cribari-Neto, Original distribucion beta parametrizacion original, 'Simplex' distribucion simplex.
#' @export
#' @examples
#' a1<-rZOIP(n=1000, mu = 0.2, sigma = 0.5, p0 = 0.2, p1 = 0.2,family='R-S')
#' a2<-rZOIP(n=1000, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='F-C')
#' a3<-rZOIP(n=1000, mu = 0.6, sigma = 2.4, p0 = 0.2, p1 = 0.2,family='Original')
#' system.time(a4<-rZOIP(n=10, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='Simplex'))
#' plot(density(a1))
#' plot(density(a2))
#' plot(density(a3))
#' plot(density(a4))


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
