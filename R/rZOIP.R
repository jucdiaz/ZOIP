#' Distribucion ZOIP
#'
#' la funcion rZOIP define la funcion generadora de numeros aleatorios para la distribucion ZOIP.
#'
#' @param n numero de observaciones. si length(n) > 1, el tamañano se tomara como el numero requerido
#' @param mu vector de parametros de localizacion
#' @param sigma vector de parametros de escala
#' @param p0 parametro de proporcion de ceros
#' @param p1 parametro de proporcion de unos
#' @param family eleccion de la parametrizacion o distribucion deseada
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
