#' Distribucion ZOIP
#'
#'La funcion dZOIP define la funcion de densidad de probabilidad de la distribucion ZOIP.
#'
#' @param x vector de cuantiles
#' @param mu vector de parametros de localizacion
#' @param sigma vector de parametros de escala
#' @param p0 parametro de proporcion de ceros
#' @param p1 parametro de proporcion de unos
#' @param family eleccion de la parametrizacion o distribucion deseada
#' @param log logico; si TRUE, las probabilidades de p estaran dadas como log(p).
#' @export
#' @examples
#' dZOIP(x=0.5, mu = 0.2, sigma = 0.5, p0 = 0.2, p1 = 0.2,family='R-S',log = FALSE)
#' dZOIP(x=0.5, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='F-C',log = FALSE)
#' dZOIP(x=0.5, mu = 0.6, sigma = 2.4, p0 = 0.2, p1 = 0.2,family='Original',log = FALSE)
#' dZOIP(x=0.5, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='Simplex',log = FALSE)
#'
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
