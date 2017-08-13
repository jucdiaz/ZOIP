#' Distribucion ZOIP
#'
#'La funcion dZOIP define la funcion de densidad de probabilidad de la distribucion ZOIP.
#'
#'x tiene distribucion ZOIP con parametros de forma "mu", de escala "sigma", de proporcion de ceros "p0" y de proporcion de unos "p1",
#'si tiene densidad: p0 si x=0, p1 si x=1, (1-p0-p1)f(x;mu,sigma)) si 0<x<1.
#'
#'donde p0 ≥ 0 representa la probabilidad que x = 0, p1 ≥ 0 representa la probabilidad
#'de que x = 1, 0 ≤ p0 + p1 ≤ 1 y f(x; μ, sigma) representa algunas de las funciones de
#'densidad de probabilidad para datos proporcionales, como la distribucion beta con sus diferentes parametrizaciones
#'y la distribucion simplex.
#'
#'Cuando family='R-S' utiliza la distribucion beta con parametrizacion beta Rigby y Stasinopoulos (2008) el cual tiene una funcion de distribucion beta.
#'mu es el parametro de media y forma, ademas sigma es el parametro de dispersion de la distribucion.
#'family='F-C' distribucion Beta parametrizacion Ferrari y Cribari-Neto (2004), donde sigma=phi, phi es un parametro de precision.
#'family='Original' distribucion beta parametrizacion original donde mu=a, a parametro de forma 1; sigma=b, b parametro de forma 2.
#'family='simplex' distribucion simplex. propuesta por Barndorff-Nielsen and Jørgensen (1991)'
#'
#' @usage dZOIP(x, mu = 0.5, sigma = 0.1, p0 = 0.08333333, p1 = 0.08333333,family='R-S', log = FALSE)
#' pZOIP(q, mu = 0.5, sigma = 0.1, p0 = 0.08333333, p1 = 0.08333333,family='R-S',lower.tail = TRUE, log.p = FALSE)
#' qZOIP(p, mu = 0.5, sigma = 0.1, p0 = 0.08333333, p1 = 0.08333333,family='R-S',lower.tail = TRUE, log.p = FALSE)
#' rZOIP(n, mu = 0.5, sigma = 0.1,p0 = 0.08333333, p1 = 0.08333333,family='R-S')
#' @param x vector de cuantiles.
#' @param mu vector de parametros de localizacion.
#' @param sigma vector de parametros de escala.
#' @param p0 parametro de proporcion de ceros.
#' @param p1 parametro de proporcion de unos.
#' @param family eleccion de la parametrizacion o distribucion deseada, family='R-S' parametrizacion distribucion beta Rigby y Stasinopoulos, 'F-C' distribucion Beta parametrizacion Ferrari y Cribari-Neto, Original distribucion beta parametrizacion original, 'Simplex' distribucion simplex.
#' @param log logico; si TRUE, las probabilidades de p estaran dadas como log(p).
#' @export
#' @examples
#'dZOIP(x=0.5, mu = 0.2, sigma = 0.5, p0 = 0.2, p1 =0.2,family='R-S',log = FALSE)
#'dZOIP(x=0.5, mu = 0.2, sigma = 3, p0 = 0.2, p1 = 0.2,family='F-C',log = FALSE)
#'dZOIP(x=0.5, mu = 0.6, sigma = 2.4, p0 = 0.2, p1 = 0.2,family='Original',log = FALSE)
#'dZOIP(x=0.5, mu = 0.2, sigma = 3, p0 = 0, p1 = 0,family='Simplex',log = FALSE)
#'dZOIP(x=0.5, mu = 0.2, sigma = 0.5, p0 = 0.2, p1 =0,family='R-S',log = FALSE)
#'dZOIP(x=0.5, mu = 0.2, sigma = 0.5, p0 = 0, p1 =0.2,family='R-S',log = FALSE)
#'dZOIP(x=0.5, mu = 0.2, sigma = 0.5, p0 = 0, p1 =0,family='R-S',log = FALSE)
#'
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
  if(family != 'Simplex'){logfy <- ifelse((x > 0 & x < 1), dbeta(x, shape1 = a, shape2 = b
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

