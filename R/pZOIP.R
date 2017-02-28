#' Distribucion ZOIP
#'
#'La funcion pZOIP define la funcion de distribucion acumulada de la distribucion ZOIP.
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
#'. mu es el parametro de media y forma, ademas sigma es el parametro de dispersion de la distribucion.
#'family='F-C' distribucion Beta parametrizacion Ferrari y Cribari-Neto (2004), donde sigma=phi, phi es un parametro de precision.
#'family='Original' distribucion beta parametrizacion original donde mu=a, a parametro de forma 1; sigma=b, b parametro de forma 2.
#'family='simplex' distribucion simplex. propuesta por Barndorff-Nielsen and Jørgensen (1991)'
#'
#' @usage dZOIP(x, mu = 0.5, sigma = 0.1, p0 = 0.08333333, p1 = 0.08333333,family='R-S', log = FALSE)
#' pZOIP(q, mu = 0.5, sigma = 0.1, p0 = 0.08333333, p1 = 0.08333333,family='R-S',lower.tail = TRUE, log.p = FALSE)
#' qZOIP(p, mu = 0.5, sigma = 0.1, p0 = 0.08333333, p1 = 0.08333333,family='R-S',lower.tail = TRUE, log.p = FALSE)
#' rZOIP(n, mu = 0.5, sigma = 0.1,p0 = 0.08333333, p1 = 0.08333333,family='R-S')
#' @param q vector de cuantiles.
#' @param mu vector de parametros de localizacion.
#' @param sigma vector de parametros de escala.
#' @param p0 parametro de proporcion de ceros.
#' @param p1 parametro de proporcion de unos.
#' @param family eleccion de la parametrizacion o distribucion deseada, family='R-S' parametrizacion distribucion beta Rigby y Stasinopoulos, 'F-C' distribucion Beta parametrizacion Ferrari y Cribari-Neto, Original distribucion beta parametrizacion original, 'Simplex' distribucion simplex.
#' @param log.p logico; si TRUE, las probabilidades de p estaran dadas como log(p).
#' @param lower.tail logico; si TRUE (default), probabilidades estaran P[X <= x], en otro caso, P[X > x].
#' @export
#' @examples
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
  if(family != 'Simplex'){cdf <- ifelse((q > 0 & q < 1), nu + pbeta(q, shape1 = a,
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
