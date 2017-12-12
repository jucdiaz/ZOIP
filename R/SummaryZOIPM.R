#' summary.ZOIPM
#'
#' Summarize a ZOIP model mixed.
#'
#' @param object An object of class \code{ZOIPM}.
#' @param ... other arguments.
#'
#' @examples
#'
#' library(ZOIP)
#' N<-15
#'
#' Times <- c(2, 10, 20)
#'
#' subject <- rep(1:N, each=length(Times))
#'
#' Days <- rep(Times, times=N)
#' b0i <- rep(rnorm(n=N,sd=1), each=length(Times))
#' b1i <- rep(rnorm(n=N,sd=0.5), each=length(Times))
#'
#' neta <- (1.6+b0i)-1.3*log(Days)
#' neta2<-(0.1+b1i)-0.8*log(Days)
#'
#'
#' mu <- 1 / (1 + exp(-neta))
#' sigma <- 1 / (1 + exp(-neta2))
#'
#' p0 <- 0.1
#' p1 <- 0.1
#'
#' mu[mu==1] <- 0.999
#' mu[mu==0] <- 0.001
#'
#' sigma[sigma==1] <- 0.999
#' sigma[sigma==0] <- 0.001
#'
#'
#' family<-'R-S'
#'
#' Y <- rZOIP(n=length(mu), mu = mu, sigma = sigma ,p0=p0,p1=p1,family=family)
#' base<-data.frame(Y,Days,subject)
#'
#'
#' formula.mu <- Y ~ log(Days)
#' formula.sigma<-~log(Days)
#' formula.p0<-~1
#' formula.p1<-~1
#'
#' formula.random=~1 | subject
#'
#' link<-c('logit','logit','identity','identity')
#'
#'
#' optimizer<-'nlminb'
#' n.points <-11
#' pruning <- TRUE
#'
#' mod<-RMM.ZOIP(formula.mu=formula.mu,formula.sigma=formula.sigma,formula.p0=formula.p0,
#'               formula.p1=formula.p1,data=base,formula.random=formula.random,link=link,
#'               family=family,optimizer=optimizer,n.points=n.points,pruning=pruning)
#' summary(mod)
#'
#' @export


summary.ZOIPM<-function(object, ...){

  estimate <- c(object$Fixed_Parameters.mu,object$Fixed_Parameters.sigma
                ,object$Fixed_Parameters.p0,object$Fixed_Parameters.p1,object$Parameters.randoms[,1])
  se       <- sqrt(diag(solve(object$HM)))
  zvalue   <- estimate / se
  pvalue   <- 2 * stats::pnorm(abs(zvalue), lower.tail=F)
  res      <- cbind(estimate=estimate, se=se, zvalue=zvalue, pvalue=pvalue)
  colnames(res) <- c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')
  res      <- as.data.frame(res)

  a <- 1:length(object$Fixed_Parameters.mu)
  b <-
    (length(object$Fixed_Parameters.mu) + 1):(length(object$Fixed_Parameters.mu) +
                                             length(object$Fixed_Parameters.sigma))
  c <-
    (length(object$Fixed_Parameters.mu) + length(object$Fixed_Parameters.sigma) +
       1):(
         length(object$Fixed_Parameters.mu) + length(object$Fixed_Parameters.sigma) + length(object$Fixed_Parameters.p0)
       )
  d <-
    (
      length(object$Fixed_Parameters.mu) + length(object$Fixed_Parameters.sigma) + length(object$Fixed_Parameters.p0) +
        1
    ):(
      length(object$Fixed_Parameters.mu) + length(object$Fixed_Parameters.sigma) + length(object$Fixed_Parameters.p0) +
        length(object$Fixed_Parameters.p1)
    )
  e <-
    (
      length(object$Fixed_Parameters.mu) + length(object$Fixed_Parameters.sigma) + length(object$Fixed_Parameters.p0) +
        length(object$Fixed_Parameters.p1) + 1
    ):(
      length(object$Fixed_Parameters.mu) + length(object$Fixed_Parameters.sigma) + length(object$Fixed_Parameters.p0) +
        length(object$Fixed_Parameters.p1) + length(object$Parameters.randoms[, 1])
    )
  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for ",
            object$link[1], "(mu) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  stats::printCoefmat(res[a,], P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for ",
            object$link[2], "(sigma) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  stats::printCoefmat(res[b,], P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for ",
            object$link[3], "(p0) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  stats::printCoefmat(res[c,], P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for ",
            object$link[4], "(p1) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  stats::printCoefmat(res[d,], P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat("---------------------------------------------------------------\n")
  cat(paste("Random effects for mu and sigma \n",sep=''))
  cat("---------------------------------------------------------------\n")
  stats::printCoefmat(res[e,], P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat("---------------------------------------------------------------\n")

}
