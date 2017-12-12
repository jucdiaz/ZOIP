#' summary.ZOIP
#'
#' Summarize a ZOIP model.
#'
#' @param object An object of class \code{ZOIP}.
#' @param ... other arguments.
#'
#' @examples
#'
#' #Test 1--------------------------------------------------
#' library(ZOIP)
#' library(boot)
#' n<-1000
#' x1<-stats::runif(n)
#' x2<-stats::runif(n)
#'
#' b1<-0.3
#' b2<--0.5
#' b3<-0.9
#' sigma_i<-boot::inv.logit(b1+b2*x1+b3*x2)
#'
#' c1<-0.2
#' c2<--1
#' c3<-0.1
#' mu_i<-boot::inv.logit(c1+c2*x1)
#'
#' d1<-0.07
#' p0_i<-rep(d1,length(n))
#'
#' e1<-0.02
#' e2<--4
#' p1_i<-boot::inv.logit(e1+e2*x2)
#'
#' param<-cbind(mu_i,sigma_i,p0_i,p1_i)
#'
#' system.time(y_i<-apply(param,1,function(x){rZOIP(1,mu=x[1],sigma=x[2],
#'                                                  p0=x[3],p1=x[4],family='R-S')}))
#' data<-as.data.frame(cbind(y_i,x1,x2))
#'
#' formula.mu=y_i~x1
#' formula.sigma=~x1+x2
#' formula.p0=~1
#' formula.p1=~x1+x2
#' link=c('logit','logit','identity','logit')
#' family='R-S'
#' mod<-RM.ZOIP(formula.mu=formula.mu,formula.sigma=formula.sigma,formula.p0=formula.p0,
#'              formula.p1=formula.p1,data=data,link=link,family=family)
#' summary(mod)
#'
#' @export

summary.ZOIP<-function(object, ...){

  estimate <- object$par
  se       <- sqrt(diag(solve(object$HM)))
  zvalue   <- estimate / se
  pvalue   <- 2 * stats::pnorm(abs(zvalue), lower.tail=F)
  res      <- cbind(estimate=estimate, se=se, zvalue=zvalue, pvalue=pvalue)
  colnames(res) <- c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')
  res      <- as.data.frame(res)


  object$nparm[object$Vec_Bool==T]<--1

  a<-cumsum(object$nparm)

  if(object$Vec_Bool[1]==F)rownames(res[seq(1,a[1]),])<-names(object$par)[seq(1,a[1])]
  if(object$Vec_Bool[2]==F)rownames(res[seq(a[1]+1,a[2]),])<-names(object$par)[seq(a[1]+1,a[2])]
  if(object$Vec_Bool[3]==F)rownames(res[seq(a[2]+1,a[3]),])<-names(object$par)[seq(a[2]+1,a[3])]
  if(object$Vec_Bool[4]==F)rownames(res[seq(a[3]+1,a[4]),])<-names(object$par)[seq(a[3]+1,a[4])]

  Aux<-data.frame(t(c(0,NA,NA,NA,NA)))
  colnames(Aux) <- c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')
  rownames(Aux)<-c('(intercept)')

  if(object$Vec_Bool[1]==FALSE){elem.mu<-res[seq(1,a[1]),]}else elem.mu<-Aux
  if(object$Vec_Bool[2]==FALSE){elem.sigma<-res[seq(a[1]+1,a[2]),]}else elem.sigma<-Aux
  if(object$Vec_Bool[3]==FALSE){elem.p0<-res[seq(a[2]+1,a[3]),]}else elem.p0<-Aux
  if(object$Vec_Bool[4]==FALSE){elem.p1<-res[seq(a[3]+1,a[4]),]}else elem.p1<-Aux

  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for ",
            object$link[1], "(mu) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  stats::printCoefmat(elem.mu, P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for ",
            object$link[2], "(sigma) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  stats::printCoefmat(elem.sigma, P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for ",
            object$link[3], "(p0) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  stats::printCoefmat(elem.p0, P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for ",
            object$link[4], "(p1) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  stats::printCoefmat(elem.p1, P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat("---------------------------------------------------------------\n")

}

