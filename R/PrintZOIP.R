#' print.ZOIP
#'
#' print a ZOIP model.
#'
#' @param x An object of class \code{ZOIP}.
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
#' mod<-RM.ZOIP(formula.mu=formula.mu,formula.sigma=formula.sigma,
#'              formula.p0=formula.p0,formula.p1=formula.p1,data=data,link=link,family=family)
#' mod
#'
#'
#' @export

print.ZOIP<-function(x, ...){

  x$nparm[x$Vec_Bool==T]<--1

  a<-cumsum(x$nparm)

  Aux<-c(0)
  names(Aux)<-c('(intercept)')

  if(x$Vec_Bool[1]==FALSE){elem.mu<-x$par[seq(1,a[1])]}else elem.mu<-Aux
  if(x$Vec_Bool[2]==FALSE){elem.sigma<-x$par[seq(a[1]+1,a[2])]}else elem.sigma<-Aux
  if(x$Vec_Bool[3]==FALSE){elem.p0<-x$par[seq(a[2]+1,a[3])]}else elem.p0<-Aux
  if(x$Vec_Bool[4]==FALSE){elem.p1<-x$par[seq(a[3]+1,a[4])]}else elem.p1<-Aux

  cat("Call:\n")
  print(x$call)
  cat("\n Results: \n")
  cat("\n Estimated coefficients for h(mu): \n")
  print(elem.mu)
  cat("\n Estimated coefficients for h(sigma): \n")
  print(elem.sigma)
  cat("\n Estimated coefficients for h(p0): \n")
  print(elem.p0)
  cat("\n Estimated coefficients for h(p1): \n")
  print(elem.p1)
  cat("\n Convergence \n")
  print(x$Convergence)
  cat("\n message \n")
  print(x$message)
  cat("\n iterations \n")
  print(x$iterations)
  cat("\n Log-likelihood \n")
  print(x$objective)

}
