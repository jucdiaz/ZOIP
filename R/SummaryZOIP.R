#' summary.ZOIP
#'
#' Summarize a ZOIP model.
#'
#' @usage summary(mod)
#' @param mod An object of class \code{ZOIP}.
#'
#' @examples
#'
#' #Test 1--------------------------------------------------
#'
#' n<-1000
#' x1<-runif(n)
#' x2<-runif(n)
#'
#' b1<-0.3
#' b2<--0.5
#' b3<-0.9
#' sigma_i<-exp(b1+b2*x1+b3*x2)
#'
#' c1<-0.2
#' c2<--1
#' c3<-0.1
#' mu_i<-exp(c1+c2*x1)
#'
#' d1<-0.07
#' p0_i<-rep(d1,length(n))
#'
#' e1<-0.02
#' e2<--4
#' p1_i<-inv.logit(e1+e2*x2)
#'
#' param<-cbind(mu_i,sigma_i,p0_i,p1_i)
#'
#' system.time(y_i<-apply(param,1,function(x){rZOIP(1,mu=x[1],sigma=x[2],p0=x[3],p1=x[4],family='R-S')}))
#' data<-as.data.frame(cbind(y_i,x1,x2))
#'
#' formula.mu=y_i~x1
#' formula.sigma=~x1+x2
#' formula.p0=~1
#' formula.p1=~x1+x2
#' link=c('logit','logit','identity','logit')
#' family='R-S'
#' mod<-RM.ZOIP(formula.mu=formula.mu,formula.sigma=formula.sigma,formula.p0=formula.p0,formula.p1=formula.p1,data=data,link=link,family=family)
#' mod
#' summary(mod)
#'
#' @export

summary.ZOIP<-function(mod){

  estimate <- mod$par
  se       <- sqrt(diag(solve(mod$HM)))
  zvalue   <- estimate / se
  pvalue   <- 2 * pnorm(abs(zvalue), lower.tail=F)
  res      <- cbind(estimate=estimate, se=se, zvalue=zvalue, pvalue=pvalue)
  colnames(res) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
  res      <- as.data.frame(res)

  mod$nparm[mod$Vec_Bool==T]<--1

  a1<-cumsum(mod$nparm)[1]
  a2<-cumsum(mod$nparm)[2]
  a3<-cumsum(mod$nparm)[3]
  a4<-cumsum(mod$nparm)[4]

  if(mod$Vec_Bool[1]==F)rownames(res[seq(1,a1+1),])<-names(mod$par)[seq(1,a1+1)]
  if(mod$Vec_Bool[2]==F)rownames(res[seq(a1+2,a2+2),])<-names(mod$par)[seq(a1+2,a2+2)]
  if(mod$Vec_Bool[3]==F)rownames(res[seq(a2+3,a3+3),])<-names(mod$par)[seq(a2+3,a3+3)]
  if(mod$Vec_Bool[4]==F)rownames(res[seq(a3+4,a4+4),])<-names(mod$par)[seq(a3+4,a4+4)]

  Aux<-data.frame(t(c(0,NA,NA,NA,NA)))
  colnames(Aux) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
  rownames(Aux)<-c('(intercept)')

  if(mod$Vec_Bool[1]==F){elem.mu<-res[seq(1,a1+1),]}else elem.mu<-Aux
  if(mod$Vec_Bool[2]==F){elem.sigma<-res[seq(a1+2,a2+2),]}else elem.sigma<-Aux
  if(mod$Vec_Bool[3]==F){elem.p0<-res[seq(a2+3,a3+3),]}else elem.p0<-Aux
  if(mod$Vec_Bool[4]==F){elem.p1<-res[seq(a3+4,a4+4),]}else elem.p1<-Aux

  cat("---------------------------------------------------------------\n")
  cat(paste("-------------- Fixed effects for ",
            link[1], "(mu) --------------------\n", sep=''))
  cat("---------------------------------------------------------------\n")
  printCoefmat(elem.mu, P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("-------------- Fixed effects for ",
            link[2], "(sigma) --------------------\n", sep=''))
  cat("---------------------------------------------------------------\n")
  printCoefmat(elem.sigma, P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("-------------- Fixed effects for ",
            link[3], "(p0) --------------------\n", sep=''))
  cat("---------------------------------------------------------------\n")
  printCoefmat(elem.p0, P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("-------------- Fixed effects for ",
            link[4], "(p1) --------------------\n", sep=''))
  cat("---------------------------------------------------------------\n")
  printCoefmat(elem.p1, P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat("---------------------------------------------------------------\n")

}

# print.ZOIP<-function(mod){
#
#   mod$nparm[mod$Vec_Bool==T]<--1
#
#   a1<-cumsum(mod$nparm)[1]
#   a2<-cumsum(mod$nparm)[2]
#   a3<-cumsum(mod$nparm)[3]
#   a4<-cumsum(mod$nparm)[4]
#
#   Aux<-c(0)
#   names(Aux)<-c('(intercept)')
#
#   if(mod$Vec_Bool[1]==F){elem.mu<-mod$par[seq(1,a1+1)]}else elem.mu<-Aux
#   if(mod$Vec_Bool[2]==F){elem.sigma<-mod$par[seq(a1+2,a2+2)]}else elem.sigma<-Aux
#   if(mod$Vec_Bool[3]==F){elem.p0<-mod$par[seq(a2+3,a3+3)]}else elem.p0<-Aux
#   if(mod$Vec_Bool[4]==F){elem.p1<-mod$par[seq(a3+4,a4+4)]}else elem.p1<-Aux
#
#   cat("Call:\n")
#   print(mod$call)
#   cat("\n Results: \n")
#   cat("\n Estimated coefficients for g(mu): \n")
#   print(elem.mu)
#   cat("\n Estimated coefficients for g(sigma): \n")
#   print(elem.sigma)
#   cat("\n Estimated coefficients for g(p0): \n")
#   print(elem.p0)
#   cat("\n Estimated coefficients for g(p1): \n")
#   print(elem.p1)
#   cat("\n Convergence \n")
#   print(mod$Convergence)
#   cat("\n message \n")
#   print(mod$message)
#   cat("\n iterations \n")
#   print(mod$iterations)
#
# }
