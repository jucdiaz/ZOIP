#' Regresion ZOIP
#'
#' la funcion RM.ZOIP ajusta un modelo de regresion ZOIP via maxima verosimilitud. el modelo puede incluir o no
#' covariables en cualquiera de sus parametros, asi como puede ser un modelo inflado bilateralmente, unilateralmente o sin parametros
#' inflados
#'
#'
#' @usage RM.ZOIP2(formula.mu,formula.sigma=~1,formula.p0=~1,formula.p1=~1,data,link=c('identity','identity','identity','identity'),family='R-S')
#' @param formula.mu Formula que define la funcion de regresion para mu, p.e y~x1+x2, es necesario definir la variable respuesta.
#' @param formula.sigma Formula que define la funcion de regresion para sigma, p.e ~x1, es necesario definir la variable respuesta.
#' @param formula.p0 Formula que define la funcion de regresion para p0, p.e ~x1, es necesario definir la variable respuesta.
#' @param formula.p1 Formula que define la funcion de regresion para p1, p.e ~x1, es necesario definir la variable respuesta.
#' @param data datos en fomato data.frame donde debe contener las nombres de las columnas tal cual como estan en las formulas.
#' @param link se debe escoger la funcion de enlace adecuada a cada parametro a estimar de acuerdo a las opciones de familia y formula, ver detalles y vignettes.
#' @param family eleccion de la parametrizacion o distribucion deseada, family='R-S' parametrizacion distribucion beta Rigby y Stasinopoulos, 'F-C' distribucion Beta parametrizacion Ferrari y Cribari-Neto, Original distribucion beta parametrizacion original, 'Simplex' distribucion simplex.
#' @examples
#'
#' #Test 1--------------------------------------------------
#' library(boot)
#' library(numDeriv)
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
#' system.time(y_i<-apply(param,1,function(x){rZOIP(1,mu=x[1],sigma=x[2],p0=x[3],p1=x[4],family='Original')}))
#' data<-as.data.frame(cbind(y_i,x1,x2))
#'
#' formula.mu=y_i~x1
#' formula.sigma=~x1+x2
#' formula.p0=~1
#' formula.p1=~x2
#' link=c('exp','exp','identity','logit')
#' family='Original'
#' mod<-RM.ZOIP(formula.mu=formula.mu,formula.sigma=formula.sigma,formula.p0=formula.p0,formula.p1=formula.p1,data=data,link=link,family=family)
#' mod
#' summary(mod)
#'
#' #Test 2--------------------------------------------------
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
#' d1<-0
#' p1_i<-rep(d1,length(n))
#'
#'
#' param<-cbind(mu_i,sigma_i,p0_i,p1_i)
#'
#' system.time(y_i<-apply(param,1,function(x){rZOIP(1,mu=x[1],sigma=x[2],p0=x[3],p1=x[4],family='Original')}))
#' data<-as.data.frame(cbind(y_i,x1,x2))
#'
#' formula.mu=y_i~x1
#' formula.sigma=~x1+x2
#' formula.p0=~1
#' formula.p1=~1
#' link=c('exp','exp','identity','identity')
#' family='Original'
#' mod<-RM.ZOIP(formula.mu=formula.mu,formula.sigma=formula.sigma,formula.p0=formula.p0,formula.p1=formula.p1,data=data,link=link,family=family)
#' mod
#' summary(mod)
#'
#' @export


RM.ZOIP2<-function(formula.mu,formula.sigma=~1,formula.p0=~1,formula.p1=~1,data,link=c('identity','identity','identity','identity'),family='R-S',optimizer='nlminb'){
  library(boot)
  library(numDeriv)
  if (any(family != 'R-S') && any(family != 'F-C') && any(family != 'Original') && any(family != 'Simplex'))
    stop(paste("family must be in R-S, F-C, Original, Simplex", "\n", ""))
  if(any(length(as.character(attr(terms(formula.mu),'variable')))==1))
    stop(paste("formula.mu must have response variable","\n",""))
  if(any(attr(terms(formula.sigma),'response')>=1))
    stop(paste("sigma can't have response variable in formula.sigma","\n",""))
  if(any(attr(terms(formula.p0),'response')>=1))
    stop(paste("p0 can't have response variable in formula.p0","\n",""))
  if(any(attr(terms(formula.p1),'response')>=1))
    stop(paste("p1 can't have response variable in formula.p1","\n",""))
  if(any(dim(data)[2])<1)
    stop(paste("Object 'data' must have data"),"\n","")
  if(any(is.data.frame(data)==F))
    stop(paste("Object 'data' must be a data frame class","\n",""))
  if(any(length(link)<4))
    stop(paste("link must have four enlace functions", "\n",""))
  if(any(link[1]!='identity' && link[1]!='logit' && link[1]!='exp'))
    stop(paste("link for mu must be in identity, logit, exp", "\n",""))
  if(any(link[2]!='identity' && link[2]!='logit' && link[2]!='exp'))
    stop(paste("link for sigma must be in identity, logit, exp", "\n",""))
  if(any(link[3]!='identity' && link[3]!='logit'))
    stop(paste("link for p0 must be in identity, logit", "\n",""))
  if(any(link[4]!='identity' && link[4]!='logit'))
    stop(paste("link for p1 must be in identity, logit", "\n",""))
  if(any(length(as.character(attr(terms(formula.mu),'variable')))==2 && link[1]!='identity'))
    stop(paste("mu don't have covariables then link must be identity", "\n",""))
  if(any(family!='Original' && length(as.character(attr(terms(formula.mu),'variable')))>2 && link[1]!='logit'))
    stop(paste("If family is diferent a Original and mu have covariables then link must be logit", "\n",""))
  if(any(family=='Original'&& length(as.character(attr(terms(formula.mu),'variable')))>2 && link[1]!='exp'))
    stop(paste("If family is Original and mu have covariables then link must be exp", "\n",""))
  if(any(length(as.character(attr(terms(formula.sigma),'variable')))==1 && link[2]!='identity'))
    stop(paste("sigma don't have covariables then link must be identity", "\n",""))
  if(any(family!='R-S' && length(as.character(attr(terms(formula.sigma),'variable')))>1 && link[2]!='exp'))
    stop(paste("If family is diferent a R-S and sigma have covariables then link must be exp", "\n",""))
  if(any(family=='R-S'&& length(as.character(attr(terms(formula.sigma),'variable')))>1 && link[2]!='logit'))
    stop(paste("If family is R-S and sigma have covariables then link must be logit", "\n",""))
  if(any(length(as.character(attr(terms(formula.p0),'variable')))==1 && link[3]!='identity'))
    stop(paste("p0 don't have covariables then link must be identity", "\n",""))
  if(any(length(as.character(attr(terms(formula.p0),'variable')))>1 && link[3]!='logit'))
    stop(paste("p0 have covariables then link must be logit", "\n",""))
  if(any(length(as.character(attr(terms(formula.p1),'variable')))==1 && link[4]!='identity'))
    stop(paste("p1 don't have covariables then link must be identity", "\n",""))
  if(any(length(as.character(attr(terms(formula.p1),'variable')))>1 && link[4]!='logit'))
    stop(paste("p1 have covariables then link must be logit", "\n",""))

  var.mu.p<-attr(terms(formula.mu),'term.labels')
  var.sigma.p<-attr(terms(formula.sigma),'term.labels')
  var.p0.p<-attr(terms(formula.p0),'term.labels')
  var.p1.p<-attr(terms(formula.p1),'term.labels')

  nparm.mu<-length(var.mu.p)
  nparm.sigma<-length(var.sigma.p)
  nparm.p0<-length(var.p0.p)
  nparm.p1<-length(var.p1.p)

  matri<-model.matrix.ZOIP(formula.mu,formula.sigma,formula.p0,formula.p1,data=data)
  opt<-fit.ZOIP2(matri,link,family,optimizer)

  nparm=c(nparm.mu+1,nparm.sigma+1,nparm.p0+1,nparm.p1+1)
  names<-c(colnames(matri$mat.mu),colnames(matri$mat.sigma),colnames(matri$mat.p0),colnames(matri$mat.p1))
  names(opt$par)<-names


  Pos_inter<-c(1,nparm.mu+2,nparm.mu+nparm.sigma+3,nparm.mu+nparm.sigma+nparm.p0+4)

  X.mu <- matri$mat.mu
  X.sigma <- matri$mat.sigma
  X.p0 <- matri$mat.p0
  X.p1 <- matri$mat.p1
  y <- matri$y

  HM<-numDeriv::hessian(func=ll.ZOIP2,x=opt$par,y=y,method='Richardson', X.mu=X.mu, X.sigma=X.sigma,X.p0=X.p0,X.p1=X.p1,link=link,family=family)

  i=1
  b=0
  Vec_Bool<-NULL
  while(i<=4){
    Vec_Bool[i]<-(nparm-1)[i]==0 && opt$par[Pos_inter[i]-b]<0.0001
    if((nparm-1)[i]==0 && opt$par[Pos_inter[i]-b]<0.0001){
      opt$par<-opt$par[-Pos_inter[i]+b]
      HM<-HM[-Pos_inter[i]+b,-Pos_inter[i]+b]
      b=b+1
    }
    i=i+1
  }


  Result<-list(par=NULL,Convergence=NULL,message=NULL,iterations=NULL,HM=NULL,nparm=NULL,Vec_Bool=NULL,objective=NULL)

  Result[[1]]<-opt$par
  Result[[2]]<-opt$convergence
  Result[[3]]<-opt$message
  Result[[4]]<-ifelse(optimizer=='nlminb',opt$iterations,NA)
  Result[[5]]<-HM
  Result[[6]]<-c(nparm.mu,nparm.sigma,nparm.p0,nparm.p1)
  Result[[7]]<-Vec_Bool
  Result[[8]]<-opt$objective


  Result$call <- match.call()
  class(Result)<-'ZOIP'

  return(Result)
}
fit.ZOIP2<-function(matri,link,family,optimizer){
  nparm.mu <- ncol(matri$mat.mu)
  nparm.sigma <- ncol(matri$mat.sigma)
  nparm.p0 <- ncol(matri$mat.p0)
  nparm.p1 <- ncol(matri$mat.p1)

  X.mu <- matri$mat.mu
  X.sigma <- matri$mat.sigma
  X.p0 <- matri$mat.p0
  X.p1 <- matri$mat.p1
  y <- matri$y


  val.inic<-rep(0.1,nparm.mu+nparm.sigma+nparm.p0+nparm.p1)

  lower.val=c(rep(ifelse((link[1]=='logit'|| link[1]=='exp'),-Inf,1e-16),nparm.mu),rep(ifelse((link[2]=='logit' || link[2]=='exp'),-Inf,1e-16),nparm.sigma),
              rep(ifelse(link[3]=='logit',-Inf,1e-16),nparm.p0),rep(ifelse(link[4]=='logit',-Inf,1e-16),nparm.p1))

  upper.mu<-if(link[1]=='logit' || ((link[1]=='exp' || link[1]=='identity') && family=='Original')){
    Inf
  }else 0.999999999

  upper.sigma<-if(link[2]=='logit' || link[2]=='exp' || (link[2]=='identity' && family!='R-S')){
    Inf
  }else 0.999999999

  upper.val=c(rep(upper.mu,nparm.mu),rep(upper.sigma,nparm.sigma),
              rep(ifelse(link[3]=='logit',Inf,0.999999999),nparm.p0),rep(ifelse(link[4]=='logit',Inf,0.999999999),nparm.p1))

  if (optimizer == 'nlminb') {
    opt <- nlminb(start=val.inic, objective=ll.ZOIP2,
                  y=y, X.mu=X.mu, X.sigma=X.sigma,X.p0=X.p0,X.p1=X.p1,
                  link=link,family=family,lower=lower.val,upper=upper.val)
    opt$objective <- -opt$objective
  }

  if (optimizer == 'optim') {

    opt <- optim(par=val.inic, fn=ll.ZOIP2,
                 y=y, X.mu=X.mu, X.sigma=X.sigma,X.p0=X.p0,X.p1=X.p1,
                 link=link,family=family,lower=lower.val,upper=upper.val)
    opt$objective <- -opt$value
  }

  return(opt)

}

ll.ZOIP2<-function(theta,y,X.mu,X.sigma,X.p0,X.p1,link,family){
  betas.mu <- matrix(theta[1:ncol(X.mu)], ncol=1)
  betas.sigma <- matrix(theta[seq(ncol(X.mu)+1,ncol(X.mu)+ncol(X.sigma))], ncol=1)
  betas.p0 <- matrix(theta[seq(ncol(X.mu)+ncol(X.sigma)+1,ncol(X.mu)+ncol(X.sigma)+ncol(X.p0))], ncol=1)
  betas.p1 <- matrix(theta[seq(ncol(X.mu)+ncol(X.sigma)+ncol(X.p0)+1,ncol(X.mu)+ncol(X.sigma)+ncol(X.p0)+ncol(X.p1))], ncol=1)

  if(link[1]=='identity'){
    mu<-X.mu %*% betas.mu
  }else if(link[1]=='logit'){
    mu<-1 / (1 + exp(- X.mu %*% betas.mu))
  }else if(link[1]=='exp'){
    mu<-exp(X.mu %*% betas.mu)
  }

  if(link[2]=='identity'){
    sigma<-X.sigma %*% betas.sigma
  }else if(link[2]=='logit'){
    sigma<-1 / (1 + exp(- X.sigma %*% betas.sigma))
  }else if(link[2]=='exp'){
    sigma<-exp(X.sigma %*% betas.sigma)
  }

  if(link[3]=='identity'){
    p0<-X.p0 %*% betas.p0
  }else if(link[3]=='logit'){
    p0<-1 / (1 + exp(- X.p0 %*% betas.p0))
  }

  if(link[4]=='identity'){
    p1<-X.p1 %*% betas.p1
  }else if(link[4]=='logit'){
    p1<-1 / (1 + exp(- X.p1 %*% betas.p1))
  }

  ll<-sum(dZOIP(x=y,mu=mu,sigma=sigma,p0=p0,p1=p1,family=family,log=TRUE))
  -ll
}


model.matrix.ZOIP <- function(formula.mu,formula.sigma,formula.p0,formula.p1, data=NULL) {
  stopifnot (class(formula.mu) == 'formula')
  stopifnot (class(formula.sigma) == 'formula')
  stopifnot (class(formula.p0) == 'formula')
  stopifnot (class(formula.p1) == 'formula')
  response <- all.vars(formula.mu)[1]
  formula.sigma <- as.formula(paste(response, paste(as.character(formula.sigma),
                                                    collapse='')))
  formula.p0 <- as.formula(paste(response, paste(as.character(formula.p0),
                                                 collapse='')))
  formula.p1 <- as.formula(paste(response, paste(as.character(formula.p1),
                                                 collapse='')))
  mat.mu <- model.matrix(formula.mu, data)
  mat.sigma <- model.matrix(formula.sigma, data)
  mat.p0 <- model.matrix(formula.p0, data)
  mat.p1 <- model.matrix(formula.p1, data)
  y <- model.frame(formula.mu, data=data)[, 1]
  matri<-list(mat.mu=mat.mu, mat.sigma=mat.sigma,mat.p0=mat.p0,mat.p1=mat.p1, y=y)
  return(matri)
}

