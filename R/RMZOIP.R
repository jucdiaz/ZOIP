#' Regresion ZOIP
#'
#' la funcion RM.ZOIP ajusta un modelo de regresion ZOIP via maxima verosimilitud. el modelo puede incluir o no
#' covariables en cualquiera de sus parametros, asi como puede ser un modelo inflado bilateralmente, unilateralmente o sin parametros
#' inflados
#'
#'
#' @usage RM.ZOIP(formula.mu,formula.sigma=~1,formula.p0=~1,formula.p1=~1,data,link=c('identity','identity','identity','identity'),family='R-S')
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

RM.ZOIP<-function(formula.mu,formula.sigma=~1,formula.p0=~1,formula.p1=~1,data,link=c('identity','identity','identity','identity'),family='R-S'){
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

  fg<-Formula.ZOIP(formula.mu,formula.sigma,formula.p0,formula.p1,link,family)
  opt<-fit.ZOIP(formula.mu,formula.sigma,formula.p0,formula.p1,data,link,family,fg)

  nparm=c(nparm.mu+1,nparm.sigma+1,nparm.p0+1,nparm.p1+1)
  names<-c('(intercept)',var.mu.p,'(intercept)',var.sigma.p,'(intercept)',var.p0.p,'(intercept)',var.p1.p)

  names(opt$par)<-names
  Pos_inter<-c(1,nparm.mu+2,nparm.mu+nparm.sigma+3,nparm.mu+nparm.sigma+nparm.p0+4)
  i=1
  while(i<=4){
    if((nparm-1)[i]==0 && opt$par[Pos_inter[i]]<0.0001){
      val<-paste0('X[',Pos_inter[i],']')
      fg<-gsub(val,'1e-16',fg,fixed = TRUE)
    }
    i=i+1
  }

  HM<-numDeriv::hessian(func=ll.ZOIP,x=opt$par,y=Yi,fg=fg,data=data,family=family)

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


  Result<-list(par=NULL,Convergence=NULL,message=NULL,iterations=NULL,HM=NULL,nparm=NULL,Vec_Bool=NULL)

  Result[[1]]<-opt$par
  Result[[2]]<-opt$convergence
  Result[[3]]<-opt$message
  Result[[4]]<-opt$iterations
  Result[[5]]<-HM
  Result[[6]]<-c(nparm.mu,nparm.sigma,nparm.p0,nparm.p1)
  Result[[7]]<-Vec_Bool


  Result$call <- match.call()
  class(Result)<-'ZOIP'

  return(Result)
}

Formula.ZOIP<-function(formula.mu,formula.sigma,formula.p0,formula.p1,link,family){

  var.mu.p<-attr(terms(formula.mu),'term.labels')
  var.sigma.p<-attr(terms(formula.sigma),'term.labels')
  var.p0.p<-attr(terms(formula.p0),'term.labels')
  var.p1.p<-attr(terms(formula.p1),'term.labels')

  nparm.mu<-length(var.mu.p)
  nparm.sigma<-length(var.sigma.p)
  nparm.p0<-length(var.p0.p)
  nparm.p1<-length(var.p1.p)

  if(nparm.mu!=0){
    var.mu<-paste0('data$',var.mu.p)
  }else{
    var.mu<-var.mu.p
  }

  if(nparm.sigma!=0){
    var.sigma<-paste0('data$',var.sigma.p)
  }else{
    var.sigma<-var.sigma.p
  }

  if(nparm.p0!=0){
    var.p0<-paste0('data$',var.p0.p)
  }else{
    var.p0<-var.p0.p
  }

  if(nparm.p1!=0){
    var.p1<-paste0('data$',var.p1.p)
  }else{
    var.p1<-var.p1.p
  }


  X.mu<-paste0('X[',seq(1,nparm.mu+1),']')
  X.sigma<-paste0('X[',seq(nparm.mu+2,nparm.mu+nparm.sigma+2),']')
  X.p0<-paste0('X[',seq(nparm.mu+nparm.sigma+3,nparm.mu+nparm.sigma+nparm.p0+3),']')
  X.p1<-paste0('X[',seq(nparm.mu+nparm.sigma+nparm.p0+4,nparm.mu+nparm.sigma+nparm.p0+nparm.p1+4),']')

  F.mu<-ifelse(nparm.mu==0,X.mu,paste0(X.mu[1],
                                       paste(paste0('+',paste0(X.mu[-1],paste0('*',var.mu))),collapse='')))
  F.sigma<-ifelse(nparm.sigma==0,X.sigma,
                  paste0(X.sigma[1],paste(paste0('+',paste0(X.sigma[-1],paste0('*',var.sigma))),collapse='')))
  F.p0<-ifelse(nparm.p0==0,X.p0,
               paste0(X.p0[1],paste(paste0('+',paste0(X.p0[-1],paste0('*',var.p0))),collapse='')))
  F.p1<-ifelse(nparm.p1==0,X.p1,
               paste0(X.p1[1],paste(paste0('+',paste0(X.p1[-1],paste0('*',var.p1))),collapse='')))

  F.mu.2<-if(link[1]=='logit'){
    paste0('inv.logit(',F.mu,')')
  }else if(link[1]=='exp'){
    paste0('exp(',F.mu,')')
  }else F.mu

  F.sigma.2<-if(link[2]=='logit'){
    paste0('inv.logit(',F.sigma,')')
  }else if(link[2]=='exp'){
    paste0('exp(',F.sigma,')')
  }else F.sigma

  F.p0.2<-ifelse(link[3]=='logit',paste0('inv.logit(',F.p0,')'),F.p0)
  F.p1.2<-ifelse(link[4]=='logit',paste0('inv.logit(',F.p1,')'),F.p1)

  F_final<-paste0('-1*sum(dZOIP(y,mu=',F.mu.2,',sigma=',F.sigma.2,
                  ',p0=',F.p0.2,',p1=',F.p1.2,',family,log = TRUE))')

  return(F_final)
}


fit.ZOIP<-function(formula.mu,formula.sigma,formula.p0,formula.p1,data,link,family,fg){
  var.mu.p<-attr(terms(formula.mu),'term.labels')
  var.sigma.p<-attr(terms(formula.sigma),'term.labels')
  var.p0.p<-attr(terms(formula.p0),'term.labels')
  var.p1.p<-attr(terms(formula.p1),'term.labels')

  nparm.mu<-length(var.mu.p)
  nparm.sigma<-length(var.sigma.p)
  nparm.p0<-length(var.p0.p)
  nparm.p1<-length(var.p1.p)

  val.inic<-rep(0.1,nparm.mu+nparm.sigma+nparm.p0+nparm.p1+4)


  lower.val=c(rep(ifelse((link[1]=='logit'|| link[1]=='exp'),-Inf,1e-16),nparm.mu+1),rep(ifelse((link[2]=='logit' || link[2]=='exp'),-Inf,1e-16),nparm.sigma+1),
              rep(ifelse(link[3]=='logit',-Inf,0),nparm.p0+1),rep(ifelse(link[4]=='logit',-Inf,0),nparm.p1+1))

  upper.mu<-if(link[1]=='logit' || ((link[1]=='exp' || link[1]=='identity') && family=='Original')){
    Inf
  }else 0.999999999

  upper.sigma<-if(link[2]=='logit' || link[2]=='exp' || (link[2]=='identity' && family!='R-S')){
    Inf
  }else 0.999999999

  upper.val=c(rep(upper.mu,nparm.mu+1),rep(upper.sigma,nparm.sigma+1),
              rep(ifelse(link[3]=='logit',Inf,1),nparm.p0+1),rep(ifelse(link[4]=='logit',Inf,1),nparm.p1+1))

  Yi<-paste0('data$',rownames(attr(terms(formula.mu),'factors'))[1])
  Yi<-eval(parse(text=Yi))

  opt<-nlminb(val.inic,objective=ll.ZOIP,y=Yi,fg=fg,data=data,family=family,lower=lower.val,upper=upper.val)

  return(opt)

}

ll.ZOIP<-function(X,y,fg,data,family){
  eval(parse(text=fg))
}


