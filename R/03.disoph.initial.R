disoph.initial=function(formula,data.all,P,Q,shape1,shape2,z2.obs,Zk,x1.obs){
  oldwarn <- options(warn=2)
  on.exit(options(oldwarn))

  #beta
  formula1=stats::as.formula(formula)
  cox.warning=try(cox.fit<-survival::coxph(formula1,data=data.all), silent=TRUE)

  if(inherits(cox.warning,"try-error")){
    data.all$DELTA=1
    cox.warning=try(cox.fit<-survival::coxph(formula1,data=data.all))
  }

  beta.hat=cox.fit$coefficients
  if(P>0){
    beta=matrix(beta.hat[-P],ncol=1) # first beta for isotonic cov
  }else{
    beta=matrix(beta.hat,ncol=1) #no iso cov
  }

  #psi
  psi=NA
  if(P>0){
    if(!is.null(shape2)){
      if(shape2=='increasing'){
        psi= abs(beta.hat[1])*(z2.obs-Zk)
      }else if(shape2=='decreasing'){
        psi=-abs(beta.hat[1])*(z2.obs-Zk)
      }
    }
  }

  #lambda
  wb.fit=survival::survreg(survival::Surv(X,DELTA)~1,dist='weibull',data=data.all)
  a=1/wb.fit$scale  #wb.shape
  b=exp(wb.fit$coefficients) #wb.scale
  #a>1 increasing
  #a=1 constant
  #a<1 decreasing
  a=abs(a)
  if(a==1) a+runif(1,0,0.01)
  lambda=a^(-b)*b*x1.obs^(a-1)

  return(list(psi=psi,beta=beta,lambda=lambda))
}
