isoph.initial=function(formula,data.all,P,Q,shape,z.obs,Zk){
  oldwarn <- options(warn=2)
  on.exit(options(oldwarn))

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

  psi=NA
  if(P>0){
    if(shape=='increasing'){
      psi= abs(beta.hat[1])*(z.obs-Zk)
    }else if(shape=='decreasing'){
      psi=-abs(beta.hat[1])*(z.obs-Zk)
    }
  }

  return(list(psi=psi,beta=beta))
}
