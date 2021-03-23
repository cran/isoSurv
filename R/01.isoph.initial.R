isoph.initial=function(formula,data.all,Q,shape,z.obs,Zk){
  oldwarn <- options(warn=2)
  on.exit(options(oldwarn))

  formula1=stats::as.formula(formula)
  cox.warning=try(cox.fit<-survival::coxph(formula1,data=data.all), silent=TRUE)

  if(class(cox.warning)=="try-error"){
    data.all$DELTA=1
    cox.warning=try(cox.fit<-survival::coxph(formula1,data=data.all))
  }
  
  beta.hat=cox.fit$coefficients
  beta=matrix(beta.hat[-1],ncol=1)
  
  if(shape=='increasing'){
    psi= abs(beta.hat[1])*(z.obs-Zk)
  }else if(shape=='decreasing'){
    psi=-abs(beta.hat[1])*(z.obs-Zk)
  }

  return(list(psi=psi,beta=beta))
}
