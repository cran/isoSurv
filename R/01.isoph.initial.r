isoph.initial=function(formula,data.all,Q,shape,z.obs,Zk){
  oldwarn <- options(warn=2)
  on.exit(options(oldwarn))

  cox.warning=try(cox.fit<-coxph(as.formula(formula),data=data.all), silent=TRUE)
  if(class(cox.warning)=="try-error"){
    beta.hat=mean(data.all$DELTA)/mean(data.all$Z.BAR)
    beta=matrix(c(0,Q),ncol=1);
  }else{
    beta.hat=cox.fit$coefficients
    beta=matrix(beta.hat[-1],ncol=1)
  }
  if(shape=='increasing'){
    psi= abs(beta.hat[1])*(z.obs-Zk)
  }else if(shape=='decreasing'){
    psi=-abs(beta.hat[1])*(z.obs-Zk)
  }

  return(list(psi=psi,beta=beta))
}
