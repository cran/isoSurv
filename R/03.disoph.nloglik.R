nloglik1=function(delta1,lambda,delta2,psi,wt2){
  lik=-sum(delta1*log(lambda))-sum(delta2*psi)+sum(wt2*exp(psi))
  return(lik)
}

nloglik2=function(Delta1,lambda,Delta2,wb,WT,n){
  lik1=-sum(Delta1*log(lambda),na.rm=TRUE)
  lik3=-sum(Delta2*wb,na.rm=TRUE)
  lik4=rep(NA,n)
  for(i in 1:n)
    lik4[i]=sum(WT[i,]*exp(wb))*lambda[i]
  lik4=sum(lik4)

  #lik4=0 #above is the same as here (useful for C++)
  #for(i in 1:n)
  #  for(j in 1:n)
  #    lik4=lik4+WT[i,j]*exp(psi.upd[j]+wb.upd[j])*lambda.upd[i]

  lik=lik1+lik3+lik4
  return(lik)
}

nloglik3=function(Delta1,lambda,Delta2,psi,wb,WT,n){
  lik1=-sum(Delta1*log(lambda),na.rm=TRUE)
  lik2=-sum(Delta2*psi,na.rm=TRUE)
  lik3=-sum(Delta2*wb,na.rm=TRUE)
  lik4=rep(NA,n)
  for(i in 1:n)
    lik4[i]=sum(WT[i,]*exp(psi+wb))*lambda[i]
  lik4=sum(lik4)
  lik.upd=lik1+lik2+lik3+lik4
}
