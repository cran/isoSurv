disoph.ti=function(TIME,STATUS,Z,ZNAME,P,Q,shape1,shape2,K, maxiter, eps){
  #Time: observed survival time
  #Status: 1 for event I(X<=T); 0 for censored
  #Z: Z1=Z[,1], ... ZP=Z[,P] for isotonic covariates; W1=Z[,P+1], W1=Z[,P+Q]: additional covariates
  #ZNAME: ZNAME[1:P] for name of isotonic covariates; ZNAME[P+1;P+Q] for name of additional covariates
  #shape1: increasing or decreasing for baseline hazard
  #shape2: increasing or decreasing for isotonic covariate
  #K: anchor
  #maxiter: max number of iterations
  #eps: convergence criteria

  #1.data set-up
  data.all=data.frame(X=TIME,DELTA=STATUS,Z)
  data1=data.all[order(data.all$X),] #sorted by time
  data1$Dt=diff(c(0,data1$X))
  data2=data1[order(data1$Z1),]
  n=nrow(data1)
  
  #2. likelihood, using upper character
  X1=data1$X; #sorted by X
  Delta1=data1$DELTA
  X2=data2$X  #sorted by Z1
  Z2=data2$Z1;
  Delta2=data2$DELTA
  Dt1=data1$Dt
  WT=matrix(NA,n,n)
  for(i in 1:n)
    for(j in 1:n)
      WT[i,j]=(X1[i]<=X2[j])*Dt1[i]

  #3. reduced likelihood, using lower character
  x1.obs=unique(X1[Delta1==1]) #X1 & Z2 are already sorted
  z2.obs=unique(Z2[Delta2==1])
  n1=length(x1.obs) #n1=n2 if no tie
  n2=length(z2.obs)
  intv1=intv2=list();

  if(shape1=="increasing"){
    for(i in 1:(n1-1))
      intv1[[i]]=c(x1.obs[i],x1.obs[i+1])
    intv1[[n1]]=c(x1.obs[n1],Inf)
  }else if(shape1=="decreasing"){
    intv1[[1]]=c(-Inf,x1.obs[1])
    for(i in 2:n1)
      intv1[[i]]=c(x1.obs[i-1],x1.obs[i])
  }
  if(shape2=="increasing"){
    for(j in 1:(n2-1))
      intv2[[j]]=c(z2.obs[j],z2.obs[j+1])
    intv2[[n2]]=c(z2.obs[n2],Inf)
  }else if(shape2=="decreasing"){
    intv2[[1]]=c(-Inf,z2.obs[1])
    for(j in 2:n2)
      intv2[[j]]=c(z2.obs[j-1],z2.obs[j])
  }

  delta1=rep(NA,n1) #compute delta1 & delta2 & wt
  delta2=rep(NA,n2)
  wt=matrix(NA,n1,n2)
  r1=r2=list() #delta1.s & delta2.s are all 1s if no ties
  if(shape1=="increasing"){ #left cont
    for(i in 1:n1){
      r1[[i]]=intv1[[i]][1]<= X1 & X1 < intv1[[i]][2]
      delta1[i]=sum(Delta1[r1[[i]]])
    }
  }else if(shape1=="decreasing"){ #right cont
    for(i in 1:n1){
      r1[[i]]=intv1[[i]][1]< X1 & X1 <= intv1[[i]][2]
      delta1[i]=sum(Delta1[r1[[i]]])
    }
  }

  if(shape2=="increasing"){ #left cont
    for(j in 1:n2){
      r2[[j]]=intv2[[j]][1]<= Z2 & Z2 < intv2[[j]][2]
      delta2[j]=sum(Delta2[r2[[j]]])
    }
  }else if(shape2=="decreasing"){ #right cont
    for(j in 1:n2){
      r2[[j]]=intv2[[j]][1]< Z2 & Z2 <= intv2[[j]][2]
      delta2[j]=sum(Delta2[r2[[j]]])
    }
  }

  for(i in 1:n1)
    for(j in 1:n2)
      wt[i,j]=sum(WT[r1[[i]],r2[[j]]])

  #4.anchor
  k2=which.max(z2.obs[z2.obs<=K])
  if(k2==0) k2=1 #k2 is set to 1 if min(z.obs) < K
  Zk=z2.obs[k2]

  #5. initial values
  data.all$Z.BAR=data.all$Z1-Zk
  formula="survival::Surv(X,DELTA)~Z.BAR"
  if(Q>0) formula=paste(c(formula,paste0("+W",1:Q)),collapse ="")

  res.initial=isoph.initial(formula=formula,data.all=data.all,Q=Q,shape=shape2,z.obs=z2.obs,Zk=Zk)
  psi=res.initial$psi
  beta=res.initial$beta

  #6. back and forth algo
  conv=iter=0
  wt1=rep(NA,n1)
  wt2=rep(NA,n2)
  lik=0
  dist=1
  
  while(dist>eps){
    iter=iter+1
    if(iter==maxiter)  break

    #6.1. update lambda
    for(i in 1:n1)
      wt1[i]=sum(wt[i,]*exp(psi))
    
    if(shape1=="increasing"){
      lambda.upd=Iso::pava(delta1/wt1,wt1)
    }else if(shape1=="decreasing"){
      lambda.upd=-Iso::pava(-delta1/wt1,wt1)
    }
    
    #6.2. update psi
    for(j in 1:n2)
      wt2[j]=sum(wt[,j]*lambda.upd)
    if(shape2=="increasing"){
      psi.upd=log(Iso::pava(delta2/wt2,wt2))
    }else if(shape2=="decreasing"){
      psi.upd=log(-Iso::pava(-delta2/wt2,wt2))
    }
    psi.upd=psi.upd-psi.upd[k2] #anchor

    #6.3. lik & dist
    lik.upd=-sum(delta1*log(lambda.upd))-sum(delta2*psi.upd)+sum(wt2*exp(psi.upd))
    dist=abs((lik.upd-lik)/lik)

    lik=lik.upd
    psi=psi.upd
    lambda=lambda.upd
  }
  conv="not converged"
  if(iter<maxiter)
    conv="converged"

  #6. BTF (back to the full rank)
  lambda.full=disoph.BTF1(n,n1,lambda.upd,X1,x1.obs,shape1)
  psi.full=disoph.BTF2(n,n2,psi.upd,Z2,z2.obs,shape2)

  #7 return
  labmda.res=data.frame(time=X1,lambda.hat=lambda.full)
  psi.res=data.frame(z=Z2,psi.hat=psi.full)
  colnames(psi.res)=c(ZNAME[1],"psi.hat")

  res=list(iso.bh=labmda.res,
           iso.cov= psi.res,
           conv=conv,iter=iter,Zk=Zk,
           shape.bh=shape1,shape.cov=shape2)

  #res=list(x=X1,               x.obs=X1[which(Delta1==1)],
  #         lambda=lambda.full, lambda.unq=lambda.full[which(Delta1==1)],
  #         z=Z2,               z.unq=Z2[which(Delta2==1)],
  #         psi=psi.full,       psi.obs=psi.full[which(Delta2==1)],
  #         conv=conv,iter=iter,Zk=Zk)
  return(res)
}
