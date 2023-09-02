#par & iso cov only
disoph.ti3=function(TIME,STATUS,Z,ZNAME,P,Q,shape1,shape2,K, maxiter, eps){
  #Time: observed survival time
  #Status: 1 for event I(X<=T); 0 for censored
  #Z: Z1=Z[,1], ... ZP=Z[,P] for isotonic covariates; W1=Z[,P+1], W1=Z[,P+Q]: additional covariates
  #ZNAME: ZNAME[1:P] for name of isotonic covariates; ZNAME[P+1;P+Q] for name of additional covariates
  #shape1: increasing or decreasing for baseline hazard
  #shape2: increasing or decreasing for isotonic covariate
  #K: anchor
  #maxiter: max number of iterations
  #eps: convergence criteria

  #0.data set-up
  data.all=data.frame(X=TIME,DELTA=STATUS,Z)
  n=nrow(data.all)

  #1. time
  #1.1. data sort
  data1=data.all[order(data.all$X),] #sorted by time
  data1$Dt=diff(c(0,data1$X))
  Dt1=data1$Dt
  X1=data1$X; #sorted by X
  Delta1=data1$DELTA

  #1.2. reduced likelihood
  x1.obs=unique(X1[Delta1==1]) #X1 is already sorted
  n1=length(x1.obs) #n1=n2 if no tie

  #1.3. interval
  intv1=list();
  if(shape1=="increasing"){
    for(i in 1:(n1-1))
      intv1[[i]]=c(x1.obs[i],x1.obs[i+1])
    intv1[[n1]]=c(x1.obs[n1],Inf)
  }else if(shape1=="decreasing"){
    intv1[[1]]=c(-Inf,x1.obs[1])
    for(i in 2:n1)
      intv1[[i]]=c(x1.obs[i-1],x1.obs[i])
  }

  #1.4. delta1
  delta1=rep(NA,n1) #compute delta1 & delta2 & wt
  r1=list() #delta1.s, all 1's if no ties
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

  #2. iso covariate
  #2.1. data sort
  data2=data1[order(data1$Z),]
  X2=data2$X  #sorted by Z1
  Z2=data2$Z1;
  Delta2=data2$DELTA

  #2.2. reduced likelihood
  z2.obs=unique(Z2[Delta2==1]) #z2 is sorted
  n2=length(z2.obs)

  #2.3. interval
  intv2=list()
  if(shape2=="increasing"){
    for(j in 1:(n2-1))
      intv2[[j]]=c(z2.obs[j],z2.obs[j+1])
    intv2[[n2]]=c(z2.obs[n2],Inf)
  }else if(shape2=="decreasing"){
    intv2[[1]]=c(-Inf,z2.obs[1])
    for(j in 2:n2)
      intv2[[j]]=c(z2.obs[j-1],z2.obs[j])
  }

  #2.4. delta2
  delta2=rep(NA,n2) #compute delta1 & delta2 & wt
  r2=list() #delta2.s, all 1's if no ties
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

  #3. Data for parametric cov
  W=as.matrix(data2[,colnames(data2)%in%paste0("W",1:Q)]) #W from data2, i.e. same order

  #4. WT & wt
  WT=matrix(NA,n,n)
  for(i in 1:n)
    for(j in 1:n)
      WT[i,j]=(X1[i]<=X2[j])*Dt1[i]

  #5.anchor & initial
  #5.1. anchor
  k2=which.max(z2.obs[z2.obs<=K])
  if(k2==0) k2=1 #k2 is set to 1 if min(z.obs) < K
  Zk=z2.obs[k2]
  data.all$Z.BAR=data.all$Z1-Zk
  formula="survival::Surv(X,DELTA)~Z.BAR"

  #5.2. coxph
  formula=paste(c(formula,paste0("+W",1:Q)),collapse ="")
  res.initial=disoph.initial(formula=formula,data.all=data.all,P=P,Q=Q,shape1=shape1,shape2=shape2,z2.obs=Z2,Zk=Zk,x1.obs=X1)
    #z2 & x1 for full rank

  psi=res.initial$psi
  beta=res.initial$beta
  lambda=res.initial$lambda
  wb=W%*%beta
  lik=nloglik3(Delta1,lambda,Delta2,psi,wb,WT,n)

  #6. back and forth algo
  conv=iter=0
  dist=1

  wt1.full=wt2.full=rep(NA,n) #full
  wt1.red=rep(NA,n1); wt2.red=rep(NA,n2) #reduced

  while(dist>eps){

    iter=iter+1
    if(iter==maxiter)  break

    #labmda
    for(i in 1:n)
      wt1.full[i]=sum(WT[i,]*exp(psi)*exp(wb))
    for(i in 1:n1)
      wt1.red[i]=sum(wt1.full[r1[[i]]])

    if(shape1=="increasing"){
      lambda.upd=Iso::pava(delta1/wt1.red,wt1.red) #delta 1 was reduced above along with r1
    }else if(shape1=="decreasing"){
      lambda.upd=-Iso::pava(-delta1/wt1.red,wt1.red)
    }
    lambda.upd=disoph.BTF1(n,n1,lambda.upd,X1,x1.obs,shape1)

    #psi
    for(j in 1:n)
      wt2.full[j]=sum(WT[,j]*lambda.upd)*exp(wb[j]) #delta 2 was reduced above along with r1
    for(j in 1:n2)
      wt2.red[j]=sum(wt2.full[r2[[j]]])

    if(shape2=="increasing"){
      psi.upd=log(Iso::pava(delta2/wt2.red,wt2.red))
    }else if(shape2=="decreasing"){
      psi.upd=log(-Iso::pava(-delta2/wt2.red,wt2.red))
    }
    psi.upd=psi.upd-psi.upd[k2] #anchor
    psi.upd=disoph.BTF2(n,n2,psi.upd,Z2,z2.obs,shape2)

    #beta
    disoph.NR=disoph.NR3(lambda.upd,psi.upd,beta,wb,W,Delta2,Q,n,WT,eps)
    beta.upd=disoph.NR$beta
    wb.upd=disoph.NR$wb

    #6.3. lik & dist
    lik.upd=nloglik3(Delta1,lambda.upd,Delta2,psi.upd,wb.upd,WT,n)
    dist=abs((lik.upd-lik)/lik)

    lik=lik.upd
    psi=psi.upd
    lambda=lambda.upd
    beta=beta.upd
    wb=wb.upd
  }

  conv="not converged"
  if(iter<maxiter)
    conv="converged"

  #7 return
  labmda.res=data.frame(t=X1,lambda.hat=lambda.upd)
  psi.res=data.frame(z=Z2,psi.hat=psi.upd)
  colnames(psi.res)=c(ZNAME[1],"psi.hat")

  res=list(iso.bh=labmda.res,
           iso.cov= psi.res,
           conv=conv,
           iter=iter,
           Zk=Zk,
           shape.bh=shape1,
           shape.cov=shape2,
           beta=beta.upd
           )

  return(res)
}
