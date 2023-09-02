isoph.ti=function(TIME, STATUS, Z, ZNAME, P, Q, shape, K, maxiter, eps){
  #Time: observed survival time
  #Status: 1 for event I(X<=T); 0 for censored
  #Z: Z1=Z[,1], ... ZP=Z[,P] for isotonic covariates; W1=Z[,P+1], W1=Z[,P+Q]: additional covariates
  #ZNAME: ZNAME[1:P] for name of isotonic covariates; ZNAME[P+1;P+Q] for name of additional covariates
  #shape: increasing or decreasing for isotonic covariate
  #K: anchor
  #maxiter: max number of iterations
  #eps: convergence criteria

  #1.sorted by z
  data.all=data.frame(X=TIME,DELTA=STATUS,Z)
  data=data.all[order(data.all$Z1),]
  n.total=nrow(data)

  #2.remove subjects whose cov is less than z^*_(1)
  DELTA.Z1=aggregate(data$DELTA,by=list(Z1=data$Z1),sum)
  if(shape=='increasing'){
    if(DELTA.Z1[1,2]==0){
      Z.star=DELTA.Z1[which(DELTA.Z1[,2]>0)[1],1]
      data=data[data$Z1>=Z.star,]
    }
  }else if(shape=='decreasing'){
    if(DELTA.Z1[nrow(DELTA.Z1),2]==0){
      Z.star=DELTA.Z1[max(which(DELTA.Z1[,2]>0)),1]
      data=data[data$Z1<=Z.star,]
    }
  }
  n=nrow(data)

  #3. redefine variables
  t=data$X
  delta=data$DELTA
  z=data$Z1
  w=data[,-c(1:3)]
  if(Q==1){; w=matrix(data[,-c(1:3)])
  }else if(Q>=2){; w=as.matrix(data[,-c(1:3)])
  }

  z.obs=unique(z[delta==1])
  m=length(z.obs)
  t.obs=sort(unique(t))
  nt=length(t.obs)

  #4. anchor
  k=sum(z.obs<K)
  if(k==0) k=1 #k is set to 1 if min(z.obs) < K
  Zk=z.obs[k]

  #5. counting process
  Y=dN=matrix(0,n,nt)       #row is the subj, col is the time corresponding z_(i);
  for(i in 1:n){
    rank.t=which(t[i]==t.obs)
    Y[i,][1:rank.t]=1
    if(delta[i]==1) dN[i,][rank.t]=1
  }

  #6. initial values
  data.all$Z.BAR=data.all$Z1-Zk
  formula="survival::Surv(X,DELTA)~Z.BAR"
  if(Q>0) formula=paste(c(formula,paste0("+W",1:Q)),collapse ="")
  res.initial=isoph.initial(formula,data.all,P,Q,shape,z.obs,Zk)
  psi=res.initial$psi
  beta=res.initial$beta

  #7. picm & beta for newton raphson algo
  if(Q==0){ #no additional covs
    #RPA (with interval), Y, dN
    rpa.Y=isoph.RPA.ti(n, nt, m, z, z.obs, Y, dN, shape)
    Y2=rpa.Y$Y2
    dN2=rpa.Y$dN2

    #picm
    dNsum=colSums(dN2)
    Delta=rowSums(dN2)

    dist=0; exp.beta=NA
    picm=isoph.picm(psi,m,z.obs,Zk,k, dN2,Y2,dNsum,Delta, eps,maxiter, shape)
    psi.new=picm$psi.new
    iter=picm$iter
    conv=picm$conv
  }else{ #additional covariates
    iter=0;  dist=1;   beta.new=rep(NA,Q)
    while(dist>=eps){
      iter=iter+1
      if(iter>maxiter) break

      #RPA (with interval), Y, dN
      Yest=matrix(NA,n,nt)
      for(j in 1:nt) Yest[,j]=Y[,j]*exp(w%*%beta)

      rpa.Y=isoph.RPA.ti(n, nt, m, z, z.obs, Yest, dN, shape)
      Y2=rpa.Y$Y2
      dN2=rpa.Y$dN2

      #picm
      dNsum=colSums(dN2)
      Delta=rowSums(dN2)

      #estimate psi
      picm=isoph.picm(psi,m,z.obs,Zk,k, dN2,Y2,dNsum,Delta, eps,maxiter, shape)
      psi.new=picm$psi.new

      psi.full=isoph.BTF(m, n, z,z.obs, psi.new,shape)
      if(picm$conv==0) stop

      #estimate beta (Y1&w1 or Y2&w2 should be the same);
      beta.new=isoph.NR(w,beta,Q,psi.full,n,nt,Y,dN, maxiter,eps)

      #update;
      dist=sqrt(sum((psi.new-psi)^2))+sqrt(sum((beta.new-beta)^2))
      #this can be reduced to: (1-picm$conv)+sqrt(sum((beta.new-beta)^2))

      psi=psi.new
      beta=beta.new
    }
    conv="not converged"
    if(dist<eps)    conv="converged"
  }

  #8. compute psi at sort(Z1) from psi at z.obs
  z.full=sort(data.all$Z1)
  psi.full=disoph.BTF2(n.total,m,psi.new,z.full,z.obs,shape)

  #return(list(est=est, exp.beta=exp.beta, conv=conv,
  #            psi=psi.obs, z=z.obs, z.range=z.range, K=K, shape=shape, n=n, nevent=sum(DELTA), njump=m,
  #            psi.full=psi.full, z.full=z.full))

  iso.cov=data.frame(z=z.full,psi.hat=psi.full)
  colnames(iso.cov)=c(ZNAME[1],"psi.hat")

  beta.res=NA
  if(Q>0){
    beta.res=data.frame(est=beta.new, HR=exp(beta.new))
    rownames(beta.res)=ZNAME[-1]
  }

  res=list(iso.cov=iso.cov,beta=beta.res,
           conv=conv,iter=iter,Zk=Zk,shape=shape)
}
