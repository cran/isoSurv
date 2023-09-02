disoph.NR2=function(lambda,beta,wb,W,Delta2,Q,n,WT,eps){
  dist2=1
  while(dist2>eps){
    U=rep(0,Q)
    H=matrix(0,Q,Q)
    for(j in 1:n){
      wt3.j=sum(WT[,j]*lambda)*exp(wb[j])
      U=U+(-Delta2[j]+wt3.j)*W[j,]
      H=H+wt3.j*W[j,]%*%t(W[j,])
    }
    #for(j in 1:n){ #above is the same as here (useful for C++)
    #  wt3=0
    #  for(i in 1:n)
    #    wt3=wt3+WT[i,j]*lambda[i]
    #  wt3.j=wt3*exp(wb[j])
    #  U=U+(-Delta2[j]+wt3.j)*W[j,]
    #  H=H+wt3.j*W[j,]%*%t(W[j,])
    #}

    beta.upd2=beta-solve(H)%*%U
    wb.upd2=W%*%beta.upd2
    dist2=sqrt(sum((beta.upd2-beta)^2))

    beta=beta.upd2
    wb=wb.upd2
  }

  return(list(beta=beta.upd2,wb=wb.upd2))
}

disoph.NR3=function(lambda,psi,beta,wb,W,Delta2,Q,n,WT,eps){
  dist2=1
  while(dist2>eps){
    U=rep(0,Q)
    H=matrix(0,Q,Q)
    for(j in 1:n){
      wt3.j=sum(WT[,j]*lambda)*exp(psi[j]+wb[j])
      U=U+(-Delta2[j]+wt3.j)*W[j,]
      H=H+wt3.j*W[j,]%*%t(W[j,])
    }
    #for(j in 1:n){ #above is the same as here (useful for C++)
    #  wt3=0
    #  for(i in 1:n)
    #    wt3=wt3+WT[i,j]*lambda.upd[i]
    #  wt3.j=wt3*exp(wb[j])
    #  U=U+(-Delta2[j]+wt3.j)*W[j,]
    #  H=H+wt3.j*W[j,]%*%t(W[j,])
    #}

    beta.upd2=beta-solve(H)%*%U
    wb.upd2=W%*%beta.upd2
    dist2=sqrt(sum((beta.upd2-beta)^2))

    beta=beta.upd2
    wb=wb.upd2
  }
  return(list(beta=beta.upd2,wb=wb.upd2))
}
