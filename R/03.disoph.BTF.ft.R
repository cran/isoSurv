disoph.BTF1=function(n,n1,lambda.upd,X1,x1.obs,shape1){
  lambda.full=rep(NA,n)
  for(i in 1:n1)
    lambda.full[which(X1 %in% x1.obs[i])]=lambda.upd[i]

  if(shape1=="increasing"){
    if(is.na(lambda.full[1]))
      lambda.full[1]=0
    which.na1=which(is.na(lambda.full))
    if(length(which.na1)>0){
      for(j in which.na1) #right continuous
        lambda.full[j]=lambda.full[j-1]
    }
  }else if(shape1=="decreasing"){
    if(is.na(lambda.full[n]))
      lambda.full[n]=0
    which.na1=which(is.na(lambda.full))
    if(length(which.na1)>0){
      for(j in sort(which.na1,decreasing=TRUE)) #left continuous
        lambda.full[j]=lambda.full[j+1]
    }
  }
  return(lambda.full)
}

disoph.BTF2=function(n,n2,psi.upd,Z2,z2.obs,shape2){
  psi.full=rep(NA,n)
  for(j in 1:n2)
    psi.full[which(Z2 %in% z2.obs[j])]=psi.upd[j]

  if(shape2=="increasing"){
    if(is.na(psi.full[1]))
      psi.full[1]=-Inf
    which.na2=which(is.na(psi.full))
    if(length(which.na2)>0){
      for(j in which.na2) #right continuous
        psi.full[j]=psi.full[j-1]
    }
  }else if(shape2=="decreasing"){
    if(is.na(psi.full[n]))
      psi.full[n]=-Inf
    which.na2=which(is.na(psi.full))
    if(length(which.na2)>0){
      for(j in sort(which.na2,decreasing=TRUE)) #left continuous
        psi.full[j]=psi.full[j+1]
    }
  }
  return(psi.full)
}
