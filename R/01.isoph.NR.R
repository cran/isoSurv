#newton-raphson algo with p>=2
isoph.NR=function(w,beta,Q,psi.full,n,nt,Y,dN, maxiter,eps){
  iter=0
  dist=1
  while(dist>=eps){
    iter=iter+1
    if(iter>maxiter) break

    wb=w%*%beta
    Yest=matrix(0,n,1); S1=matrix(0,Q,1)
    U=matrix(0,Q,1);  H=matrix(0,Q,Q)

    for(i in 1:n){
      for(j in which(dN[i,]>=1)){
        Yest=Y[,j]*exp(psi.full+wb) #Y*estimated exp(psi.hat+wb)
        S0=sum(Yest)
        S1=matrix(t(Yest) %*% w)
        E1=S1/S0
        U=U+(w[i,]-E1)*dN[i,j]

        S2=matrix(0,Q,Q)
        for(m in 1:n)
          S2=S2+Yest[m]*(w[m,]%*%t(w[m,]))

        H=H+(-S2/S0+E1%*%t(E1))*dN[i,j]
      }
    }

    beta.new = beta - solve(H)%*%U

    #distance
    dist=sqrt(sum((beta.new-beta)^2))
    beta=beta.new
  }
  return(beta.new)
}
