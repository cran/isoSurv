isoph.picm=function(psi,m,z.obs,zk,k, dN2, Y2, dNsum, Delta, eps,maxiter, shape){
  #picm
  iter=0
  d.e=1
  while(d.e>=eps){
    iter=iter+1
    if(iter>maxiter) break

    #picm
    den=colSums(Y2*exp(psi))
    index.zero=which(den>0) #0/0=0

    weight=c()
    for(s in 1:m)
      weight[s]=sum( (Y2[s,]*dNsum/den)[index.zero] )

    if(sum(is.na(weight))+sum(is.infinite(weight))>=1)
      break;

    if(min(weight)==0) #For trt only: Y2[,s] can be zero due to x*beta is zero
      break;

    if(shape=='increasing'){
      exp.psi.new=pava(Delta/weight, weight)
    }else if(shape=='decreasing'){
      exp.psi.new=-pava(-Delta/weight, weight)
    }
    psi.new=log(exp.psi.new)

    #distance
    d.e=sum(abs(exp(psi.new)-exp(psi)))
    psi=psi.new
  }

  #impose the anchor
  psi.new=psi.new-psi.new[k] #psi is the same as psi.new;

  conv="not converged"
  if(d.e<eps) conv="converged"

  return(list(psi.new=psi.new, conv=conv, iter=iter))
}
