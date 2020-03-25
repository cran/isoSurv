isoph.BTF=function(m, n, z,z.obs, psi.new,shape){
  if(shape=="increasing"){
    psi.full=rep(-Inf,n)
    for(h in 1:m)
      psi.full[z>=z.obs[h]]=psi.new[h] #right continuous
    return(psi.full);
  }else if(shape=="decreasing"){
    psi.full=rep(-Inf,n);
    for(h in m:1)
      psi.full[z<=z.obs[h]]=psi.new[h] #left continuous
    return(psi.full);
  }
}
