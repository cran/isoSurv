isoph.RPA.ti=function(n, nt, m, z, z.obs, Y, dN,shape){
  Y2=matrix(0,m,nt)
  dN2=matrix(0,m,nt)

  if(shape=="increasing"){
    z.obs=c(z.obs,Inf) #right continuous
    for(h in 1:m){
      idx=which(z.obs[h]<=z & z<z.obs[h+1]) #right continuous
      if(length(idx)==1){
        Y2[h,]=Y[idx,]
        dN2[h,]=dN[idx,]
      }else{
        Y2[h,]=colSums(Y[idx,])
        dN2[h,]=colSums(dN[idx,])
      }
    }
  }else if(shape=="decreasing"){
    z.obs=c(-Inf,z.obs) #left continuous
    for(h in 1:m){
      idx=which(z.obs[h]<z & z<=z.obs[h+1]) #left continuous
      if(length(idx)==1){
        Y2[h,]=Y[idx,]
        dN2[h,]=dN[idx,]
      }else{
        Y2[h,]=colSums(Y[idx,])
        dN2[h,]=colSums(dN[idx,])
      }
    }
  }
  return(list(Y2=Y2,dN2=dN2));
}
