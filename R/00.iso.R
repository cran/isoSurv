iso=function(z,shape="increasing",K=NA){
  if(shape=="inc") shape="increasing"
  if(shape=="dec") shape="decreasing"
  zname=deparse(substitute(z))
  if(is.na(K)) K=median(z) #median anchor

  attributes(z) = c(list(name=zname,shape=shape,K=K))
  return(z=z)
}
