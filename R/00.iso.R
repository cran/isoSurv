iso=function(z,shape){
  shape=tolower(shape)
  if(shape=="inc") shape="increasing"
  if(shape=="dec") shape="decreasing"
  zname=deparse(substitute(z))
  K=stats::median(z)

  attributes(z) = c(list(name=zname,shape=shape,K=K))
  class(z)="iso covariate"
  return(z=z)
}
