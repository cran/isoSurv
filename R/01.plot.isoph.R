plot.isoph=function(x,...){
  y=exp(x$iso.cov[,2]) #HR
  z=x$iso.cov[,1]
  zlab=colnames(x$iso.cov)[1]
  ylab="Isotonic Covariate HR"

  if(x$shape=="increasing"){
    type='s'
  }else if(x$shape=="decreasing"){
    type='S'
  }

  plot(y=y,x=z, xlab=zlab, ylab=ylab, type=type)
}
