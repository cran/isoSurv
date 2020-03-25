plot.bivisoph=function(x, ...){
  x1=x$iso.bh
  x2=x$iso.cov

  y1=exp(x1[,2]) #HR
  z1=x1[,1]
  zlab1="time"
  ylab1="Isotonic Baseline HR"

  y2=exp(x2[,2]) #HR
  z2=x2[,1]
  zlab2=colnames(x2)[1]
  ylab2="Isotonic Covariate HR"

  if(x$shape.bh=="increasing"){
    type1='s'
  }else if(x$shape.bh=="decreasing"){
    type1='S'
  }

  if(x$shape.cov=="increasing"){
    type2='s'
  }else if(x$shape.cov=="decreasing"){
    type2='S'
  }

  opar=par(mfrow=c(1,2))
  on.exit(par(opar))
  plot(y=y1,x=z1, xlab=zlab1, ylab=ylab1, type=type1)
  plot(y=y2,x=z2, xlab=zlab2, ylab=ylab2, type=type2)
}
