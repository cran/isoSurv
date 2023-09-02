disoph=function(formula, bshape, data=NULL, maxiter=10^4, eps=10^-3){

  bshape=tolower(bshape)
  if(bshape=="inc") bshape="increasing"
  if(bshape=="dec") bshape="decreasing"

  #1. Surv response and cov
  mf=stats::model.frame(formula=formula, data=data) #data must be data.frame
  surv.y=stats::model.response(mf)
  if(ncol(surv.y)==2){ #for time-independent
    COV.TYPE='time.indep'
    TIME=surv.y[,1]; STATUS=surv.y[,2]
  }else if(ncol(surv.y)==3){ #for time-dependent (or time-independent)
    COV.TYPE='time.dep'
    START=surv.y[,1]; STOP=surv.y[,2];  STATUS=surv.y[,3]
  }

  #2. COV.TYPE is time.indep or time.dep
  Z.P=Z.Q=list() #Z.P for isotonic covariates; Z.Q for additional covariates
  P=Q=0 #number of isotonic & additioal covariates
  shape=K=ZNAME.P=ZNAME.Q=NULL

  if(COV.TYPE=='time.indep'){
    Z.mf=mf[,-1] #mf[,1]: surv outcome; mf[,2] cov
    if(ncol(mf)==2){ #single covariate
      P.mf=1
      if(inherits(Z.mf,"iso covariate")){ #is there a better approach?
        P=1
        Z.P[[P]]=Z.mf

        Z.mf.attr=attributes(Z.mf);
        shape[P]=Z.mf.attr$shape
        K[P]=Z.mf.attr$K
        ZNAME.P[P]=Z.mf.attr$name
      }else{
        Q=Q+1
        Z.Q[[Q]]=Z.mf

        ZNAME.Q[Q]=colnames(Z.mf)
      }
    }else{ #multiple covariate
      P.mf=ncol(Z.mf)
      for(j in 1:P.mf){
        if(inherits(Z.mf[,j],"iso covariate")){
          P=P+1
          Z.P[[P]]=Z.mf[,j]

          Z.mf.attr=attributes(Z.mf[,j]);
          shape[P]=Z.mf.attr$shape
          K[P]=Z.mf.attr$K
          ZNAME.P[P]=Z.mf.attr$name
        }else{
          Q=Q+1
          Z.Q[[Q]]=Z.mf[,j]
          ZNAME.Q[Q]=colnames(Z.mf)[j]
        }
      }
    }
    if(P>0){
      Z.P=do.call("cbind",Z.P)
      colnames(Z.P)=paste0("Z",1:P)
    }
    if(Q>0){
      Z.Q=do.call("cbind",Z.Q)
      colnames(Z.Q)=paste0("W",1:Q)
    }

    shape1=bshape #shape for baseline hazard function
    shape2=shape  #shape for covarate effect function (it can be NULL, if no z is in the formla object)

    if(P>0&Q==0){ #iso covs only
      Z=Z.P
      ZNAME=ZNAME.P
      res=disoph.ti1(TIME=TIME, STATUS=STATUS, Z=Z, ZNAME=ZNAME, P=P, Q=Q, shape1=bshape, shape2=shape, K=K, maxiter=maxiter, eps=eps)
    }else if(P==0&Q>0){ #par cov only
      #Z=Z.Q
      #ZNAME=ZNAME.Q
      #res=disoph.ti2(TIME=TIME, STATUS=STATUS, Z=Z, ZNAME=ZNAME, P=P, Q=Q, shape1=bshape, shape2=shape, K=K, maxiter=maxiter, eps=eps)
    }else if(P>0&Q>0){ #iso & par covs only
      Z=cbind(Z.P,Z.Q)
      ZNAME=c(ZNAME.P,ZNAME.Q)
      res=disoph.ti3(TIME=TIME, STATUS=STATUS, Z=Z, ZNAME=ZNAME, P=P, Q=Q, shape1=bshape, shape2=shape, K=K, maxiter=maxiter, eps=eps)
    }

    res$call=match.call()
    res$formula=formula
    class(res)="disoph"

  }else if(COV.TYPE=='time.dep'){
    stop("time-depdent cov is not supported for the current version of the isoph package")
  }

  return(res)
}
