bivisoph=function(formula, bshape="increasing", data=NULL, maxiter=10^4, eps=10^-3){

  if(bshape=="inc") bshape="increasing"
  if(bshape=="dec") bshape="decreasing"

  #1. Surv response and cov
  mf=model.frame(formula=formula, data=data) #data must be data.frame
  surv.y=model.response(mf)
  if(ncol(surv.y)==2){ #for time-indepnedent
    COV.TYPE='time.indep'
    TIME=surv.y[,1]; STATUS=surv.y[,2]
  }else if(ncol(surv.y)==3){ #for time-depnedent (or time-independent)
    COV.TYPE='time.dep'
    START=surv.y[,1]; STOP=surv.y[,2];  STATUS=surv.y[,3]
  }

  #2. COV.TYPE is time.indep or time.dep
  if(COV.TYPE=='time.indep'){
    Z.mf=mf[,-1] #mf[,1]: surv outcome; mf[,2] cov
    P.mf=ncol(mf)-1 #number of covariates

    #2.1. attributes
    Z.mf.attr=list()
    if(P.mf==1){ #single cov
      if(is.null(attributes(Z.mf))){; Z.mf.attr[[1]]=NA;
      }else{; Z.mf.attr[[1]]=attributes(Z.mf);
      }
      Z.mf=matrix(Z.mf)
    }else if(P.mf>=2){ #more than 2 cov
      for(j in 1:P.mf){
        if(is.null(attributes(Z.mf[,j]))){; Z.mf.attr[[j]]=NA;
        }else{; Z.mf.attr[[j]]=attributes(Z.mf[,j]);
        }
      }
    }

    #2.2. Seperate covariates to isotonic & additional covariates
    Z.P=Z.Q=list() #Z.P for isotonic covariates; Z.Q for additional covariates
    P=Q=0 #number of isotonic & additioal covariates
    shape=K=ZNAME.P=ZNAME.Q=NULL
    for(j in 1:P.mf){
      if(is.na(Z.mf.attr[[j]][1])){
        Q=Q+1
        Z.Q[[Q]]=Z.mf[,j]
        ZNAME.Q[Q]=colnames(Z.mf)[j]
      }else{
        P=P+1
        shape[P]=Z.mf.attr[[j]]$shape
        K[P]=Z.mf.attr[[j]]$K
        ZNAME.P[P]=Z.mf.attr[[j]]$name
        Z.P[[P]]=Z.mf[,j]
      }
    }
    Z.P=do.call("cbind",Z.P)
    Z.Q=do.call("cbind",Z.Q)
    Z=cbind(Z.P,Z.Q)
    ZNAME=c(ZNAME.P,ZNAME.Q)
    if(P>0 & Q>0){; colnames(Z)=c(paste0("Z",1:P),paste0("W",1:Q))
    }else if(Q==0){; colnames(Z)=paste0("Z",1:P)
    }
    #colnames(Z)=c(paste0("iso",".",ZNAME.P),ZNAME.Q)

    if(P==0){
      stop("there are no isotonic covariates in the formula")
    }else if(Q>0){
      stop("Additional covariates are not supported for the current version of the isoSurv package")
    }else if(P==1){ #single iso cov
      res=bivisoph.ti(TIME=TIME, STATUS=STATUS, Z=Z, ZNAME=ZNAME, P=P, Q=Q, shape1=bshape, shape2=shape, K=K, maxiter=maxiter, eps=eps)
      res$call=match.call()
      res$formula=formula
      class(res)="bivisoph"
    }else if(P>=2){ #double iso covs for additive iso
      stop("More than two isotonic covariates are not supported for the current version of the isoSurv package")
    }
  }else if(COV.TYPE=='time.dep'){
    stop("time-depdent cov is not supported for the current version of the isoph package")
  }

  return(res)
}
