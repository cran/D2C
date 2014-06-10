pcor1<-function(x,y,z){
  ## partial correlation cor(x,y|z)
  ## 14/10/11
  
  if (is.numeric(z)){
    rho.xy<-cor(x,y,"pairwise.complete.obs")
    rho.xz<-cor(x,z,"pairwise.complete.obs")
    rho.yz<-cor(z,y,"pairwise.complete.obs")
    if (is.na(rho.xz+rho.yz+rho.xy))
      return(0)
    if (rho.xz==1 | rho.yz==1)
      return(0)
    rho<-(rho.xy-rho.xz*rho.yz)/(sqrt(1-min(rho.xz^2,0.99))*sqrt(1-min(rho.yz^2,0.99)))
    return(rho)
  } else {
    stop("z should be numeric")
    
  }
  
  
  
}


corDC<-function(X,Y){
  ## correlation continuous matrix and discrete vector
  ## NB: the notion of sign has no meaning in this case. Mean of absolute values is taken
  ## 14/11/2011
  
  if (!is.factor(Y))
    stop("This is not the right function. Y is not a factor !!")
  
  N<-NROW(X)
  L<-levels(Y)
  
  if( length(L)==2)
    lL<-1
  else
    lL<-length(L)
  
  cxy<-NULL
  for (i in 1:lL){
    yy<-numeric(N)
    ind1<-which(Y==L[i])
    ind2<-setdiff(1:N,ind1)
    yy[ind1]<-1
    cxy<-cbind(cxy,abs(cor(X,yy)))
  }
  
  apply(cxy,1,mean)
}


Icond<-function(x,y=NULL,z,lambda=0){
  ## conditional  information cor(x,y|z)
  ## 14/10/11
  ## added x=matrix case on 25/10/11 
  
  
  
  ## numeric z
  if (is.numeric(z)){
    if (is.vector(x))
      return(cor2I2(pcor1(x,y,z)))
    X<-x
    n<-NCOL(X)
    Ic<-array(0,c(n,n))
    for (i in 1:(n-1))
      for (j in (i+1):n){
        Ic[i,j]<-Icond(X[,i],X[,j],z)
        Ic[j,i]<-Ic[i,j]
      }
    return(Ic)
    
  }
  ## factor z and vectors x and y 
  if (! is.null(y)){
    L<-levels(z)
    lL<-length(L)
    w<-numeric(lL)
    for (i in 1:lL)
      w[i]<-length(which(z==L[i]))
    w<-w/sum(w)
    
    Ic<-NULL
    for (i in 1:lL){
      ind1<-which(z==L[i])
      Ic<-c(Ic,cor2I2(cor(x[ind1],y[ind1])))
    }
    
    return(as.numeric(w*Ic))
  }
  
  ## factor z and matrix x
  X<-x
  n<-NCOL(X)
  L<-levels(z)
  lL<-length(L)
  w<-numeric(lL)
  for (i in 1:lL)
    w[i]<-length(which(z==L[i]))
  w<-w/sum(w)
  
  Ic<-array(0,c(n,n))
  W<-0
  for (i in 1:lL){
    ind1<-which(z==L[i])
    
    
    if (length(ind1)>8){
      
      
      Ic<-Ic+w[i]*cor2I2(cor.shrink(X[ind1,],lambda=lambda,verbose=F))
      W<-W+w[i]
    }
  }
  
  
  return(Ic/W)
  
}







  
  
ppears<-function(r.hat,N,S=0){
  n<-length(r.hat)
  p<-numeric(n)
  
  for (i in 1:n){
    z<-abs(0.5*(log(1+r.hat[i])-log(1-r.hat[i])))*sqrt(N[i]-S-3)
    
    p[i]<-pnorm(z,lower.tail=F)
    
    
  }
  p
}

corXY<-function(X,Y){
  ## correlation continuous matrix and continuous/discrete vectormatrix
  ## 14/11/2011
  
  n<-NCOL(X)
  N<-NROW(X)
  m<-NCOL(Y)
  
  cXY<-array(NA,c(n,m))
  
  for (i in 1:m){
    if (m==1)
      YY<-Y
    else
      YY<-Y[,i]
    if (is.numeric(YY)){
      cXY[,i]<-cor(X,YY,use="pairwise.complete.obs")
    } else {
      cXY[,i]<-corDC(X,YY)
    }
  }
  cXY
}


cor2I2<-function(rho){
  rho<-pmin(rho,1-1e-5)
  rho<-pmax(rho,-1+1e-5)
  -1/2*log(1-rho^2)
  
  
}


lazy.pred<- function(X,Y,X.ts,class=FALSE,return.more=FALSE,
                     conPar=3,linPar=5,cmbPar=10){

  n<-NCOL(X)
  N<-NROW(X)
  
  if (class){ ## classification
    l.Y<-levels(Y)
    L<-length(l.Y)
    u<-unique(Y)
    
    if (length(u)==1){
      P<-array(0,c(NROW(X.ts),L))
      colnames(P)<-l.Y
      P[,u]<-1
      out.hat<-factor(rep(as.character(u),length(X.ts)),levels=l.Y)
      return(list(pred=out.hat,prob=P))
    }
    
    if (L==2) {
      
      
      stop("not supported")
    } else {
      algo="lazy"
      
      stop("not supported")
      
    }
  } else { ## regression
    d<-data.frame(cbind(Y,X))
    names(d)[1]<-"Y"
    names(d)[2:(n+1)]<-paste("x",1:n,sep="")
    
    
    
    mod<-lazy(Y~.,d,control=lazy.control(distance="euclidean",
                                         conIdPar=conPar,
                                         linIdPar=linPar,
                                         cmbPar=cmbPar))
    if (is.vector(X.ts) & n>1)
      X.ts<-array(X.ts,c(1,n))
    d.ts<-data.frame(X.ts)
    
    names(d.ts)<-names(d)[2:(n+1)]
    
    if (!return.more){
      ll<- predict(mod,d.ts)
      return(ll$h)
      
    } else {
      ll<- predict(mod,d.ts,S.out=TRUE,k.out=FALSE)
      return(ll)
    }
  }
  
}




mimr<-function(X,Y,nmax=5,
               init=FALSE,lambda=0.5,
               spouse.removal=TRUE,
               caus=1){
  ## mimr filter
  # if caus =1 it searches for causes otherwise if caus=-1 it searches for effects
  
  
  NMAX<-nmax
  m<-NCOL(Y) # number of outputs
  n<-NCOL(X)
  orign<-n
  N<-NROW(X)
  H<-apply(X,2,var)
  HY<-var(Y)
  CY<-corXY(X,Y)
  Iy<-cor2I2(CY)  
  subset<-1:n
  pv.rho<-ppears(c(CY),N+numeric(n))
  if (spouse.removal){
    pv<-ppears(c(CY),N+numeric(n))
    s<-sort(pv,decreasing=TRUE,index.return=TRUE)  
    hw<-min(n-nmax,length(which(s$x>0.05)))
    spouse<-s$ix[1:hw]
    subset<-setdiff(1:n,s$ix[1:hw])
    X<-X[,subset]
    Iy<-Iy[subset]
    n<-NCOL(X)
  }
  
  
  CCx<-cor(X) 
  Ix<-cor2I2(CCx)
  Ixx<-Icond(X,z=Y,lambda=0.02)
  Inter<-array(NA,c(n,n))
  
  if (init){    
    max.kj<--Inf
    for (kk in 1:(n-1)){
      for (jj in (kk+1):n){
        Inter[kk,jj]<- (1-lambda)*(Iy[kk]+Iy[jj])+caus*lambda*(Ixx[kk,jj]-Ix[kk,jj])
        Inter[jj,kk]<-Inter[kk,jj]
        if (Inter[kk,jj]>max.kj){
          max.kj<-Inter[kk,jj]
          subs<-c(kk,jj)
        }
      }
    }
  } else {
    subs<-which.max(Iy)
  }
  
  if (nmax>length(subs)){
    last.subs<-0
    for (j in length(subs):min(n-1,NMAX-1)){
      mrmr<-numeric(n)-Inf
      if (length(subs)<(n-1)){
        if (length(subs)>1){
          mrmr[-subs]<- (1-lambda)*Iy[-subs]+caus*lambda*apply(-Ix[subs,-subs]+Ixx[subs,-subs],2,mean)
        } else {
          mrmr[-subs]<- (1-lambda)*Iy[-subs]+caus*lambda*(-Ix[subs,-subs]+Ixx[subs,-subs])
        }
      } else {
        mrmr[-subs]<-Inf
      }
      s<-which.max(mrmr)
      subs<-c(subs,s)  
    }
    
    
    
    
  }
  
  ra<-subset[subs]
  
  if (nmax>length(ra))
    ra<-c(ra,setdiff(1:orign,ra))
  
  ra
  
}



mrmr<-function(X,Y,nmax=5,first=NULL,back=FALSE){
  ## mRMR filter
  # 17/10/11
  
  
  n<-NCOL(X)
  N<-NROW(X)
  m<-NCOL(Y)
  
  if (is.factor(X[,1]) && is.factor(Y)){
    Iy<-numeric(n)
    Ix<-array(0,c(n,n))
    for (i in 1:n){
      for (j in 1:n){
        Ix[i,j]<-mutinformation(X[,i],X[,j])
      }
      Iy[i]<-mutinformation(X[,i],Y)
    }
  }else {
    
    X<-scale(X)
    Iy<-cor2I2(corXY(X,Y))
    
    CCx<-cor(X,use="pairwise.complete.obs")
    Ix<-cor2I2(CCx)
    
  }
  
  subs<-which.max(Iy)
  for (j in length(subs):min(n-1,nmax)){
    mrmr<-numeric(n)-Inf
    if (length(subs)<(n-1)){
      if (length(subs)>1){
        mrmr[-subs]<- Iy[-subs]+apply(-Ix[subs,-subs],2,mean)
      } else {
        
        mrmr[-subs]<- Iy[-subs]+(-Ix[subs,-subs])
        
      }
    } else {
      mrmr[-subs]<-Inf
    }
    
    s<-which.max(mrmr)
    subs<-c(subs,s)
    
  }
  
  
  if (back){  ## backward reordering based on linear regression
    nsubs<-NULL
    while (length(subs)>1){
      pd<-numeric(length(subs))
      
      for (ii in 1:length(subs))
        pd[ii]<-regrlin(X[,setdiff(subs,subs[ii])],Y)$MSE.emp
      
      nsubs<-c(subs[which.min(pd)],nsubs)
      subs<-setdiff(subs,subs[which.min(pd)])
      
      
    }
    subs<-c(subs,nsubs)
  }
  
  
  
  subs[1:nmax]
  
  
}


assoc <-function(x,y){
  c(abs(cor(x,y)),cor.test(x,y)$p.value)
  ##c(cor.test(x,y,method="kendall")$p.value)
}




