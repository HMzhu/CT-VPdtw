
CTWOS <- function (MX,Y,SX,ytest,width,order,deriv,span,A0) {
  result <- savgol (MX,width,order,deriv)
  D=result$D
  my.hat=MX%*%D
  my.hat=as.matrix(my.hat)
  sy.hat=SX%*%D
  sy.hat=as.matrix(sy.hat)
  my1.hat=my.hat
  sy1.rest=sy.hat
  #######################################################################################################
  n=nrow(my1.hat)
  m=ncol(my1.hat)
  hb<-rbind(sy1.rest,my1.hat)
  a=colMeans(hb)-min(hb)
  a = a
  h<-matrix(data=NA,n,m-20)
  for (i in 1:n){
    z<-matrix(data=0,span,m-20)
    distance<-vector(mode="numeric",length=span)
    for (e in 1:span){
      b = my1.hat[i,]-min(hb)
      b = b
      x = seq(1, length(a))
      result <- VPdtw(reference=a,query=b, penalty=dilation(a,e)/1,maxshift=80)
      distance[e]<-sum(abs((result$warpedQuery[11:(m-10)])-a[11:(m-10)]))
      z[e,]<-result$warpedQuery[11:(m-10)]+rep(min(hb),times=m-20)
    }
    ind<-which.min(distance)
    query<-z[ind,]
    h[i,]<-query
  }
  h=as.matrix(h)
  #########################################################################3333
  h1<-snv(h)
  y11=matrix(Y,nrow=n)
  result <- pretreat (h1)
  mx=result$mx
  xp1=result$para1
  xp2=result$para2
  result <- pretreat (y11)
  my=result$mx
  yp1=result$para1
  yp2=result$para2
  result<-plscv (h1,y11,A0,10)
  cv<-result$RMSECV
  A=result$Optlv
  res <- pls1_nipals(mx,my,A)
  P=res$P
  W=res$W
  Wstar=W %*% solve(t(P) %*% W)
  Q=res$C
  Q=matrix(Q,nrow=A)
  ##########??????ģ????###################################################################
  xn2=nrow(sy1.rest)
  distancea<-vector(mode="numeric",length=0)
  g<-matrix(data=NA,xn2,m-20)
  for (i in 1:xn2){
    z<-matrix(data=0,span,m-20)
    distance<-vector(mode="numeric",length=span)
    for (e in 1:span){
      b = sy1.rest[i,]-min(hb)
      b = b
      x = seq(1, length(a))
      result <- VPdtw(reference=a,query=b, penalty=dilation(a,e)/1,maxshift=80)
      distance[e]<-sum(abs((result$warpedQuery[11:(m-10)])-a[11:(m-10)]))
      z[e,]<-result$warpedQuery[11:(m-10)]+rep(min(hb),times=m-20)
    }
    ind<-which.min(distance)
    query<-z[ind,]
    g[i,]<-query
  }
  g1<-snv (g)
  Xtext<-g1
  xpara1=t(xp1)
  xpara2=t(xp2)
  ypara1=t(yp1)
  ypara2=t(yp2)
  ypred <- plspredtest (Wstar,Q,g1,xpara1, xpara2, ypara1, ypara2,A)
  ytest=ytest
  result <- RMSEP (ytest, ypred)
  RMSEP=result$RMSEP
  Q2=result$Q2
  result<-list(RMSEP=RMSEP,Q2=Q2,ypred=ypred)
  return(result)
}
