# Basic function=============================
# The weight matrix W_n
W.n=function(n){
  
  # Create an adjacency matrix
  # A=(a_{ij})_{n*n}, a_{ii}=0, a_{ij}=a{ji}, j \neq i
  p=10/n
  A=matrix(rbinom(n*n,1,p),n,n)
  for(i in 1:(n-1))
  {
    A[(i+1):n,i]=A[i,(i+1):n]
  }
  diag(A)=0

  i=1
  while (i<=nrow(A))
    if(sum(A[i,])==0){
      A=A[,-i];
      A=A[-i,];
    } else{
      i=i+1}
  n=nrow(A)
  
  # Rownormalization
  Wn=A
  for(i in 1:n){
    Wn[i,]=A[i,]/sum(A[i,])}
  return(Wn)
}

# The mean squared errors
MSE=function(g,c){
  N=ncol(g)
  d=nrow(g)
  se=rep(0,N)
  for (i in 1:N) {
    se[i]=crossprod(g[,i]-c)
  }
  return(MSE=sum(se)/N)
}

# Print results
RESULT=function(result.oracle,
                result.lasso,
                result.mcp,
                result.scad,
                result.oga,
                real){
  
  sub0=which(real==0)
  real0=real[-sub0]
  d=length(real)
  N=ncol(result.oracle)
  
  # Mean
  m=rep(0,5)
  m[1]=mean(result.oracle[1,])
  m[2]=mean(result.lasso[1,])
  m[3]=mean(result.mcp[1,])
  m[4]=mean(result.scad[1,])
  m[5]=mean(result.oga[1,])
  m=round(m,3)
  
  # MSE
  mse0=MSE(result.oracle,real0)
  mse1=MSE(result.lasso,real)
  mse2=MSE(result.mcp,real)
  mse3=MSE(result.scad,real)
  mse4=MSE(result.oga,real)
  mse=round(c(mse0,mse1,mse2,mse3,mse4),3)
  
  
  # Count the number of zero correctly for T
  c=c(8,rep(0,4))
  c[2]=length(which(result.lasso[sub0,]==0))/(length(sub0)*N)*8
  c[3]=length(which(result.mcp[sub0,]==0))/(length(sub0)*N)*8
  c[4]=length(which(result.scad[sub0,]==0))/(length(sub0)*N)*8
  c[5]=length(which(result.oga[sub0,]==0))/(length(sub0)*N)*8
  
  # Count the number of zero erroneously for T
  ic=rep(0,5)
  ic[2]=length(which(result.lasso[-sub0,]==0))/((16-length(sub0))*N)*8
  ic[3]=length(which(result.mcp[-sub0,]==0))/((16-length(sub0))*N)*8
  ic[4]=length(which(result.scad[-sub0,]==0))/((16-length(sub0))*N)*8
  ic[5]=length(which(result.oga[-sub0,]==0))/((16-length(sub0))*N)*8
  
  
  result=cbind(m,mse,c,ic)
  row.names(result)=c("oracle","lasso","mcp","scad","oga")
  return(result)
}


# AIC
AIC<-function(xn,y,W,x0,coef){
  
  n=nrow(y)
  T=ncol(y)
  d=ncol(xn)/T
  
  y1.vec=as.vector(y[,-1])
  y2.vec=as.vector(y[,-T])
  Wy1.vec=as.vector(W%*%y[,-1])
  Wy2.vec=as.vector(W%*%y[,-T])
  x0.vec=kronecker(rep(1,T-1),x0)
  wx0=W%*%x0[,-1]
  wx0.vec=kronecker(rep(1,T-1),wx0)
  p.coef=NULL
  xn.vec=NULL
  for (t in 2:T) {
    xn.vec=rbind(xn.vec,xn[,(d*(t-1)+1):(d*t)])
  }
  
  wxn.vec=NULL
  for (t in 2:T) {
    wxn.vec=rbind(wxn.vec,W%*%xn[,(d*(t-1)+1):(d*t)])
  }
  
  
  X.vec=cbind(Wy1.vec,y2.vec,Wy2.vec,
              x0.vec,wx0.vec,
              xn.vec,wxn.vec)
  n1=length(y1.vec)
  error=crossprod(y1.vec-X.vec%*%coef)/n1
  aic=2*length(which(coef!=0))/n1+log(error)
  return(aic=aic)
}

# BIC
BIC<-function(xn,y,W,x0,coef){
  
  n=nrow(y)
  T=ncol(y)
  d=ncol(xn)/T
  
  y1.vec=as.vector(y[,-1])
  y2.vec=as.vector(y[,-T])
  Wy1.vec=as.vector(W%*%y[,-1])
  Wy2.vec=as.vector(W%*%y[,-T])
  x0.vec=kronecker(rep(1,T-1),x0)
  wx0=W%*%x0[,-1]
  wx0.vec=kronecker(rep(1,T-1),wx0)
  p.coef=NULL
  xn.vec=NULL
  for (t in 2:T) {
    xn.vec=rbind(xn.vec,xn[,(d*(t-1)+1):(d*t)])
  }
  
  wxn.vec=NULL
  for (t in 2:T) {
    wxn.vec=rbind(wxn.vec,W%*%xn[,(d*(t-1)+1):(d*t)])
  }
  
  
  X.vec=cbind(Wy1.vec,y2.vec,Wy2.vec,
              x0.vec,wx0.vec,
              xn.vec,wxn.vec)
  n1=length(y1.vec)
  error=crossprod(y1.vec-X.vec%*%coef)/n1
  bic=log(n1)*length(which(coef!=0))/n1+log(error)
  return(bic=bic)
}

#ELS==============================================
ELSVS<-function(xn,y,W,x0,method=c("mc+","scad","lasso","oga"))
{
  #Y_t=\rho W Y_{t}+beta_0+\beta_1 X_{1,t}+...+\beta_d X_{d,t}
  #      +\delta_0 Y_{t-1}+\delta_1 W Y_{t-1}+\eps_t
  #xn=[X_{1,t},...,X_{d,t},t=1,...,T]
  #x0=[1,X_1,...,X_{d0-1}]
  
  # Validate input dimensions
  n=nrow(y)
  T=ncol(y)
  d0=ncol(x0)
  d=ncol(xn)/T
  
  # Multiply by eigenvectors
  eg=eigen(t(W))
  egv=eg$vectors
  ege=eg$values
  
  Xn<-t(egv)%*%xn
  X0<-t(egv)%*%x0[,-1]
  Y=t(egv)%*%y
  
  
  # Create variables by ELS
  Yv=NULL
  Xv=NULL
  for(i in 1:n)
  {  
    Xi=cbind(Y[i,1:(T-1)],rep(1,T-1),matrix(Xn[i,-c(1:d)],(T-1),d,byrow=T))
    hbeta=lm(Y[i,2:T]~Xi-1)$coefficients
    YXiw=chol(t(Xi)%*%Xi/(T-1))
    Yvi=YXiw%*%hbeta
    Tvi=matrix(0,(d+2),(2*d+2*d0+1))
    Tvi[1,1:2]=c(1,ege[i])
    Tvi[2,3:(d0+2)]=t(egv[,i])%*%x0
    Tvi[2,(d0+3):(2*d0+1)]=ege[i]*t(egv[,i])%*%x0[,-1]
    Tvi[3:(d+2),(2*d0+2):(2*d0+d+1)]=diag(rep(1,d))
    Tvi[3:(d+2),(2*d0+d+2):(2*d0+2*d+1)]=ege[i]*diag(rep(1,d))
    Xvi=cbind(ege[i]*Yvi,YXiw%*%Tvi)
    Yv=c(Yv,Yvi)
    Xv=rbind(Xv,Xvi)
  }
  
  
  if(method=='mc+'){
    object <- plus(Xv,Yv,method="mc+",gamma=2.4,intercept = F,normalize =F,eps = 1e-50)
    bic=log(dim(Xv)[1])*object$dim+dim(Xv)[1]*log(as.vector((1-object$r.square)*sum(Yv^2))/length(Yv))
    step.bic=which.min(bic)  
    p.coef <- coef(object, lam =object$lam[step.bic])

  }else if (method=='scad'){
    object <- plus(Xv,Yv,method="scad",gamma=2.4,intercept = F,normalize =F,eps = 1e-50)
    bic=log(dim(Xv)[1])*object$dim+dim(Xv)[1]*log(as.vector((1-object$r.square)*sum(Yv^2))/length(Yv))
    step.bic=which.min(bic)  
    p.coef <- coef(object, lam =object$lam[step.bic])
    
  }else if (method=='lasso'){
    object <- plus(Xv,Yv,method="lasso",gamma=2.4,intercept = F,normalize =F,eps = 1e-50)
    bic=log(dim(Xv)[1])*object$dim+dim(Xv)[1]*log(as.vector((1-object$r.square)*sum(Yv^2))/length(Yv))
    step.bic=which.min(bic)  
    p.coef <- coef(object, lam =object$lam[step.bic])

  }else if (method=='oga'){
    ohit=Ohit(Xv, Yv, Kn = NULL, c1 = 5, HDIC_Type = "HDHQ", c2 = 2, c3 = 2.01,intercept = FALSE)
    p.coef=rep(0,2*d0+2*d+2)
    p.coef[ohit$J_Trim]=ohit$betahat_Trim$coefficients[,1]
  }
  selected = which(p.coef!= 0)
  structure(list(Yv=Yv,Xv=Xv,coef=p.coef,selected = selected))
  
}

# Oracle
ELS<-function(xn1,y,W,x0)
{
  #Y_t=\rho W Y_{t}+ wx0 beta_0+\beta_1 X_{1,t}+...+\beta_d X_{d,t}
  #      +\delta_0 Y_{t-1}+\delta_1 W Y_{t-1}+\eps_t
  #xn=[X_{1,t},...,X_{d,t},t=1,...,T]
  #x0=[1,X_1,...,X_{d0-1}]
  n=nrow(y)
  T=ncol(y)
  d0=ncol(x0)
  d=ncol(xn1)/T
  
  eg=eigen(t(W))
  egv=eg$vectors
  ege=eg$values
  
  Xn<-t(egv)%*%xn1
  Y=t(egv)%*%y
  
  
  
  Yv=NULL
  Xv=NULL
  p.coef=NULL
  
  for(i in 1:n)
  {  
    Xi=cbind(Y[i,1:(T-1)],rep(1,T-1),matrix(Xn[i,-c(1:d)],(T-1),d,byrow=T))
    hbeta=lm(Y[i,2:T]~Xi-1)$coefficients
    YXiw=chol(t(Xi)%*%Xi/(T-1))
    Yvi=YXiw%*%hbeta
    
    Tvi=matrix(0,(d+2),(d+d0+2))
    Tvi[1,1:2]=c(1,ege[i])
    Tvi[2,3:(d0+2)]=cbind(t(egv[,i])%*%x0[,1],ege[i]*t(egv[,i])%*%x0[,-1])
    Tvi[3:(d+2),(d0+d+1):(d0+d+2)]=ege[i]*diag(rep(1,d))
    Xvi=cbind(ege[i]*Yvi,YXiw%*%Tvi)
    Yv=c(Yv,Yvi)
    Xv=rbind(Xv,Xvi)
  }
  coef=lm(Yv~Xv-1)$coefficients
  return(list(Yv=Yv,Xv=Xv,coef=coef))
  
}
#IV=================================
IV<-function(xn,y,W,x0,method=c("mc+","scad","lasso")){
  
  n=nrow(y)
  T=ncol(y)
  d0=ncol(x0)
  d=ncol(xn)/T

  y1.vec=as.vector(y[,-1])
  y2.vec=as.vector(y[,-T])
  Wy1.vec=as.vector(W%*%y[,-1])
  Wy2.vec=as.vector(W%*%y[,-T])
  x0.vec=kronecker(rep(1,T-1),x0)
  wx0=W%*%x0[,-1]
  wx0.vec=kronecker(rep(1,T-1),wx0)
  p.coef=NULL
  xn.vec=NULL
  for (t in 2:T) {
    xn.vec=rbind(xn.vec,xn[,(d*(t-1)+1):(d*t)])
  }
  
  wxn.vec=NULL
  for (t in 2:T) {
    wxn.vec=rbind(wxn.vec,W%*%xn[,(d*(t-1)+1):(d*t)])
  }
  
  
  X.vec=cbind(Wy1.vec,y2.vec,Wy2.vec,
              x0.vec,wx0.vec,
              xn.vec,wxn.vec)

  rhohat=summary(lm(y1.vec~X.vec-1))$coefficients[1]
  A=W%*%solve(diag(n)-rhohat*W)

   
  xx=cbind(y2.vec,Wy2.vec,
           x0.vec,wx0.vec,
           xn.vec,wxn.vec)
  H=NULL
  h1=A%*%xx[1:n,]
  H1=h1%*%solve(t(h1)%*%h1)%*%t(h1)
  H=bdiag(H1)
  for (t in 2:(T-1)) {
    h1=A%*%xx[(n*(t-1)+1):(n*t),]
    H1=h1%*%solve(t(h1)%*%h1)%*%t(h1)
    H=bdiag(H,H1)
  }
  Z=as.matrix(H%*%Wy1.vec)
  Zstar=cbind(Z,xx)
  
  Yv=y1.vec
  Xv=Zstar

  if(method=='mc+'){
    object <- plus(Xv,Yv,method="mc+",gamma=2.4,intercept = F,normalize =F,eps = 1e-50)
    bic=log(dim(Xv)[1])*object$dim+dim(Xv)[1]*log(as.vector((1-object$r.square)*sum(Yv^2))/length(Yv))
    step.bic=which.min(bic)  
    p.coef <- coef(object, lam =object$lam[step.bic])
    
  }else if (method=='scad'){
    object <- plus(Xv,Yv,method="scad",gamma=2.4,intercept = F,normalize =F,eps = 1e-50)
    bic=log(dim(Xv)[1])*object$dim+dim(Xv)[1]*log(as.vector((1-object$r.square)*sum(Yv^2))/length(Yv))
    step.bic=which.min(bic)  
    p.coef <- coef(object, lam =object$lam[step.bic])
    
  }else if (method=='lasso'){
    object <- plus(Xv,Yv,method="lasso",gamma=2.4,intercept = F,normalize =F,eps = 1e-50)
    bic=log(dim(Xv)[1])*object$dim+dim(Xv)[1]*log(as.vector((1-object$r.square)*sum(Yv^2))/length(Yv))
    step.bic=which.min(bic)  
    p.coef <- coef(object, lam =object$lam[step.bic])
    
  }else if (method=='oga'){
    ohit=Ohit(Xv, Yv, Kn = NULL, c1 = 5, HDIC_Type = "HDHQ", c2 = 2, c3 = 2.01,intercept = FALSE)
    p.coef=rep(0,2*d0+2*d+2)
    p.coef[ohit$J_Trim]=ohit$betahat_Trim$coefficients[,1]
  }
  selected = which(p.coef!= 0)
  structure(list(Yv=Yv,Xv=Xv,coef=p.coef,selected = selected))
}



# Oracle
IVS<-function(xn1,y,W,x0){
  n=dim(y)[1]
  T=dim(y)[2]
  d=dim(xn1)[2]/T
  
  y1.vec=as.vector(y[,-1])
  y2.vec=as.vector(y[,-T])
  Wy1.vec=as.vector(W%*%y[,-1])
  Wy2.vec=as.vector(W%*%y[,-T])
  x0.vec=kronecker(rep(1,T-1),x0)
  
  xn1.vec=NULL
  for (t in 2:T) {
    xn1.vec=rbind(xn1.vec,xn1[,(2*(t-1)+1):(2*t)])
  }
  
  
  X1.vec=cbind(Wy1.vec,y2.vec,Wy2.vec,
               x0.vec,xn1.vec)
  
  rhohat=summary(lm(y1.vec~X1.vec-1))$coefficients[1]
  A=W%*%solve(diag(n)-rhohat*W)
  
  
  xx=cbind(y2.vec,Wy2.vec,
           x0.vec,xn1.vec)
  H=NULL
  h1=A%*%xx[1:n,]
  H1=h1%*%solve(t(h1)%*%h1)%*%t(h1)
  H=bdiag(H1)
  for (t in 2:(T-1)) {
    h1=A%*%xx[(n*(t-1)+1):(n*t),]
    H1=h1%*%solve(t(h1)%*%h1)%*%t(h1)
    H=bdiag(H,H1)
  }
  Z=as.matrix(H%*%Wy1.vec)
  Zstar=cbind(Z,xx)
  
  Yv=y1.vec
  Xv=Zstar
  
  coef=lm(Yv~Xv-1)$coefficients
  return(list(Yv=Yv,Xv=Xv,coef=coef))
}