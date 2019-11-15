source("function.R")
library(plus)
library(Ohit)
library(Matrix)
library(xtable)

################################
# Simulation
set.seed(1111)
n=50
T=60

W=W.n(n)
x0=cbind(rep(1,n),matrix(rnorm(2*n),n,2))
xv=matrix(c(1,0.5,0.5,1),2,2)
x=matrix(rnorm(2*n*T),n*T,2)%*%chol(xv)

xn<-NULL
for(t in 2:T)
  xn<-cbind(xn,x[((t-1)*n+1):(t*n),],x[((t-2)*n+1):((t-1)*n),])
xn=xn[,37:(4*(T-1))]

xn1<-NULL
for(t in 2:T)
  xn1<-cbind(xn1,x[((t-1)*n+1):(t*n),])
xn1=xn1[,19:(2*(T-1))]

# Parameter
rho=0.2
del0=0.3
del1=0.5
gam0=c(0.5,0,0)
gam1=c(-1.5,2.5)
beta1=rep(0,2)
beta2=rep(0,2)
beta3=c(-1,2)
beta4=rep(0,2)

real=c(rho,del0,del1,gam0,gam1,beta1,beta2,beta3,beta4)
#Y_t=\rho W Y_{t}+\delta_0 Y_{t-1}+\delta_1 W Y_{t-1}
#     +\gamma_0 x0 +\gamma_1 W x0
#     +\beta_1 X_{1,t}+...+\beta_d X_{d,t}
#     +\beta_11 W X_{1,t}+...+\beta_dd W X_{d,t}+\eps_t
#xn=[X_{1,t},...,X_{d,t},t=1,...,T]
#x0=[1,X_1,...,X_{d0-1}]

#M1##################################################
N=1000
els.oracle=matrix(0,8,N)
els.lasso=els.mcp=els.scad=els.oga=matrix(0,16,N)
for (i in 1:N) {
  
  y=matrix(0,n,T)
  y[,1]=rnorm(n)
  
  # veps=rnorm(n)
  # veps=sqrt(1/3)*rt(n,3)
  for(t in 2:T)
  {
    y[,t]=solve(diag(n)-rho*W)%*%(del0*y[,t-1]+del1*W%*%y[,t-1]+x0%*%gam0+W%*%x0[,-1]%*%gam1
                    +x[((t-2)*n+1):((t-1)*n),]%*%beta2+x[((t-1)*n+1):(t*n),]%*%beta1
                    +W%*%x[((t-1)*n+1):(t*n),]%*%beta3+W%*%x[((t-2)*n+1):((t-1)*n),]%*%beta4
                    +sqrt(1/3)*rt(n,3))
  }
  
  y=y[,11:T]
  
  elsvs0=ELS(xn1,y,W,x0)
  elsvs1=ELSVS(xn,y,W,x0,"lasso")
  elsvs2=ELSVS(xn,y,W,x0,"mc+")
  elsvs3=ELSVS(xn,y,W,x0,"scad")
  elsvs4=ELSVS(xn,y,W,x0,"oga")
  
  els.oracle[,i]=elsvs0$coef
  els.lasso[,i]=elsvs1$coef
  els.mcp[,i]=elsvs2$coef
  els.scad[,i]=elsvs3$coef
  els.oga[,i]=elsvs4$coef
}

result1=RESULT(els.oracle,els.lasso,els.mcp,els.scad,els.oga,real)
result1

