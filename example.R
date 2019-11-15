source("function.R")
library(plus)
library(Ohit)
library(Matrix)
library(xtable)
#################################
#'@examples
set.seed(1111)
n=100; T=100; d0=10; d=10; k=5
k<2*d0+4*d

W=W.n(n)
x0<-cbind(rep(1,n),matrix(rnorm(d0*n),n,d0))
xv=matrix(0.5,d,d); diag(xv)=1
x=matrix(rnorm(d*n*T),n*T,d)%*%chol(xv)

xn<-NULL
for(t in 2:T)
  xn<-cbind(xn,x[((t-1)*n+1):(t*n),],x[((t-2)*n+1):((t-1)*n),])

rho=0.3
del0=0.2
del1=0.5
alp=0.5

nonzero = sample(2*d0+4*d, k)
beta=rep(0,2*d0+4*d)
beta[nonzero]=round(runif(length(nonzero),-3,3),1)
#beta = 1.5 * (1:(2*d0+4*d) %in% nonzero)

gam0=beta[1:d0]
gam1=beta[(d0+1):(2*d0)]
beta1=beta[(1+2*d0):(2*d0+d)]
beta2=beta[(1+2*d0+d):(2*d0+2*d)]
beta3=beta[(1+2*d0+2*d):(2*d0+3*d)]
beta4=beta[(1+2*d0+3*d):(2*d0+4*d)]
real=c(rho,del0,del1,alp,beta)

y=matrix(0,n,T)
y[,1]=round(rnorm(n))

# veps=rnorm(n)
# veps=sqrt(1/3)*rt(n,3)
for(t in 2:T)
{
  y[,t]=round(solve(diag(n)-rho*W)
              %*%(del0*y[,t-1]+del1*W%*%y[,t-1]+x0%*%c(alp,gam0)+W%*%x0[,-1]%*%gam1
                  +x[((t-2)*n+1):((t-1)*n),]%*%beta2+x[((t-1)*n+1):(t*n),]%*%beta1
                  +W%*%x[((t-1)*n+1):(t*n),]%*%beta3+W%*%x[((t-2)*n+1):((t-1)*n),]%*%beta4
                  +rnorm(n)))
}

y=y[,11:T]
xn=xn[,(9*2*d+1):((2*d)*(T-1))]

# ELS 
elsvs=ELSVS(xn,y,W,x0,"lasso")
elsvs1=ELSVS(xn,y,W,x0,"mc+")
elsvs2=ELSVS(xn,y,W,x0,"scad")
elsvs3=ELSVS(xn,y,W,x0,"oga")

elsvs$coef
elsvs1$coef
elsvs2$coef
elsvs3$coef
real

elsvs$selected
elsvs1$selected
elsvs2$selected
elsvs3$selected
which(real!=0)



# IV
ive=IV(xn,y,W,x0,"lasso")
ive1=IV(xn,y,W,x0,"mc+")
ive2=IV(xn,y,W,x0,"scad")
ive3=IV(xn,y,W,x0,"oga")


ive$coef
ive1$coef
ive2$coef
ive3$coef
real

ive$selected
ive1$selected
ive2$selected
ive3$selected
which(real!=0)

