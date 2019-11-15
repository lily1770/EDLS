library(Ohit)
library(Matrix)
library(plus)
library(xtable)
source(function.R)

###########################################
# Part I and Part II crimes
par1<-read.table("part1.txt")[,-1]
par1<-as.matrix(par1)
par2<-read.table("part2.txt")[,-1]
par2<-as.matrix(par2)

# Wn
Ws<-read.table("W.txt")[-1,-1]
Ws<-as.matrix(Ws)

# x0
socio<-read.table("socio.txt")[,-1]
socio<-as.matrix(socio)

sna=apply(socio,2,median,na.rm=T)
socio1=socio
for(i in 1:dim(socio)[2])
  socio1[which(is.na(socio[,i])),i]=sna[i]

x0=as.matrix(socio1)
x0=cbind(1,x0)

# log(1+c) transformed counts of Part I and Part II
ln_par1=log(par1+1)
ln_par2=log(par2+1)

# Part I without part II#########################################
y=ln_par1
n=nrow(y)
T=ncol(y)

# xn=(cos(2¦Ðt/12), sin(2¦Ðt/12))
xn<-NULL
for(t in 1:T)
  xn<-cbind(xn,rep(cos(2*pi*t/12),n),rep(sin(2*pi*t/12),n))

# Print the results of ELS
els.lasso=ELSVS(xn,y,Ws,x0,"lasso")
els.mcp=ELSVS(xn,y,Ws,x0,"mc+")
els.scad=ELSVS(xn,y,Ws,x0,"scad")
els.oga=ELSVS(xn,y,Ws,x0,"oga")

# Coefficients of lasso by the results of lm
coef=as.vector(summary(lm(els.lasso$Yv~els.lasso$Xv[,els.lasso$selected]-1))$coefficients[,1])
els.lasso$coef[els.lasso$selected]=coef

data.frame(els.lasso$coef,els.mcp$coef,els.scad$coef,els.oga$coef)

els.lasso$selected
els.mcp$selected
els.scad$selected
els.oga$selected

# The significance of the regression coefficient
Yv=els.lasso$Yv
Xv=els.lasso$Xv
summary(lm(Yv~Xv[,els.lasso$selected]-1))
summary(lm(Yv~Xv[,els.mcp$selected]-1))
summary(lm(Yv~Xv[,els.scad$selected]-1))
summary(lm(Yv~Xv[,els.oga$selected]-1))

# AIC and BIC
AIC(xn,y,Ws,x0,els.lasso$coef)
AIC(xn,y,Ws,x0,els.mcp$coef)
AIC(xn,y,Ws,x0,els.scad$coef)
AIC(xn,y,Ws,x0,els.oga$coef)

BIC(xn,y,Ws,x0,els.lasso$coef)
BIC(xn,y,Ws,x0,els.mcp$coef)
BIC(xn,y,Ws,x0,els.scad$coef)
BIC(xn,y,Ws,x0,els.oga$coef)

# K-S test
res1=lm(Yv~Xv[,which(els.lasso$coef!=0)]-1)$residuals
res2=lm(Yv~Xv[,which(els.mcp$coef!=0)]-1)$residuals
res3=lm(Yv~Xv[,which(els.scad$coef!=0)]-1)$residuals
res4=lm(Yv~Xv[,which(els.oga$coef!=0)]-1)$residuals

ks.test(res1,"pnorm",mean(res1),sd(res1))
ks.test(res2,"pnorm",mean(res2),sd(res2))
ks.test(res3,"pnorm",mean(res3),sd(res3))
ks.test(res4,"pnorm",mean(res4),sd(res4))


# residual test (LASSO)
par(mfrow=c(1,3))
fit=lm(Yv~Xv[,els.lasso$selected]-1)
hist(fit$residuals,main = "(a) Histogram of residual",xlab="residuals")
qqnorm(fit$residuals,main="(b) Q-Q normal plot")
qqline(fit$residuals)
plot(fit$residuals,ylab = "residual",main = "(c) Scatter plot of residuals")
abline(h=0)


#Part I with Part II#########################################
y=ln_par1
xn1=ln_par2
# xn=(par2_t,cos(2¦Ðt/12), sin(2¦Ðt/12))
xn<-NULL
for(t in 1:T)
  xn<-cbind(xn,xn1[,t],rep(cos(2*pi*t/12),n),rep(sin(2*pi*t/12),n))

# Print the results of ELS
els.lasso=ELSVS(xn,y,Ws,x0,"lasso")
els.mcp=ELSVS(xn,y,Ws,x0,"mc+")
els.scad=ELSVS(xn,y,Ws,x0,"scad")
els.oga=ELSVS(xn,y,Ws,x0,"oga")

# Coefficients of lasso by the results of lm
coef=as.vector(summary(lm(els.lasso$Yv~els.lasso$Xv[,els.lasso$selected]-1))$coefficients[,1])
els.lasso$coef[els.lasso$selected]=coef

data.frame(els.lasso$coef,els.mcp$coef,els.scad$coef,els.oga$coef)

els.lasso$selected
els.mcp$selected
els.scad$selected
els.oga$selected

# The significance of the regression coefficient
Yv=els.lasso$Yv
Xv=els.lasso$Xv
summary(lm(Yv~Xv[,els.lasso$selected]-1))
summary(lm(Yv~Xv[,els.mcp$selected]-1))
summary(lm(Yv~Xv[,els.scad$selected]-1))
summary(lm(Yv~Xv[,els.oga$selected]-1))

# AIC and BIC
AIC(xn,y,Ws,x0,els.lasso$coef)
AIC(xn,y,Ws,x0,els.mcp$coef)
AIC(xn,y,Ws,x0,els.scad$coef)
AIC(xn,y,Ws,x0,els.oga$coef)

BIC(xn,y,Ws,x0,els.lasso$coef)
BIC(xn,y,Ws,x0,els.mcp$coef)
BIC(xn,y,Ws,x0,els.scad$coef)
BIC(xn,y,Ws,x0,els.oga$coef)

# K-S test
res1=lm(Yv~Xv[,which(els.lasso$coef!=0)]-1)$residuals
res2=lm(Yv~Xv[,which(els.mcp$coef!=0)]-1)$residuals
res3=lm(Yv~Xv[,which(els.scad$coef!=0)]-1)$residuals
res4=lm(Yv~Xv[,which(els.oga$coef!=0)]-1)$residuals

ks.test(res1,"pnorm",mean(res1),sd(res1))
ks.test(res2,"pnorm",mean(res2),sd(res2))
ks.test(res3,"pnorm",mean(res3),sd(res3))
ks.test(res4,"pnorm",mean(res4),sd(res4))


# residual test (LASSO)
par(mfrow=c(1,3))
fit=lm(Yv~Xv[,els.lasso$selected]-1)
hist(fit$residuals,main = "(a) Histogram of residual",xlab="residuals")
qqnorm(fit$residuals,main="(b) Q-Q normal plot")
qqline(fit$residuals)
plot(fit$residuals,ylab = "residual",main = "(c) Scatter plot of residuals")
abline(h=0)

