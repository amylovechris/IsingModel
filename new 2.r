library(MASS)
library(ncvreg)
sigma<-rep(1,2500)
r<-mvrnorm(n=200, rep(0,2500),diag(sigma))

cov<-matrix(0,nrow=500,ncol=500)
#define the variance-covariance matrix
for(i in 1:500){
  for(j in 1:500){
  cov[i,j]=0.6^(abs(i-j))
}}
z<-mvrnorm(n=200, rep(0,500), cov)
x<-matrix(0,nrow = 200, ncol = 2500)
for(j in 1:500){
  for(k in 1:5){
  x[,5*(j-1)+k]<-(z[,j]+r[,5*(j-1)+k])/sqrt(2)
}}
e<-rnorm(200,0,1)
beta_true<-c(0,0.5, 1, 1.5, 1, 0.5,1, 1, 1, 1, 1,-1, 0, 1, 2, 1.5,-1.5, 1, 0.5, 0.5, 0.5,rep(0,1030),rep(0,1450))#前21个
y<-x%*%beta_true[-1]+e # generate y
#beta_esi<-ncvreg(x,y,penalty = "SCAD",lambda=0.001,family="gaussian")$beta[-1]
beta_esi<-ncvreg(x, y,lambda=3, family="gaussian",penalty="lasso")$beta[-1]