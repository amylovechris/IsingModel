b=2
alpha=0
rho=c(0,4)
phi=c(6,2)

sampler1=function(n,x,y){
  mat=matrix(,ncol=2,nrow=n)
  mat[1,]=c(x,y)
  for(i in 1:n){
  x=rnorm(1,rho[y+1]/b,1/sqrt(b))
  y=rbinom(1,1,exp(rho[2]*x+phi[2])/(exp(rho[1]*x+phi[1])+exp(rho[2]*x+phi[2])))
  mat[i,]=c(x,y)
  }
  mat
}

sampler2=function(n,x,y){
  mat=matrix(,ncol=2,nrow=n)
  mat[1,]=c(x,y)
  for(i in 1:n){
    y=rbinom(1,1,exp(rho[2]*x+phi[2])/(exp(rho[1]*x+phi[1])+exp(rho[2]*x+phi[2])))
    x=rnorm(1,rho[y+1]/b,1/sqrt(b))
    mat[i,]=c(x,y)
  }
  mat
}

par(mfrow=c(2,1))

sa_1=sampler1(100000,0,0)
sa_2=sampler2(100000,0,0)
sa_1=sa_1[seq(1000,100000,100),]
sa_2=sa_2[seq(1000,100000,100),]
plot(sa_1,main="order x y",xlab="x",ylab="y")
plot(sa_2,main="order y x",xlab="x",ylab="y")