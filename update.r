rm(list=ls())#remove (almost) everything in the working environment
library(igraph)
library(MASS)
library(grpreg)
library(ncvreg)

p<-20   #p is the no. of covariates
rho<-0.2 #rho controls the sparsity of the graph
beta_<-2 #beta_ controls the signal strength
q_<-10  #q_ is the no. of variables 
n<-200 #n is the sample size
iter<-100 #iter is the number of iteration of y
lambda<-0.001
##################################################################################################
#1 simulation
##################################################################################################
#1.1 randomly generated scale-free networks

test.graph<-barabasi.game(q_,power=1,directed=FALSE)
degree<-degree(test.graph,v=V(test.graph))
v<-as.vector(V(test.graph))
G<-matrix(0,length(v),length(v))#adjacent matrix
for(k in 1:length(v))
{
    nei<-neighbors(test.graph, v=k, mode = 1)
	G[k,nei]<-1
}
g0 <- graph.adjacency(G,mode="undirected")
V(g0)$color<-"lightgreen"
plot.igraph(g0,layout=layout.fruchterman.reingold,vertex.size=12,vertex.label.cex=.9,edge.arrow.size=0.3,main="True scale-free network")

#1.2 update nonzero theta[j,k] in G
theta<-rep(0,q_*q_*p)
dim(theta)<-c(q_,q_,p)
set.seed(1)
for(j in 1:q_){
	for(k in 1:q_){
		for(l in 1:p){
			if(G[j,k]==0){
			theta[j,k,l]<-0}
			else{
			a<-rbinom(1,1,rho)
				if(a==0){
				theta[j,k,l]<-0}
				else{
				theta[j,k,l]<-2*beta_*rbinom(1,1,0.5)-beta_
				}
			}
			
		}
	}
}

#1.3 symmetrizing(separate-max)
for(j in 1:q_){
	for(k in 1:q_){
		for(l in 1:p){
			theta[j,k,l]<-ifelse(abs(theta[j,k,l])>abs(theta[k,j,l]),theta[j,k,l],theta[k,j,l])
			theta[k,j,l]<-theta[j,k,l]	
		}
	}
}

#1.3 生成x
x<-mvrnorm(n,rep(0,p), diag(rep(1,p)))

#1.4 生成y(Ising)
y <- matrix(sample(c(0, 1), n*q_, replace = TRUE), nrow = n, ncol = q_, byrow = TRUE)
eta <- matrix(0,n,q_)
b<-0
for(i in 1:n){
uniformRV <-runif(q_)
	for(it in 1:iter){
		for(j in 1:q_){
			bb<-0
			for(k in 1:q_){
				if(k!=j) 
				b<-t(theta[j,k,])%*%x[i,]*y[i,k]
				bb<-bb+b
				k=k+1   ##q_-1个相加
			}
			part_1<-t(theta[j,j, ])%*%x[i, ]
			part_2<-bb
			eta[i,j]<-part_1+part_2
			s =1- 1/(1+exp(eta[i,j]))
			if(uniformRV[j]<=s) y[i,j]<-1 
			else y[i,j] <-0 
		}
	it<-it+1
	}
}

#1.4 生成y(Gaussian)
y <- matrix(sample(c(0, 1), n*q_, replace = TRUE), nrow = n, ncol = q_, byrow = TRUE)
#y <- matrix(0,nrow = n, ncol = q_, byrow = TRUE)
for(i in 1:n){
	yt<-rnorm(q_,0,1)
	for(iter in 1:200){
		for(j in 1:q_){
			u<-0
			for(k in 1:q_){
			if(k!=j) 
			u<-u+t(theta[j,k,])%*%x[i,]*y[i,k]
			#u<-u+t(theta[j,k,])%*%x[i,]*yt[k]
			k<-k+1   ##q_个k相加	
		}
		yt[j]<-mean(mvrnorm(200,u,1))
	}
	iter<-iter+1}
	y[i,]<-yt
}
###################################################################################################
#2 the estimated theta
###################################################################################################
#2.1 kronecker product
phi<-function(p,q_,n,x,y_)
{
    phi<-matrix(0,n,p*(q_-1))
	for(i in 1:nrow(phi))
	{
		phi[i,]<-kronecker(y_[i,],x[i,])  #dim(kronecker(y[i,],x[i,]))=200,length(phi[i,])=180
	}
	    return(t(phi))
}

#2.2 the estimated theta (symmetric) for a single lambda
theta_hat<-function(p,q_,n,x,y,lambda,penalty="SCAD")
{
    theta_hat<-rep(0,q_*q_*p)
    dim(theta_hat)<-c(q_,q_,p)
	for(j in 1:q_)
	{
		eta<-y[,j] #the response  200*1
	    y_<-y[,-j] #200*9
		phi_<-phi(p,q_,n,x,y_) #design matrix x
		esi<-ncvreg(phi_,eta,penalty = "SCAD",family="gaussian")$beta[-1]
		dim(esi)<-c(p,q_-1)
		theta_hat[j,-j,]<-t(esi)	
	}
}

theta_hat_symetrize<-function(p,q_,theta_hat,rule=c("and","or"))
{
    if(rule=="and")
	{
	    for(j in 1:(q_-1))
		{
			for(k in (j+1):q_)
			{
				for(l in 1:p)
				{
					theta_hat[j,k,l]<-ifelse(abs(theta_hat[j,k,l])>=abs(theta_hat[k,j,l]),theta_hat[j,k,l],theta_hat[k,j,l])
				    theta_hat[k,j,l]<-theta_hat[j,k,l]	  
				}
			}
		}
	}
	
	if(rule=="or")
	{
		for(j in 1:(q_-1))
		{
			for(k in (j+1):q_)
			{
				for(l in 1:p)
				{
					theta_hat[j,k,l]<-ifelse(abs(theta_hat[j,k,l])<=abs(theta_hat[k,j,l]),theta_hat[j,k,l],theta_hat[k,j,l])
					theta_hat[k,j,l]<-theta_hat[j,k,l]	  
				}
			}
		}
	}		
	return(theta_hat)
}			   
			   
theta_hat<-theta_hat_symetrize(p,q_,theta_hat,"and")
list(theta_hat=theta_hat,lambda=lambda)