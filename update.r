rm(list=ls())#remove (almost) everything in the working environment
library(igraph)
library(MASS)
library(grpreg)
library(ncvreg)

p<-20   #p is the no. of covariates
rho<-0.2 #rho controls the sparsity of the graph
beta<-2 #beta controls the signal strength
q<-10  #q is the no. of variables 
n<-200 #n is the sample size

#1 simulation
#################################################################################################################
#1.1 randomly generated scale-free networks

test.graph<-barabasi.game(q,power=1,directed=FALSE)
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
theta<-rep(0,q*q*p)
dim(theta)<-c(q,q,p)
for(j in 1:q)
{
	for(k in 1:q)
	{
		for(l in 1:p)
		{
			if(G[j,k]==0) theta[j,k,l]<-0
			else
			{
				a<-rbinom(1,1,rho)
				if(a==0) theta[j,k,l]<-0
				else
				theta[j,k,l]<-2*beta*rbinom(1,1,0.5)-beta
			}
		}
	}
}

#1.3 symmertrizing theta[j,k,l]
for(j in 1:q)
{
	for(k in 1:q)
	{
		for(l in 1:p)
		{
			if(abs(theta[j,k,l])>abs(theta[k,j,l])) theta[k,j,l]<-theta[j,k,l]
			else
			if(abs(theta[j,k,l])>abs(theta[k,j,l])) theta[j,k,l]<-theta[k,j,l]
			#theta[j,k,l]<-max(theta[j,k,l],theta[k,j,l])
			#theta[k,j,l]<-max(theta[j,k,l],theta[k,j,l])			
		}
	}
}

# 1.4 generated x
x<-mvrnorm(n,rep(0,p), diag(rep(1,p)))

# 1.5 given x and theta[j,k,l],we use Gibbs sampling to generate the y
y <- matrix(0,nrow = n, ncol = q)
eta <- matrix(0,n,q)

for(i in 1:n)
{
	uniformRV<-runif(n)
	for(j in 1:q)
	{
		b<-0
		bb<-0
		for(k in 1:q)
		{
			if(k!=j)
			{
				bb<-bb+t(theta[j,k,])%*%x[i,]*y[i,k]
				k<-k+1
				bb
			}	
		}
	}
	part_1<-t(theta[j,j, ])%*%x[i, ]
	part_2<-bb
	eta[i,j]<-part_1+part_2
	s =1- 1/(1+exp(eta[i,j]))
	if(uniformRV[j]<= s) y[i,j]<-1
	else y[i,j]<-0 	
}