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
#iter<-100 #iter is the number of iteration of y
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
x<-cbind(1,x)


#1.4 生成y
y<-matrix(sample(c(0, 1), n*q_, replace = TRUE), nrow = n, ncol = q_, byrow = TRUE) #初始化
sampler<-function(iter){
mat<-matrix(,nrow=iter,ncol=q_)  #用于放每次迭代得到的y
	for(m in 1:iter){
	yt<-rep(1,q_)
		for(i in 1:n){
			for(j in 1:q_){
			    u<-0
				for(k in 1:q_){
					if(k!=j) 
				    u<-u+theta[j,k,]%*%x[i,]*yt[k]
				    k<-k+1   ##q-1个相加
			}
			yt[j]<-mean(rnorm(1000,u,abs(yt[j]-u)))
		    yt[j]
		}
		mat[i,]<-yt[j]
	}
	return(mat)
}				  
y<-sampler(10000)
#y<-y[seq(1000,10000,by=100)]
