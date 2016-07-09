rm(list=ls())#remove (almost) everything in the working environment
library(igraph)
library(MASS)
library(grpreg)
library(ncvreg)


p<-20  #p is the no. of covariates
rho<-0.8#rho controls the sparsity of the graph
beta<-2#beta controls the signal strength
q<-10  #q is the no. of variables 
n<-200#n is the sample size

#y = matrix(nrow = n, ncol = q^2, 0)
x<-mvrnorm(n,rep(0,p), diag(rep(0.2,p)))

#generate the scale-free netwrok and thus getting the adjaency matrix G
test.graph<-barabasi.game(q,power=1,directed=FALSE)
degree<-degree(test.graph,v=V(test.graph))
v<-as.vector(V(test.graph))
G<-matrix(0,length(v),length(v))#adjacent matrix
for(k in 1:length(v)){
    nei<-neighbors(test.graph, v=k, mode = 1)
	G[k,nei]<-1
	}
	
g0 <- graph.adjacency( G,mode="undirected")
V(g0)$color<-"lightgreen"
plot.igraph(g0,layout=layout.fruchterman.reingold,vertex.size=12,vertex.label.cex=.9,edge.arrow.size=0.3,main="True scale-free network")
	
theta<-rep(0,q*q*p)
dim(theta)<-c(q,q,p)#theta is a tensor



set.seed(1)
for(j in 1:(q-1)){
  for(k in (j+1):q){
     for(l in 1:p){
	 if(G[j,k]==0){
	 theta[j,k,l]<-0}
	 else{
	 a<-rbinom(1,1,rho)
	 if(a==0)
	 theta[j,k,l]<-0
	 else
	 theta[j,k,l]<-2*beta*rbinom(1,1,0.5)-beta}
	 theta[k,j,l]<-theta[j,k,l]}}}
	 
	 
	 
	 
	 
y <- matrix(sample(c(0, 1), n*q, replace = TRUE), nrow = n, ncol = q, byrow = TRUE)
eta <- matrix(0,n,q)
b=0
for(i in 1:n){
uniformRV <-runif(q)
yt <-runif(q)
   for(t in 1:100){
	for(j in 1:q){
		bb<-0
		for(k in 1:q){
			if(k!=j) 
			b<-t(theta[j,k,])%*%x[i,]*yt[k]
			bb<-bb+b
			k=k+1   ##q-1个相加
			bb
		}
		part_1<-t(theta[j,j, ])%*%x[i, ]
		part_2<-bb
		eta[i,j]<-part_1+part_2
		s =1- 1/(1+exp(eta[i,j]))
		if(uniformRV[j]<=s) yt[j]<-1 
		else yt[j] <- 0 
	}
	t<-t+1}
	y[i,]<-yt
	
	
}

#list(y)







