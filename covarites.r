rm(list=ls())#remove (almost) everything in the working environment
library(igraph)
library(MASS)
library(grpreg)
library(ncvreg)
phi<-function(p,q,n,x,y_){
                phi<-matrix(0,n,p*(q-1))
				for(i in 1:nrow(phi)){
				phi[i,]<-kronecker(y_[i,],x[i,])}
				return(phi)}
#the estimated theta (symmetric) for a single lambda
theta_hat<-function(p,q,n,x,y,lambda,rule=c("and","or")){
          
           theta_hat<-rep(0,q*q*p)
           dim(theta_hat)<-c(q,q,p)
		   for(j in 1:q){
		       eta<-y[,j]#the response
			   y_<-y[,-j]
			   phi_<-phi(p,q,n,x,y_)#design matrix
			   esi<-ncvreg(phi_,eta,penalty="SCAD",family="gaussian",lambda=lambda,gamma=3.7)$beta[-1]
			   dim(esi)<-c(p,q-1)
			   theta_hat[j,-j,]<-t(esi)			   
			   }
			   theta_hat<-theta_hat_symetrize(p,q,theta_hat,rule)
			   list(theta_hat=theta_hat,lambda=lambda)
           		   }
#symetrize the estimated theta by "and"	rule or "or" rule			   
theta_hat_symetrize<-function(p,q,theta_hat,rule=c("and","or")){
                    if(rule=="and"){
					for(j in 1:(q-1)){
					   for(k in (j+1):q){
					      for(l in 1:p){
						  theta_hat[j,k,l]<-ifelse(theta_hat[j,k,l]>=theta_hat[k,j,l],theta_hat[j,k,l],theta_hat[k,j,l])
						  theta_hat[k,j,l]<-theta_hat[j,k,l]	  
						  
					}}}}
					if(rule=="or"){
					for(j in 1:(q-1)){
					   for(k in (j+1):q){
					      for(l in 1:p){
						  theta_hat[j,k,l]<-ifelse(theta_hat[j,k,l]<=theta_hat[k,j,l],theta_hat[j,k,l],theta_hat[k,j,l])
						  theta_hat[k,j,l]<-theta_hat[j,k,l]	  
						  
					}}}}
					
					return(theta_hat)}

#compute the fpr and tpr for a single estimator					
fpr_tpr<-function(theta,theta_hat){
                TP<-0
				FP<-0
				T<-sum(theta!=0)
				F<-sum(theta==0)
				for(j in 1:(q-1)){
				   for(k in (j+1):q){
				      for(l in 1:p){
					  if(theta[j,k,l]!=0&theta_hat[j,k,l]!=0)
					  TP<-TP+1
					  if(theta[j,k,l]==0&theta_hat[j,k,l]!=0)
					  FP<-FP+1
					  }}}
				TPR<-TP/T
				FPR<-FP/F
				list(TPR=TPR,FPR=FPR,T=T,F=F,TP=TP,FP=FP,P=TP+FP)
                }					  
					

					
					


					
					

p<-10  #p is the no. of covariates
rho<-0.8#rho controls the sparsity of the graph
beta<-2#beta controls the signal strength
q<-20  #q is the no. of variables 
n<-200#n is the sample size

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
	
#generate the true parameters theta 

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



#generate the covariates x
x<-mvrnorm(n,rep(0,p), diag(rep(0.2,p)))

#generate y
# a<-0.245#a is the signal strength 
# sigma_inverse<-diag(rep(1,q))
# for(j in 1:(q-1)){
  # for(k in (j+1):q){
  # if(all(theta[j,k,]==0)){
  # sigma_inverse[j,k]<-0}
  # else{
  # sigma_inverse[j,k]<-a}
  # sigma_inverse[k,j]<-sigma_inverse[j,k]
  
  # }}

# y<-mvrnorm(n,rep(0,q), solve(sigma_inverse))
#有问题的部分
y<-matrix(0,n,q)
for(i in 1:n){
   sigma_inverse<-diag(rep(0,q))
   for(j in 1:q){
    for(l in setdiff(1:q,j)){
   sigma_inverse[l,j]<-(-1)*theta[j,l,]%*%x[i,]}}
   for(s in 1:q){sigma_inverse[s,s]<-sum(abs(sigma_inverse[s,]))+1}
   y[i,]<-rnorm(1,rep(0,q),solve(sigma_inverse))  ##solve可以改为广义逆
   }
######################				   


fpr<-NULL
tpr<-NULL				   
from<-0.000001
to<-0.1
seg<-0.01
for(lambda in (seq(from,to,seg))){
theta_hat1<-theta_hat(p,q,n,x,y,lambda,rule="and")$theta_hat
fpr<-c(fpr,fpr_tpr(theta,theta_hat1)$FPR)
tpr<-c(tpr,fpr_tpr(theta,theta_hat1)$TPR)

}

theta_hat1<-theta_hat(p,q,n,x,y,lambda,rule="and")$theta_hat
sum(theta_hat1)
fpr_tpr(theta,theta_hat1)$FPR
fpr_tpr(theta,theta_hat1)$TPR


