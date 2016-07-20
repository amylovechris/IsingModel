rm(list=ls())#remove (almost) everything in the working environment
library(igraph)
library(MASS)
library(Matrix)
library(grpreg)
library(ncvreg)
library(glmnet)
library(foreach)
phi<-function(p,q,n,x,y_){
                phi<-matrix(0,n,(p+1)*(q-1))
				for(i in 1:nrow(phi)){
				phi[i,]<-kronecker(y_[i,],x[i,])}
				return(phi)}
#the estimated theta (symmetric) for a single lambda
theta_hat<-function(p,q,n,x,y,lambda,rule=c("and","or"),penalty=c("lasso","SCAD")){
          
           theta_hat<-rep(0,q*q*(p+1))
           dim(theta_hat)<-c(q,q,(p+1))
		   for(j in 1:q){
		       eta<-y[,j]#the response
			   y_<-y[,-j]
			   phi_<-phi(p,q,n,x,y_)#design matrix
			   esi<-ncvreg(phi_,eta,penalty,family="gaussian",lambda=lambda,gamma=3.7)$beta[-1]
			   #esi<-glmnet(phi_,eta,family="gaussian",lambda=lambda,intercept = "FALSE")$beta
			   dim(esi)<-c(p+1,q-1)
			   theta_hat[j,-j,]<-t(esi)			   
			   }
			   theta_hat<-theta_hat_symetrize(p,q,theta_hat,rule)
			   for(s in 1:(q-1)){
					   for(k in (s+1):q){
					      for(l in 1:(p+1)){
						  theta_hat[s,k,l]<-ifelse(abs(theta_hat[j,k,l])>=0.2,theta_hat[s,k,l],0)
						  theta_hat[k,s,l]<-theta_hat[s,k,l]	  
						  
					}}}
			   list(theta_hat=theta_hat,lambda=lambda)
           		   }
				   
#symetrize the estimated theta by "and"	rule or "or" rule			   
theta_hat_symetrize<-function(p,q,theta_hat,rule=c("and","or")){
                    if(rule=="and"){
					for(j in 1:(q-1)){
					   for(k in (j+1):q){
					      for(l in 1:(p+1)){
						  theta_hat[j,k,l]<-ifelse(abs(theta_hat[j,k,l])>=abs(theta_hat[k,j,l]),theta_hat[j,k,l],theta_hat[k,j,l])
						  theta_hat[k,j,l]<-theta_hat[j,k,l]	  
						  
					}}}}
					if(rule=="or"){
					for(j in 1:(q-1)){
					   for(k in (j+1):q){
					      for(l in 1:(p+1)){
						  theta_hat[j,k,l]<-ifelse(abs(theta_hat[j,k,l])<=abs(theta_hat[k,j,l]),theta_hat[j,k,l],theta_hat[k,j,l])
						  theta_hat[k,j,l]<-theta_hat[j,k,l]	  
						  
					}}}}
					
					return(theta_hat)}

#compute the fpr and tpr for a single estimator					
fpr_tpr<-function(theta,theta_hat){
                TP<-0
				FP<-0
				T<-sum(theta!=0)/2
				F<-sum(theta==0)/2
				for(j in 1:(q-1)){
				   for(k in (j+1):q){
				      for(l in 1:(p+1)){
					  if(theta[j,k,l]!=0&theta_hat[j,k,l]!=0)
					  TP<-TP+1
					  if(theta[j,k,l]==0&theta_hat[j,k,l]!=0)
					  FP<-FP+1
					  }}}
				TPR<-TP/T
				FPR<-FP/F
				list(TPR=TPR,FPR=FPR,T=T,F=F,TP=TP,FP=FP,P=TP+FP)
                }					  
					
fpr_tpr_edge<-function(G,theta_hat){
                TP<-0
				FP<-0
				T<-sum(G!=0)/2
				F<-sum(G==0)/2
				for(j in 1:(q-1)){
				   for(k in (j+1):q){
				      
					  if(G[j,k]!=0&crossprod(theta_hat[j,k,]!=0))
					  TP<-TP+1
					  if(G[j,k]==0&crossprod(theta_hat[j,k,]!=0))
					  FP<-FP+1
					  }}
				TPR<-TP/T
				FPR<-FP/F
				list(TPR=TPR,FPR=FPR,T=T,F=F,TP=TP,FP=FP,P=TP+FP)
                }					
				
p<-20 #p is the no. of covariates
rho<-0.2#rho controls the sparsity of the graph
beta<-2#beta controls the signal strength
q<-10  #q is the no. of variables 
n<-500#n is the sample size
lambda<-0.01

#generate the scale-free netwrok and thus getting the adjaency matrix G
test.graph<-barabasi.game(q,power=1,directed=FALSE)
#test.graph<-erdos.renyi.game(q, 0.1)
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

theta<-rep(0,q*q*(p+1))
dim(theta)<-c(q,q,(p+1))#theta is a tensor

set.seed(1)
for(j in 1:(q-1)){
  for(k in (j+1):q){
     for(l in 1:(p+1)){
	 if(G[j,k]==0){
	 theta[j,k,l]<-0}
	 else{
	 a<-rbinom(1,1,rho)
	 if(a==0)
	 theta[j,k,l]<-0
	 else
	 theta[j,k,l]<-2*beta*rbinom(1,1,0.5)-beta}
	 theta[k,j,l]<-theta[j,k,l]}}}

 timestart<-Sys.time()
 l<-1#l is the no. of replications
 fpr<-NULL
 tpr<-NULL
for(m in 1:l){
cat("iteration=",m,"\n")                              	 
#generate the covariates x
x<-mvrnorm(n,rep(0,p), diag(rep(1,p)))
x<-cbind(1,x)
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
# y<-matrix(0,n,q)
# for(i in 1:n){
   # sigma_inverse<-diag(rep(0,q))
   # for(j in 1:q){
    # for(l in setdiff(1:q,j)){
   # sigma_inverse[l,j]<-(-1)*theta[j,l,]%*%x[i,]}}
   # for(s in 1:q){sigma_inverse[s,s]<-sum(abs(sigma_inverse[s,]))+1}
   # y[i,]<-rnorm(1,rep(0,q),ginv(sigma_inverse))
   # }
     # one <- rep(1, nrow(x))
			  # meanx <- drop(one %*% x)/nrow(x)
			   # x <- scale(x, meanx, FALSE)
			   # normx <- sqrt(drop(one %*% (x^2)))
			   # x <- scale(x, FALSE, normx)
			   
y <- matrix(0,nrow = n, ncol = q)

for(i in 1:n){
yt <-rnorm(q)
yt<-rep(0,q)
yt_new <-rnorm(q)
   for(t in 1:500){
	for(j in 1:q){
		u<-0
		for(k in 1:q){
			if(k!=j) 
			
			u<-u+(theta[j,k,]%*%x[i,])*y[i,k]
			k=k+1   ##q-1个相加
			
		}
		yt[j]<-mean(rnorm(500,u,1))
	}
	t<-t+1}
	y[i,]<-yt
	
	
}

config.mean = 0;
config.sd = 1;
config.number = q;
config.inter = 2000;
config.samples = n;



##################################################################################				   

 # one <- rep(1, nrow(y))
			    # meany <- drop(one %*% y)/nrow(y)
			   # y <- scale(y, meany, FALSE)
			  # # normy <- sqrt(drop(one %*% (y^2)))
			  # y <- scale(y, FALSE, normy)

				   
from<-0
to<-0.001
seg<-0.0001
for(lambda in (seq(from,to,seg))){
theta_hat1<-theta_hat(p,q,n,x,y,lambda,rule="and",penalty="SCAD")$theta_hat
fpr<-c(fpr,fpr_tpr_edge(G,theta_hat1)$FPR)
tpr<-c(tpr,fpr_tpr_edge(G,theta_hat1)$TPR)

}
}
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime) 	

fit1<-loess(tpr~fpr,span=0.8,degree=2)
TPR<-predict(fit1,seq(range(fpr)[1], range(fpr)[2], 0.005), se=TRUE)$fit
FPR<-seq(range(fpr)[1], range(fpr)[2], 0.005)
plot(FPR,TPR,col="green",xlim=c(0,1),ylim=c(0,1),xlab="False positive rate",ylab="True positive rate",lty=1,type="l",lwd=3)

theta_hat1<-theta_hat(p,q,n,x,y,lambda,rule="and",penalty="SCAD")$theta_hat
sum(theta_hat1)
fpr_tpr(theta,theta_hat1)$FPR
fpr_tpr(theta,theta_hat1)$TPR


