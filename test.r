range(y)
i<-3
j<-6
yt<-rep(0.01,q_)
u<-0
for(k in 1:q_){
    if(k!=j) 
    u<-u+theta[j,k,]%*%x[i,]*yt[k]
	k<-k+1   ##q-1个相加
}
u
y[3,6]