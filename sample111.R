neighbours<-function(spin.conf,x,y,size){  #counter[1] is black,counter[2] is white,(x,y)is the site.
  counter<-vector("numeric",2) # Count positive and negative neighbours
  if(x+1<=size) { # Check right neighbour
    if(spin.conf[x+1,y]==1) counter[1]<-1 else counter[2]<-1
  }
  if(x>=2) { # Check left neighbor
    if(spin.conf[x-1,y]==1) counter[1]<-counter[1]+1 else counter[2]<-counter[2]+1
  }
  if(y+1<=size) { # Check lower neighbor
    if(spin.conf[x,y+1]==1) counter[1]<-counter[1]+1 else counter[2]<-counter[2]+1
  }
  if(y>=2) { # Check upper neighbor
    if(spin.conf[x,y-1]==1) counter[1]<-counter[1]+1 else counter[2]<-counter[2]+1
  }
  return(list(counter[1],counter[2]))
}

size=10;beta = 0.8;num=200;
# 初始化
spin.conf <- matrix(sample(c(0,1),size^2,replace = TRUE), nrow = size, ncol = size)
counter <- vector("numeric",2)
index <- sample(c(1:size^2), size^2, replace = FALSE)
samples = matrix(nrow=num,ncol=size^2,0)

for(j in 1:num){
	n=length(spin.conf)
	uniformRV<-runif(n)
	for(i in 1:n){
		x<-trunc((index[i]-1)/size)+1
		y<-(index[i]-1)%%size+1
		counter = neighbours(spin.conf,x,y,size)# 计算黑色邻居数
		counter
		exp.top<-exp(beta*counter[[1]])
		exp.bottom<-exp(beta*counter[[1]])+exp(beta*counter[[2]])
		p=exp.top/exp.bottom
		if(uniformRV[i]<=p) spin.conf[x,y]<-1 else spin.conf[x,y]<-0 
	}
	samples[j,]=spin.conf[1:size^2]
}