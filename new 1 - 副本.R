library(igraph)
library(MASS)
library(Matrix)

config.samples_num = 200 # 样本数
config.node_num = 10 # 节点数q
config.co_num = 20 # 属性信息个数p
config.power = 1 # 生成scale-free network 幂率指数
config.beta = 4 # signal strength
config.rho = 0.8 # proportion of relevant covariates
G = vector('list', config.samples_num) # 图list
M = vector('list', config.samples_num) # 矩阵list
Th = vector('list', config.samples_num) # theta矩阵list

GenerateVectorTheta = function( config ){
    len = config.co_num +1
    v_theta = rep(0, len)
    rand = runif(len, min = 0, max = 1)
    for (j in 1:len){
        if( rand[j]>config.rho_now ){
            v_theta[j] = 0
        }else if( rand[j]<= config.rho_now/2 ){
            v_theta[j] = config.beta
        }else{
            v_theta[j] = (-config.beta)
        }
    }
    return (v_theta)
}

GibbsSampling = function( config ,theta ){
    y = matrix(rep(rep(0,config.node_num),config.iter),nrow = config.iter, ncol = config.node_num)
    for (i in 1:config.iter){
        for (j in 1:config.node_num){
            uniformRV<-runif(1)
            # if(y[i,j] == 1){
            #     A = theta[j,] %*% t(y[i,])
            # }else{
            #     A = theta[j,] %*% t(y[i,]) + theta[j,j]
            # }
            # exp_A = exp(A)
            if( uniformRV < config.p_theta ){
                y[i,j] = 1
            }else{
                y[i,j] = 0
            }
        }
    }
    return (y)
}

if( config.rho >0 && config.rho <1 && config.p_theta >0 && config.p_theta<1 ){
    for (i in 1:config.samples_num){
        g <- barabasi.game(config.node_num , power = config.power, directed = FALSE) #scale-free network 生成图
        G[[i]] = g
        M[[i]] = as_adjacency_matrix(g)
        Th[[i]] = M[[i]]

        theta = matrix( rep(0, config.node_num*config.node_num), config.node_num, config.node_num )
        for (ir in 1:config.node_num){
            for (ic in ir:config.node_num){
                if( ir == ic){
                    config.rho_now = 1
                }else{
                    config.rho_now = config.rho
                }
                x = rnorm ( config.co_num )
                x = cbind (1, t(x))
                v_theta = GenerateVectorTheta(config)
                if( ir == ic){
                    theta[ir,ic] = (x %*% v_theta)/2
                }else{
                    theta[ir,ic] = (x %*% v_theta)
                }
            }
        }
        theta = theta + t(theta)
        theta_sparse =  matrix( rep(0, config.node_num*config.node_num), config.node_num, config.node_num )
        info = summary(Th[[i]])
        for (io in 1:nrow(info)){
            Th[[i]]@x[io] = theta[info[io,1],info[io,2]]
            for (ir in 1:config.node_num){
                for (ic in 1:config.node_num){
                    if( ir == info[io,1] && ic == info[io,2]){
                        theta_sparse[ir,ic] = theta[info[io,1],info[io,2]]
                    }
                }
            }
        }
        for (ir in 1:config.node_num){
            theta_sparse[ir,ir] = theta[ir,ir]
        }

        y_binary = GibbsSampling(config, theta_sparse)

    }
}else{
    cat("rho is error!\n");
}

