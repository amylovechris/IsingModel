library(Matrix)
config.mean = 0;
config.sd = 1;
config.number = 100;
config.inter = 2000;
config.samples = 1000;

X = rnorm(config.number, config.mean, config.sd );
Z = matrix( rep(0, (config.inter+1)*config.number), (config.inter+1), config.number )
Mu = matrix( rep(0, config.inter+1), 1, config.inter+1 )
 
Mu[1,1] = 1;
 
for (i in 1:(config.inter+1)){
    current_mean = config.sd * sum(X) / (config.number + 1);
    current_sd = sqrt(config.sd / (config.number +1));
    Mu[1,i] = rnorm(1, current_mean , current_sd);
}
 
Y = Mu[,(config.inter-config.samples+1):config.inter];
hist(Y) 