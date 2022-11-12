
## prerequisites ##

library(rootSolve)
path.dat <- "/nas/longleaf/home/haolin/adane/s8/data/"
path.results <- "/nas/longleaf/home/haolin/adane/s8/results/"

## simulation setting ##

n=2000
ntilde = 364
beta1 = matrix(c(1, 0.5), nrow=2)
beta2 = matrix(c(-0.5, 0.5), nrow=2)
prob=0.4
nsim=500
event.prop = rep(0, nsim)
censor.prop = rep(0, nsim)
ratio = rep(0, nsim)

## simulation ##

for (j in 1:nsim){
  cat(j)
  
  # covariates
  Z1 = rbinom(n, 1, 0.5)
  Z2 = runif(n, 1, 2)
  Z = cbind(Z1, Z2)
  
  # event type
  prob1 = 1-(1-prob)*exp(-Z %*% beta1)
  epsilon = matrix(0, nrow = n)
  for (i in 1:n){
    epsilon[i] = sample(c(1:2), 1, prob = c(prob1[i],1-prob1[i]))
  }
  
  # survival time
  unif1 = matrix(runif(n), nrow = n)
  unif2 = matrix(runif(n), nrow = n)
  T1 = matrix(0, nrow = n)
  T2 = matrix(0, nrow = n)
  for (i in 1:n){
    fun1 <- function(t){
      (1-(1-prob*(1-exp(-t)))*exp(as.numeric(-Z[i,]%*%beta1)*(1-exp(-t))))/(1-(1-prob)*exp(as.numeric(-Z[i,]%*%beta1)))-unif1[i]
    } 
    fun2 <- function(t){
      1-exp(-t-as.numeric(Z[i,]%*%beta2)*(1-exp(-t)))-unif2[i]
    }
    T1[i] = uniroot(fun1, c(0,100))$root
    T2[i] = uniroot(fun2, c(0,100))$root
  }
  Tlat = matrix(0, nrow = n)
  for (i in 1:n){
    if(epsilon[i] == 1){
      Tlat[i] = T1[i]
    }else if (epsilon[i] == 2){
      Tlat[i] = T2[i]
    }
  }
  
  # censoring
  C = runif(n, min = 0, max = 0.11)
  min = 0.105
  
  # observed data
  obs.T = matrix(0, nrow = n)
  for (i in 1:n){
    obs.T[i] = min(Tlat[i], C[i], min)
  }
  delta = as.numeric(obs.T == Tlat)
  delta.epsilon = epsilon*delta
  subid = 1:n
  dat = data.frame(subid, Z1, Z2, obs.T, delta, delta.epsilon)
  
  # event and censoring proportion
  event.prop[j] = table(delta.epsilon)[2]/n
  censor.prop[j] = (n-sum(delta))/n
  
  # take case-cohort samples
  ccind = sample(c(1:n), ntilde, replace = F)
  dat$xi = as.numeric(dat$subid %in% ccind)
  for (i in 1:n){
    if ((dat$subid[i] %in% ccind)|(dat$delta.epsilon[i] ==1)){
      dat$Z1[i] = dat$Z1[i]
      dat$Z2[i] = dat$Z2[i]
    } else{
      dat$Z1[i] = NA
      dat$Z2[i] = NA
    }
  }
  
  # check ratios
  ratio[j] = table(dat[(is.na(dat$Z1)==F),]$delta.epsilon)[2]/(table(dat[(is.na(dat$Z1)==F),]$delta.epsilon)[1]+table(dat[(is.na(dat$Z1)==F),]$delta.epsilon)[3])
  
  # write data set in csv
  write.csv(dat, file=paste0(path.dat, 'dat', j, '.csv'), row.names = F)
}

mean(event.prop)
mean(censor.prop)
mean(ratio)







