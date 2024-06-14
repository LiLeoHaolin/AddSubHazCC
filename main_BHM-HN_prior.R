###########################
### model specification ###
###########################

# observed data #
numGroups=10
x =c(2,0,1,6,7,3,5,1,0,3)
n =c(15,13,12,28,29,29,26,5,2,20)

# prior information #
mean.Mu=-1.735
perc.Mu=0.146
tau.prec=1

# analysis information #
targetResp=0.30
params = c("p")
n.iter=40000
n.burnin=10000

# directories #
setwd("C:/Users/Haolin Li/Desktop/intern/04_code/BHM-HN_prior")

################################
### analysis (DO NOT MODIFY) ###
################################

library(R2WinBUGS)
sink("mod1.txt")        
cat("model
 {
 for (i in 1:numGroups)
 {
 x[i] ~ dbin(p[i],n[i]);
 logit(p[i]) <- rho[i];
 rho[i] ~ dnorm(mu,prec.tau)
 pg[i] <- step(p[i]- targetResp)
 }
 mu ~ dnorm(mean.Mu, perc.Mu)
 tau ~ dnorm(0,tau.prec)I(0,)
 prec.tau<-pow(tau,-2)
 }", fill = TRUE)
sink()
data = list(x=x, n=n, numGroups=numGroups, targetResp=targetResp, mean.Mu=mean.Mu, perc.Mu=perc.Mu, tau.prec=tau.prec)
inits <- function () {list(mu=1, tau=0.1)}
bugs.out <- bugs(data=data, inits=inits, parameters.to.save=params, model.file="mod1.txt", n.chains=2, n.iter=n.iter, n.burnin=n.burnin, debug=F, DIC=TRUE, bugs.directory = "C:\\Program Files\\WinBUGS14", working.directory=getwd())
print(bugs.out, digits = 3)



