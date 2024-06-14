###########################
### model specification ###
###########################

# observed data #
numGroups=10
x =c(2,0,1,6,7,3,5,1,0,3)
n =c(15,13,12,28,29,29,26,5,2,20)

# prior information #
pex = 0.5
pnex = 0.5
ex.mu.mean = -1.735
ex.mu.prec=0.146
HN.scale = 1
nex.mean=-1.734
nex.prec=0.128

# analysis information #
targetResp=0.30
params = c("p")
n.iter=40000
n.burnin=10000

# directories #
setwd("C:/Users/Haolin Li/Desktop/intern/04_code/EXNEX")

################################
### analysis (DO NOT MODIFY) ###
################################

library(R2WinBUGS)
sink("mod1.txt")        
cat("model{
  for(jj in 1:Nexch){
    mu[jj]~dnorm(mu.mean[jj],mu.prec[jj])
    prior.tau.prec[jj]<-pow(tau.HN.scale[jj],-2)
    tau[jj] ~ dnorm(0,prior.tau.prec[jj])I(0,)
    prec.tau[jj]<-pow(tau[jj],-2)
  }
  for(jj in 1:Nexch){
    for(j in 1:Nstrata){
      re[jj,j]~dnorm(0,prec.tau[jj])
      LogOdds[jj,j]<-mu[jj]+re[jj,j]
    }
  }
  for(j in 1:Nstrata) {
    LogOdds[Nmix,j]~dnorm(nex.mean,nex.prec)
  }
  for(j in 1:Nstrata){
    exch.index[j]~dcat(pMix[1:Nmix])
    for(jj in 1:Nmix) {
      exch[j,jj]<-equals(exch.index[j],jj)
    }
  }
  for(j in 1:Nstrata) {
    theta[j]<-LogOdds[exch.index[j],j]
  }
  for(i in 1:Nstrata) {
    logit(p[i])<-theta[i]
    pg[i]<-step(p[i]-p.cut)
    r[i]~dbin(p[i],n[i])
  }
}", fill = TRUE)
sink()
data = list(Nexch=2,Nmix=3, pMix=c(pex,0,pnex), Nstrata=numGroups, n=n, r=x, mu.mean=c(ex.mu.mean,0),mu.prec=c(ex.mu.prec,1), tau.HN.scale=c(HN.scale,1), nex.mean=nex.mean,nex.prec=nex.prec, p.cut=targetResp)
inits <- function () {list(mu=c(-0.5,0.5))}
bugs.out <- bugs(data=data, inits=inits, parameters.to.save=params, model.file="mod1.txt", n.chains=2, n.iter=n.iter, n.burnin=n.burnin,debug=F, DIC=F, bugs.directory = "C:\\Program Files\\WinBUGS14", working.directory=getwd())
print(bugs.out, digits = 3)




