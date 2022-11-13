library(rms)
library(zoo)

nset = 1

path.dat <- paste0("/nas/longleaf/home/haolin/adane/s", nset, "/data/")
path.results <- paste0("/nas/longleaf/home/haolin/adane/s", nset, "/results/")

nsim = 500
seq = 1:nsim

beta.hat.cum = matrix(0, nrow = 2, ncol=nsim)
se.cum = matrix(0, nrow = 2, ncol=nsim)
cov.cum = matrix(0, nrow = 2, ncol=nsim)
beta.true = c(0.5, 0.5)

for (k in 1:nsim){
  
  ## necessary input ##
  
  cat(k)
  dat = read.csv(paste0(path.dat, 'dat', k, '.csv'))
  time = dat$obs.T
  delta = dat$delta
  delta.epsilon = dat$delta.epsilon
  cov = as.matrix(cbind(dat$Z1, dat$Z2))
  ccind = as.numeric(is.na(dat$Z1)==F)
  subind = dat$xi
  ncov = ncol(cov)
  ssize = nrow(dat)
  
  ## point estimate of regression coefficients ##
  
  time.pt = unique(sort(time))
  Y = matrix(0, nrow=ssize, ncol=length(time.pt))
  N = matrix(0, nrow=ssize, ncol=length(time.pt))
  dN = matrix(0, nrow=ssize, ncol=length(time.pt))
  
  for (i in 1:length(time.pt)){
    N[,i] = as.numeric((time <= time.pt[i])&(delta.epsilon==1))
    dN[,i] = as.numeric((time == time.pt[i])&(delta.epsilon==1))
    Y[,i] = 1-N[,i]
  }
  
  alpha.hat = (t(as.matrix(subind)*as.matrix(as.numeric(delta.epsilon!=1))) %*% Y) / (t(as.matrix(as.numeric(delta.epsilon!=1))) %*% Y)
  rho = as.matrix(as.numeric(delta.epsilon==1)) %*% matrix(1, nrow=1, ncol=length(time.pt)) + (as.matrix(as.numeric(delta.epsilon!=1))*as.matrix(subind)) %*% (as.matrix(1/alpha.hat))
  
  cen.fit <- survfit(Surv(time, (1-delta)) ~ 1)
  time.mid = rollmean(c(0,time.pt), 2)
  r = matrix(0, nrow=ssize, ncol=length(time.pt))
  r.mid = matrix(0, nrow=ssize, ncol=length(time.mid))
  
  G.hat.num = matrix(1, nrow=ssize, ncol=1) %*% t(as.matrix(summary(cen.fit, times = time.pt)$surv))
  G.hat.num.mid = matrix(1, nrow=ssize, ncol=1) %*% t(as.matrix(summary(cen.fit, times = time.mid)$surv))
  
  if (dim(G.hat.num)[2]<length(time.pt)){
    G.hat.num = cbind(G.hat.num, matrix(0, nrow = ssize, ncol = length(time.pt)-dim(G.hat.num)[2]))
    G.hat.num.mid = cbind(G.hat.num.mid, matrix(0, nrow = ssize, ncol = length(time.mid)-dim(G.hat.num.mid)[2]))
  }else{
    G.hat.num = G.hat.num
    G.hat.num.mid = G.hat.num.mid
  }
  
  G.hat.den = matrix(0, nrow=ssize, ncol=length(time.pt))
  G.hat.den.mid = matrix(0, nrow=ssize, ncol=length(time.pt))
  for (i in 1:length(time.pt)){
    r[,i] = as.numeric((delta==1)|((delta==0)&(time>=time.pt[i])))
    r.mid[,i] = as.numeric((delta==1)|((delta==0)&(time>=time.mid[i])))
    G.hat.den[,i] = summary(cen.fit, times = min(time.pt[i], time))$surv
    G.hat.den.mid[,i] = summary(cen.fit, times = min(time.mid[i], time))$surv
  }
  omega.mid = r.mid*G.hat.num.mid/G.hat.den.mid
  omega = r*G.hat.num/G.hat.den
  
  Zbar.c.mid = matrix(0, nrow=ncov, ncol=length(time.pt))
  Zbar.c = matrix(0, nrow=ncov, ncol=length(time.pt))
  Zbar.c[,-1] = Zbar.c[,-2]
  for (i in 1:ncov){
    Zbar.c.mid[i,] = colSums(rho*omega.mid*Y*(as.matrix(cov[,i]) %*% t(as.matrix(exp(-time.mid)))), na.rm=T)/colSums(rho*omega.mid*Y, na.rm=T)
    Zbar.c[i,] = colSums(rho*omega*Y*(as.matrix(cov[,i]) %*% t(as.matrix(exp(-time.pt)))), na.rm=T)/colSums(rho*omega*Y, na.rm=T)
  }
  
  time.diff = diff(c(0,time.pt))
  beta.num = matrix(0, nrow = ncov, ncol=1)
  for (i in 1:ncov){
    beta.num[i] = sum((((cov[,i] %*% t(as.matrix(exp(-time.pt)))) - (matrix(1, nrow=ssize, ncol=1) %*% Zbar.c[i,]))*rho*omega*dN) , na.rm=T)
  }
  
  beta.den = matrix(0, nrow = ncov, ncol=ncov)
  for (i in 1:ncov){
    for (j in 1:ncov){
      prelim = ((cov[,i] %*% t(as.matrix(exp(-time.mid))) - matrix(1, nrow=ssize, ncol=1) %*% Zbar.c.mid[i,])*(cov[,j] %*%t(as.matrix(exp(-time.mid))) - matrix(1, nrow=ssize, ncol=1) %*% Zbar.c.mid[j,])*rho*omega.mid*Y)
      prelim[is.na(prelim)==T] = 0
      beta.den[i,j] = sum((prelim %*% (as.matrix(time.diff))), na.rm=T)
      }
  }
  
  beta.hat =solve(beta.den) %*% beta.num
  
  # collect the point estimates
  
  beta.hat.cum[,k] = beta.hat
  
  ## variance estimate of regression coefficients ##
  
  Omega.hat = beta.den/ssize
  
  alpha.tilde = sum(subind)/ssize
  rho.tilde = (delta.epsilon == 1) + (delta.epsilon != 1)*subind/alpha.tilde
  
  # eta_i part
  
  dLambda1 = t(as.matrix(colSums(rho*omega*dN)/colSums(rho*omega*Y)))
  dLambda1[length(dLambda1)] = 0
  Z.beta.hat = as.matrix(cov %*% beta.hat)
  Z.t.beta.hat.mid = Z.beta.hat %*% t(as.matrix(exp(-time.mid)))
  Z.t.beta.hat = Z.beta.hat %*% t(as.matrix(exp(-time.pt)))
  dLambda2 = t(as.matrix(colSums(rho*omega.mid*Y*Z.t.beta.hat.mid,na.rm=T)/colSums(rho*omega.mid*Y)))
  
  eta.hat = matrix(0, nrow = ssize, ncol = ncov)
  for (i in 1:ncov){
    eta.hat[,i] = rowSums(omega*dN*(cov[,i] %*% t(as.matrix(exp(-time.pt))) - matrix(1, nrow=ssize, ncol=1) %*% Zbar.c[i,]), na.rm=T) - 
      rowSums(omega*Y*(matrix(1, nrow=ssize, ncol=1) %*% dLambda1)*(cov[,i] %*% t(as.matrix(exp(-time.pt))) - matrix(1, nrow=ssize, ncol=1) %*% Zbar.c[i,]), na.rm=T) +
      rowSums(omega.mid*Y*(matrix(1, nrow=ssize, ncol=1) %*% dLambda2)*(cov[,i] %*% t(as.matrix(exp(-time.mid))) - matrix(1, nrow=ssize, ncol=1) %*% Zbar.c.mid[i,])*(matrix(1, nrow=ssize, ncol=1) %*% t(as.matrix(time.diff))), na.rm=T) - 
      rowSums(omega.mid*Y*Z.t.beta.hat.mid*(cov[,i] %*% t(as.matrix(exp(-time.mid))) - matrix(1, nrow=ssize, ncol=1) %*% Zbar.c.mid[i,])*(matrix(1, nrow=ssize, ncol=1) %*% t(as.matrix(time.diff))), na.rm=T)
  }
  
  # psi_i part
  
  X.less.u = ((as.matrix(time) %*% matrix(1, nrow = 1, ncol=length(time.pt))) < (matrix(1, nrow = ssize, ncol=1) %*% t(as.matrix(time.pt))))
  X.geq.u = ((as.matrix(time) %*% matrix(1, nrow = 1, ncol=length(time.pt))) >= (matrix(1, nrow = ssize, ncol=1) %*% t(as.matrix(time.pt))))
  
  q1 = matrix(0, nrow = ncov, ncol = length(time.pt))
  for (i in 1:ncov){
    work1 = omega*dN*(cov[,i] %*% t(as.matrix(exp(-time.pt))) - matrix(1, nrow=ssize, ncol=1) %*% Zbar.c[i,]) -
      omega*Y*(matrix(1, nrow=ssize, ncol=1) %*% dLambda1)*(cov[,i] %*% t(as.matrix(exp(-time.pt))) - matrix(1, nrow=ssize, ncol=1) %*% Zbar.c[i,])
    work2 = omega.mid*Y*(matrix(1, nrow=ssize, ncol=1) %*% dLambda2)*(cov[,i] %*% t(as.matrix(exp(-time.mid))) - matrix(1, nrow=ssize, ncol=1) %*% Zbar.c.mid[i,])*(matrix(1, nrow=ssize, ncol=1) %*% t(as.matrix(time.diff)))-
      omega.mid*Y*Z.t.beta.hat.mid*(cov[,i] %*% t(as.matrix(exp(-time.mid))) - matrix(1, nrow=ssize, ncol=1) %*% Zbar.c.mid[i,])*(matrix(1, nrow=ssize, ncol=1) %*% t(as.matrix(time.diff)))
    work1[is.na(work1)] <- 0
    work2[is.na(work2)] <- 0
    work2 = cbind(work2[,-1], matrix(0, nrow=ssize, ncol=1))
    work = (work1+work2)
    work <- t(apply(work, 1, rev))
    work = t(apply(work, 1, cumsum))
    work <- t(apply(work, 1, rev))
    q1[i,] = -colSums(X.less.u*work)/ssize
  }
  
  pi.hat = colSums(rho*X.geq.u)/ssize
  
  dNc = matrix(0, nrow=ssize, ncol=length(time.pt))
  for (i in 1:length(time.pt)){
    dNc[,i] = as.numeric((time == time.pt[i])&(delta.epsilon==0))
  }
  dLambdac = colSums(rho*dNc)/colSums(rho*X.geq.u)
  
  psi.hat = matrix(0, nrow = ssize, ncol = ncov)
  for (i in 1:ncov){
    rat = q1[i,]/pi.hat
    rat[rat==Inf] = 0
    rat[rat==-Inf] = 0
    rat[is.na(rat)==T] = 0
    psi.hat[,i] = rowSums((matrix(1, nrow=ssize, ncol=1)%*%(rat))*dNc, na.rm=T) -
     rowSums((matrix(1, nrow=ssize, ncol=1)%*%(rat*dLambdac))*X.geq.u, na.rm=T)
  }
  
  # mu_i part
  
  E.hat.denom = colSums((as.matrix(subind)*(1-(delta.epsilon==1))/alpha.tilde) %*% matrix(1, nrow = 1, ncol = length(time.pt))*Y)/ssize
  E.hat.num = matrix(0, nrow = ncov, ncol = length(time.pt))
  for (i in 1:ncov){
    E.hat.num[i,] = colSums((as.matrix(subind)*(1-(delta.epsilon==1))/alpha.tilde)%*% matrix(1, nrow = 1, ncol = length(time.pt))*
      (cov[,i] %*% t(as.matrix(exp(-time.mid))) - matrix(1, nrow=ssize, ncol=1) %*% Zbar.c.mid[i,])*omega.mid*Y*Z.t.beta.hat.mid, na.rm = T)/ssize
  }
  
  mu.hat.3 = matrix(0, nrow = ssize, ncol = ncov)
  for (i in 1:ncov){
    mu.hat.3[,i] = rowSums((1-(delta.epsilon==1)) %*% matrix(1, nrow = 1, ncol = length(time.pt)) *Y * (matrix(1, nrow = ssize, ncol = 1) %*% (E.hat.num[i,]/E.hat.denom)) * (matrix(1, nrow=ssize, ncol=1) %*% t(as.matrix(time.diff)))) 
  }
  
  mu.hat.1.and.2 = matrix(0, nrow = ssize, ncol = ncov)
  for (i in 1:ncov){
    mu.hat.1.and.2[,i] = rowSums(((1-(delta.epsilon==1)) %*% matrix(1, nrow = 1, ncol = length(time.pt)))*omega*Y*(matrix(1, nrow=ssize, ncol=1) %*% dLambda1)*(cov[,i] %*% t(as.matrix(exp(-time.pt))) - matrix(1, nrow=ssize, ncol=1) %*% Zbar.c[i,]), na.rm=T) -
      rowSums(((1-(delta.epsilon==1)) %*% matrix(1, nrow = 1, ncol = length(time.pt)))*omega.mid*Y*(matrix(1, nrow=ssize, ncol=1) %*% dLambda2)*(cov[,i] %*% t(as.matrix(exp(-time.mid))) - matrix(1, nrow=ssize, ncol=1) %*% Zbar.c.mid[i,])*(matrix(1, nrow=ssize, ncol=1) %*% t(as.matrix(time.diff))), na.rm=T) + 
      rowSums(((1-(delta.epsilon==1)) %*% matrix(1, nrow = 1, ncol = length(time.pt)))*omega.mid*Y*Z.t.beta.hat.mid*(cov[,i] %*% t(as.matrix(exp(-time.mid))) - matrix(1, nrow=ssize, ncol=1) %*% Zbar.c.mid[i,])*(matrix(1, nrow=ssize, ncol=1) %*% t(as.matrix(time.diff))), na.rm=T)
  }
  
  mu.hat = mu.hat.1.and.2-mu.hat.3
  
  # expression for Sigma.hat
  
  Sigma.hat = matrix(0, nrow = ncov, ncol = ncov)
  for (i in 1:ncov){
    for (j in 1:ncov){
      Sigma.hat[i,j] = (sum((as.matrix(eta.hat[,i]+psi.hat[,i])) * as.matrix(eta.hat[,j]+psi.hat[,j]) *as.matrix(rho.tilde)) +
        (1-alpha.tilde)/alpha.tilde * sum((as.matrix(mu.hat[,i])) * as.matrix(mu.hat[,j]) *as.matrix(rho.tilde)))/ssize
    }
  }
  
  Var = solve(Omega.hat)%*% Sigma.hat %*% solve(Omega.hat)/ssize
  
  # collect the variance estimates and coverages
  
  se.cum[,k] = sqrt(diag(Var))
  cov.cum[,k] = (beta.hat+1.96*sqrt(diag(Var)) >= beta.true)&(beta.hat-1.96*sqrt(diag(Var)) <= beta.true)
  
  # point estimation of baseline cumulative subdistribution hazard
  
  Lambda0 = dLambda1 - dLambda2*t(as.matrix(time.diff))
  Lambda10 = t(apply(Lambda0, 1, cumsum))
  
  # variance estimate of estimated baseline cumulative subdistribution hazard
  
  e.hat = (matrix(1, nrow = ncov, ncol = 1) %*% time.diff) * Zbar.c.mid
  e.hat = -t(apply(e.hat, 1, cumsum))
  
  Y.bar = colSums(rho*omega*Y)/ssize
  Y.bar.mid = colSums(rho*omega.mid*Y)/ssize
  
  # W1 component
  
  W1.2 = (omega*dN*(matrix(1, nrow=ssize, ncol=1) %*% (1/Y.bar))) - 
    (omega*Y*(matrix(1, nrow=ssize, ncol=1) %*% dLambda1)*(matrix(1, nrow=ssize, ncol=1) %*% (1/Y.bar))) +
    (omega.mid*Y*(matrix(1, nrow=ssize, ncol=1) %*% dLambda2)*(matrix(1, nrow=ssize, ncol=1) %*% (1/Y.bar.mid))*(matrix(1, nrow=ssize, ncol=1) %*% t(as.matrix(time.diff)))) - 
    (omega.mid*Y*Z.t.beta.hat.mid*(matrix(1, nrow=ssize, ncol=1) %*% (1/Y.bar.mid))*(matrix(1, nrow=ssize, ncol=1) %*% t(as.matrix(time.diff))))
  W1.2 = t(apply(W1.2, 1, cumsum))
  
  W1.3 = matrix(0, nrow = ssize, ncol = length(time.pt))
  
  q2 = matrix(0, nrow = 1, ncol = length(time.pt))
  work1 = omega*dN*(matrix(1, nrow=ssize, ncol=1) %*% (1/Y.bar)) -
    omega*Y*(matrix(1, nrow=ssize, ncol=1) %*% dLambda1)*(matrix(1, nrow=ssize, ncol=1) %*% (1/Y.bar))
  work2 = omega.mid*Y*(matrix(1, nrow=ssize, ncol=1) %*% dLambda2)*(matrix(1, nrow=ssize, ncol=1) %*% (1/Y.bar.mid))*(matrix(1, nrow=ssize, ncol=1) %*% t(as.matrix(time.diff)))-
    omega.mid*Y*Z.t.beta.hat.mid*(matrix(1, nrow=ssize, ncol=1) %*% (1/Y.bar.mid))*(matrix(1, nrow=ssize, ncol=1) %*% t(as.matrix(time.diff)))
  work1[is.na(work1)] <- 0
  work2[is.na(work2)] <- 0
  work2 = cbind(work2[,-1], matrix(0, nrow=ssize, ncol=1))
  work = (work1+work2)
  
  for (j in 1:length(time.pt)){
    work.int = cbind(matrix(1, nrow = ssize, ncol = j), matrix(0, nrow = ssize, ncol = length(time.pt)-j))
    work.new = work*work.int
    work.new <- t(apply(work, 1, rev))
    work.new = t(apply(work, 1, cumsum))
    work.new <- t(apply(work, 1, rev))
    q2 = -colSums(X.less.u*work.new)/ssize
    rat2 = q2/pi.hat
    rat2[rat2==Inf] = 0
    rat2[rat2==-Inf] = 0
    rat2[is.na(rat2)==T] = 0
    W1.3[,j] = rowSums((matrix(1, nrow=ssize, ncol=1)%*%(rat2))*dNc) -
      rowSums((matrix(1, nrow=ssize, ncol=1)%*%(rat2*dLambdac))*X.geq.u)
  }
  
  W1 = (eta.hat+psi.hat)%*%solve(Omega.hat)%*%e.hat + W1.2 + W1.3
  
  # W2 component
  
  E.hat.num.new = colSums((as.matrix(subind)*(1-(delta.epsilon==1))/alpha.tilde) %*% matrix(1, nrow = 1, ncol = length(time.pt))*omega.mid*Y*Z.t.beta.hat.mid, na.rm = T)/ssize
  b.hat = omega.mid*Y*Z.t.beta.hat.mid - Y*(matrix(1, nrow = ssize, ncol=1)%*%(E.hat.num.new/E.hat.denom))
  
  W2.2 = b.hat*(matrix(1, nrow=ssize, ncol=1) %*% (1/Y.bar.mid))*(matrix(1, nrow=ssize, ncol=1) %*% t(as.matrix(time.diff)))
  W2.2 = t(apply(W2.2, 1, cumsum))
  W2.2 = ((1-(delta.epsilon==1))%*%matrix(1, nrow = 1, ncol = length(time.pt)))*W2.2
  
  W2 = (mu.hat)%*%solve(Omega.hat)%*%e.hat + W2.2
  
  # expression for Sigma hat of estimated baseline cumulative subdistribution hazard
  
  Sigma.hat.cum.sub.h = (colSums((rho.tilde %*% matrix(1, nrow=1, ncol=length(time.pt)))*W1*W1, na.rm = T) + (1-alpha.tilde)/alpha.tilde * colSums((rho.tilde %*% matrix(1, nrow=1, ncol=length(time.pt)))*W2*W2, na.rm = T))/ssize
  
  upper = Lambda10 + 1.96*sqrt(Sigma.hat.cum.sub.h/ssize)
  lower = Lambda10 - 1.96*sqrt(Sigma.hat.cum.sub.h/ssize)
  
  info.cum.sub.h = data.frame(time = time.pt, est = as.numeric(Lambda10), lower = as.numeric(lower), upper = as.numeric(upper))
  
  write.csv(info.cum.sub.h, file=paste0(path.results, 'info_cum_sub_h_', k, '.csv'), row.names = F)
  
}

write.csv(beta.hat.cum, file=paste0(path.results, 'beta_hat_cum_', nset, '.csv'), row.names = F)
write.csv(se.cum, file=paste0(path.results, 'se_cum_', nset, '.csv'), row.names = F)
write.csv(cov.cum, file=paste0(path.results, 'cov_cum_', nset, '.csv'), row.names = F)







