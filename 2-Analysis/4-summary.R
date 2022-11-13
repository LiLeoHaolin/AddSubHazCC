# these need to be changed:
nset = 1
tau = 0.12
true = c(0.5, 0.5)

# set up the working directory and number of simulations
nsim = 500
path.results = paste0('C:/Users/Haolin Li/Desktop/adane/02_simulation/02_results/s', nset, '/')

# summary for regression coefficients

est = read.csv(paste0(path.results, 'beta_hat_cum_', nset, '.csv'))
se = read.csv(paste0(path.results, 'se_cum_', nset, '.csv'))
cov = read.csv(paste0(path.results, 'cov_cum_', nset, '.csv'))

rowMeans(est, na.rm = T) -true
rowMeans(se, na.rm = T)
sd(est[1,], na.rm = T)
sd(est[2,], na.rm = T)
rowMeans(cov, na.rm = T)

# summary for baseline cumulative subdistribution hazard

t1 = tau/3
t2 = 2*tau/3
Lambda1 = -log(1-0.4*(1-exp(-t1)))
Lambda2 = -log(1-0.4*(1-exp(-t2)))

est.cum = matrix(NA, nrow = 2, ncol = nsim)
se.cum = matrix(NA, nrow = 2, ncol = nsim)
cov.cum = matrix(NA, nrow = 2, ncol = nsim)

for (i in 1:nsim){
  try({
    result = read.csv(paste0(path.results, 'info_cum_sub_h_', i, '.csv'))
    
    result1 = result[result$time <= t1, ]
    result2 = result[result$time <= t2, ]
    
    est.cum[1, i] = result1$est[nrow(result1)]
    est.cum[2, i] = result2$est[nrow(result2)]
    
    cov.cum[1, i] = as.numeric((result1$lower[nrow(result1)]<=Lambda1)&(result1$upper[nrow(result1)]>=Lambda1))
    cov.cum[2, i] = as.numeric((result2$lower[nrow(result2)]<=Lambda2)&(result2$upper[nrow(result2)]>=Lambda2))
    
    se.cum[1, i] = (result1$est[nrow(result1)] - result1$lower[nrow(result1)])/1.96
    se.cum[2, i] = (result2$est[nrow(result2)] - result2$lower[nrow(result2)])/1.96
  })
}

rowMeans(est.cum, na.rm = T) - c(Lambda1, Lambda2)
rowMeans(se.cum, na.rm = T)
sd(est.cum[1,], na.rm = T)
sd(est.cum[2,], na.rm = T)
rowMeans(cov.cum, na.rm = T)

# summary for covariate specific CIF

t1 = tau/3
t2 = 2*tau/3
Lambda1 = -log((1-0.4*(1-exp(-t1)))*exp(-(1-exp(-t1))))
Lambda2 = -log((1-0.4*(1-exp(-t2)))*exp(-(1-exp(-t2))))
#Lambda1 = -log((1-0.4*(1-exp(-t1)))*exp(-(1-exp(-1.75*t1))))
#Lambda2 = -log((1-0.4*(1-exp(-t2)))*exp(-(1-exp(-1.75*t2))))

CIF1 = 1-exp(-Lambda1)
CIF2 = 1-exp(-Lambda2)

est.cum = matrix(NA, nrow = 2, ncol = nsim)
se.cum = matrix(NA, nrow = 2, ncol = nsim)
cov.cum = matrix(NA, nrow = 2, ncol = nsim)

for (i in 1:nsim){
  try({
    result = read.csv(paste0(path.results, 'info_csh_', i, '.csv'))
    
    result1 = result[result$time <= t1, ]
    result2 = result[result$time <= t2, ]
    
    est.cum[1, i] = 1-exp(-result1$CIF[nrow(result1)])
    est.cum[2, i] = 1-exp(-result2$CIF[nrow(result2)])
    
    se.cum[1, i] = exp(-2*(1-exp(-result1$CIF[nrow(result1)])))*(result1$CIF[nrow(result1)] - result1$lower[nrow(result1)])/1.96
    se.cum[2, i] = exp(-2*(1-exp(-result1$CIF[nrow(result1)])))*(result2$CIF[nrow(result2)] - result2$lower[nrow(result2)])/1.96
    
    cov.cum[1, i] = as.numeric((est.cum[1, i]-1.96*se.cum[1, i]<=CIF1)&(est.cum[1, i]+1.96*se.cum[1, i]>=CIF1))
    cov.cum[2, i] = as.numeric((est.cum[2, i]-1.96*se.cum[2, i]<=CIF2)&(est.cum[2, i]+1.96*se.cum[2, i]>=CIF2))
  })
}

rowMeans(est.cum, na.rm = T) - c(Lambda1, Lambda2)
rowMeans(se.cum, na.rm = T)
sd(est.cum[1,], na.rm = T)
sd(est.cum[2,], na.rm = T)
rowMeans(cov.cum, na.rm = T)

# summary for covariate-specific confidence band

cov.EP.un = matrix(NA, nrow = 1, ncol = nsim)
cov.EP.log = matrix(NA, nrow = 1, ncol = nsim)
cov.EP.loglog = matrix(NA, nrow = 1, ncol = nsim)
cov.HW.un = matrix(NA, nrow = 1, ncol = nsim)
cov.HW.log = matrix(NA, nrow = 1, ncol = nsim)
cov.HW.loglog = matrix(NA, nrow = 1, ncol = nsim)

for (i in 1:nsim){
  try({
    result = read.csv(paste0(path.results, 'info_CB_', i, '.csv'))
    
    #result$true = -log((1-0.4*(1-exp(-result$time)))*exp(-(1-exp(-1.75*result$time))))
    result$true = -log((1-0.4*(1-exp(-result$time)))*exp(-(1-exp(-result$time))))
    
    result$cov.EP.un = as.numeric((result$CB.EP.un.lb<=result$true)&(result$CB.EP.un.ub>= result$true))
    result$cov.EP.log = as.numeric((result$CB.EP.log.lb<=result$true)&(result$CB.EP.log.ub>= result$true))
    result$cov.EP.loglog = as.numeric((result$CB.EP.loglog.lb<=result$true)&(result$CB.EP.loglog.ub>= result$true))
    
    result$cov.HW.un = as.numeric((result$CB.HW.un.lb<=result$true)&(result$CB.HW.un.ub>= result$true))
    result$cov.HW.log = as.numeric((result$CB.HW.log.lb<=result$true)&(result$CB.HW.log.ub>= result$true))
    result$cov.HW.loglog = as.numeric((result$CB.HW.loglog.lb<=result$true)&(result$CB.HW.loglog.ub>= result$true))
    
    cov.EP.un[i] = prod(result$cov.EP.un)
    cov.EP.log[i] = prod(result$cov.EP.log)
    cov.EP.loglog[i] = prod(result$cov.EP.loglog)
    
    cov.HW.un[i] = prod(result$cov.HW.un)
    cov.HW.log[i] = prod(result$cov.HW.log)
    cov.HW.loglog[i] = prod(result$cov.HW.loglog)
  })
}

mean(cov.EP.un, na.rm = T)
mean(cov.EP.log, na.rm = T)
mean(cov.EP.loglog, na.rm = T)
mean(cov.HW.un, na.rm = T)
mean(cov.HW.log, na.rm = T)
mean(cov.HW.loglog, na.rm = T)





