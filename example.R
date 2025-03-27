library(flexsurv)
library(survival)
library(boot)
library(parallel)
library(dplyr)
library(tidyr)
library(backports)
library(discSurv)
library(withr)
library(latex2exp)

source("source_codes/function_source_code.r")

set.seed(4)

#///////////////////////////////////////////////////////////////////////////////
# 1 - Simulate Data 
#///////////////////////////////////////////////////////////////////////////////

load('simulated_data.Rdata')

svec = seq(0, 20, 1) ## time-points at which to evaluate survival

## approximate true survival curve by taking large sample from DGP 
## under sequence (1,1) while eliminating censoring 
temp = sim_data(n = 2000000,sim_true = T)
true_curve = sapply(svec, function(s) mean(temp$Y>s) )

#///////////////////////////////////////////////////////////////////////////////
# 2 - Summary Statistics
#///////////////////////////////////////////////////////////////////////////////

km = survfit(data=data, formula = Surv(Y, delta) ~ 1)
plot(km, xlim=c(0,30))


## distributions of event times times
print(km)

## time to course 2
survfit(data=data, formula = Surv(W1, next_trt) ~ 1) 
## time from course 2 to death
survfit(data=data, formula = Surv(W2, death2) ~ 1) 

png(filename = 'wtdist.png', width = 750, height = 750 )
par(mar=c(5,6,4,1)+.1)
hist(data$W1[data$next_trt], breaks=100, freq=F, ann=F,cex.axis=2)
title(xlab = TeX('Waiting Time, $W_1$'),
      ylab = 'Density',
      main = 'Distribution of Observed Waiting time to Course 2', 
      cex.lab=2, cex.main=2)

lines(density(data$W1[data$next_trt]), col='black', lwd=2)
dev.off()

## Summary Statistics table
n1 = nrow(data)
n2 = sum(data$next_trt)
n1
n2

table(data$A1)
table(data$A1) / n1
table(data$A2)
table(data$A2) / n2

table(data$L1)
table(data$L1) / n1
table(data$L2)
table(data$L2) / n2

#///////////////////////////////////////////////////////////////////////////////
# 3 - Continuous-Time Adjustment Models
#///////////////////////////////////////////////////////////////////////////////

###-- IPTW model adjusted for informative timing (adj=1)
## continuous-time IPTW adjustment implemented in run_iptw_boot

boot_ct_iptw_adj = boot(data, run_iptw_boot, sval=svec, adj=1, R=300)  

lwr_ct_adj = apply(boot_ct_iptw_adj$t, 2, quantile, p=.025)
upr_ct_adj = apply(boot_ct_iptw_adj$t, 2, quantile, p=.975)


###-- IPTW model *not* adjusted for informative timing (adj=2)

boot_ct_iptw_unadj = boot(data, run_iptw_boot, sval=svec, adj=2, R=300)  

lwr_ct_unadj = apply(boot_ct_iptw_unadj$t, 2, quantile, p=.025)
upr_ct_unadj = apply(boot_ct_iptw_unadj$t, 2, quantile, p=.975)

#///////////////////////////////////////////////////////////////////////////////
# 4 - Discrete-Time Adjustment Models
#///////////////////////////////////////////////////////////////////////////////

###-- adjusted (adj=1) discrete-time IPTW model,
## implemented in run_iptw_discrete_boot
boot_dt_iptw_adj = boot(data, run_iptw_discrete_boot, sval=svec, adj=1, R=300)  

lwr_dt_adj = apply(boot_dt_iptw_adj$t, 2, quantile, p=.025)
upr_dt_adj = apply(boot_dt_iptw_adj$t, 2, quantile, p=.975)

#///////////////////////////////////////////////////////////////////////////////
# 5 - Plot Estimated Curves
#///////////////////////////////////////////////////////////////////////////////

## plot true curve
png(filename = 'surv_plots.png', width = 750, height = 750 )

par(mar=c(8,6,4,1)+.1)

plot(svec, true_curve, col='red', pch=20, ylim=c(0,1), type='l',lty=2,ann = F, 
     axes = F, lwd=2 )

title(ylab = TeX('$P(T^{1,1}>t)$'), 
      main = "IPTW Survival Estimates Adjusted for Informative Timing",
      cex.main=2, cex.lab=2)

title(xlab = TeX('time, $t$$'),
      line = 6, cex.main=2, cex.lab=2)

## plot continuous-time adjusted IPTW estimate
points(svec, boot_ct_iptw_adj$t0, col='blue', pch=20, type='o', cex=2)
segments(svec, lwr_ct_adj, svec, upr_ct_adj, col='blue', lwd=2)

## plot continuous-time unadjusted IPTW estimate
points(svec+.1, boot_ct_iptw_unadj$t0, col='darkgreen', pch=20, type='o',cex=2)
#segments(svec+.1, lwr_ct_unadj, svec+.1, upr_ct_unadj, col='darkgreen')

## plot discrete-time adjusted IPTW estimate
points(svec+.2, c(1,boot_dt_iptw_adj$t0), col='black', pch=17, type='o',cex=2,)
segments(svec + .2, c(1,lwr_dt_adj), svec+.2, c(1,upr_dt_adj), col='black', pch=17, lwd=2)


legend('topright', 
       legend = c( 'Unadjusted Continuous-Time IPTW',
                   'Adjusted Continuous-Time IPTW',
                   'Adjusted Discrete-Time IPTW',
                   'True Curve'), 
       col=c('darkgreen','blue', 'black','red'), 
       lty=c(1,1,1,2), pch=c(20,20,20,NA), bty='n', cex=2, pt.cex = 2)

axis(2, at=seq(0, 1, .25), labels = seq(0, 1, .25), cex.axis=2)

atvec =seq(0,20,5)
nvec = sapply(atvec, function(x) sum(data$Y > x))
axis(1, at=atvec,padj = .65, 
     labels = paste0(atvec," \n (n=",nvec,")" ), cex.axis=2)
dev.off()