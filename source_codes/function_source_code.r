#------------------------------------------------------------------------------#
## Functions for Generating Synthetic Data  ----
#------------------------------------------------------------------------------#
sim_data = function(n = 1000, sim_true=F){
  
  L1 = rbinom(n, 1, .5)
  if(sim_true){
    A1 = rep(1, n)
  }else{
    A1 = rbinom(n, 1, plogis(1 - 1*L1))
  }
  
  
  W1_next_trt = rexp(n = n, rate = exp(-3 + A1 + L1) )
  W1_death = rexp(n = n, rate = exp(-3 - A1 + L1) )
  
  if(sim_true){
    C1 = rep(Inf, n)
  }else{
    C1 = rexp(n = n, rate = exp(-4 - A1 + L1) )  
  }
  
  
  W1 = pmin(W1_next_trt, W1_death, C1)
  next_trt = W1_next_trt < pmin( W1_death, C1)
  death1  = W1_death < pmin(C1, W1_next_trt)
  
  L2 = rbinom(n, 1, .5 + .15*L1 - .15*(W1>15) )
  if(sim_true){
    A2 = rep(1, n)
  }else{
    A2 = rbinom(n, 1, plogis( 1 - 1*L2 + 1.5*(W1>15) ) )
  }
  
  if(sim_true){
    C2 = rep(Inf, n)
  }else{
    C2 = rexp(n = n, rate = exp( -4 + A2 + L2 - 2*(W1>15) ) )
  }
  
  W2_death = rexp(n = n, rate = exp( -3 + A2 + L2 - 2*(W1>15) ) )
  
  W2 = pmin(W2_death,C2)
  death2 = W2_death < C2
  
  Y = W1 + W2*next_trt
  
  d = data.frame(L1=L1, A1=A1, W1=W1, 
                 next_trt=next_trt, death1 = death1,
                 L2=L2, A2=A2, W2=W2, death2=death2,
                 Y=Y, delta = ifelse(death1 | death2,1,0))
  
  d$L2[!next_trt] = NA
  d$A2[!next_trt] = NA
  d$W2[!next_trt] = NA
  d$death2[!next_trt] = NA
  
  return(d)
}

#------------------------------------------------------------------------------#
## Continuous-Time Methods  ----
#------------------------------------------------------------------------------#

run_iptw = function(s, d, adj=1){
  d$W1_15 = 1*(d$W1>15)
  ## time point 1
  a1_logit = glm(data=d, formula = A1 ~ L1, family=binomial('logit'))
  d$pa1_hat = predict(a1_logit, d, 'response')  
  
  d$cen1 = ifelse(!d$next_trt & !d$death1, 1, 0)
  d$cen2[d$next_trt] = 1 - d$death2[d$next_trt]
  
  c1_formula = survival::Surv(W1, cen1) ~ A1 + L1
  c1_reg = flexsurv::flexsurvreg(data=d, 
                                 formula = c1_formula, 
                                 dist='exponential')
  
  eta = cbind(1, d$A1, d$L1) %*% c1_reg$coefficients
  d$pc1_hat = exp(-1*exp(eta)*d$W1 )
  
  ## time point 2
  if(adj==1){
    ## stage 2 propensity score model
    a2_logit = glm(data=d[d$next_trt==1, ], formula = A2 ~ L2 + W1_15, family=binomial('logit'))
    
    ## stage 2 censoring model
    c2_formula = survival::Surv(W2, cen2) ~ A2 + L2 + W1_15
    c2_reg = flexsurv::flexsurvreg(data=d[d$next_trt, ], formula = c2_formula, dist='exponential')
    
    eta = cbind(1, d$A2, d$L2, d$W1_15 ) %*% c2_reg$coefficients
    d$pc2_hat = exp(-1*exp(eta)*d$W2 )
    d$pc2_hat[is.na(d$pch2_hat)] = 0
    
  }else if(adj==2){
    ## stage 2 propensity score model
    a2_logit = glm(data=d[d$next_trt==1, ], formula = A2 ~ L2, family=binomial('logit'))
    
    ## stage 2 censoring model
    c2_formula = survival::Surv(W2, cen2) ~ A2 + L2
    c2_reg = flexsurv::flexsurvreg(data=d[d$next_trt, ], formula = c2_formula, dist='exponential')
    
    eta = cbind(1, d$A2, d$L2) %*% c2_reg$coefficients
    d$pc2_hat = exp(-1*exp(eta)*d$W2 )
    d$pc2_hat[is.na(d$pch2_hat)] = 0
  }else if(adj==3){
    a2_logit = glm(data=d[d$next_trt==1, ], formula = A2 ~ 1, family=binomial('logit'))
    
    c2_formula = survival::Surv(W2, cen2) ~ A2 
    c2_reg = flexsurv::flexsurvreg(data=d[d$next_trt, ], formula = c2_formula, dist='exponential')
    
    eta = cbind(1, d$A2) %*% c2_reg$coefficients
    d$pc2_hat = exp(-1*exp(eta)*d$W2 )
    d$pc2_hat[is.na(d$pch2_hat)] = 0
  }
  
  d$pa2_hat = predict(a2_logit, d, 'response')  
  
  d$pscore = ifelse(d$next_trt==0, d$pa1_hat, d$pa1_hat*d$pa2_hat)
  
  d$cscore = ifelse(d$next_trt==0, d$pc1_hat, d$pc1_hat*d$pc2_hat)
  
  d$wght = 1 / ( d$pscore*d$cscore )
  
  EY_num = sum( (d$Y[d$next_trt==0]>s)*d$A1[d$next_trt==0]*(1-d$cen1[d$next_trt==0])*d$wght[d$next_trt==0]) + 
    sum( (d$Y[d$next_trt==1]>s)*d$A1[d$next_trt==1]*d$A2[d$next_trt==1]*(1-d$cen2[d$next_trt==1])*d$wght[d$next_trt==1])
  
  EY_denom = sum( d$A1[d$next_trt==0]*(1-d$cen1[d$next_trt==0])*d$wght[d$next_trt==0] ) + 
    sum( d$A1[d$next_trt==1]*d$A2[d$next_trt==1]*(1-d$cen2[d$next_trt==1])*d$wght[d$next_trt==1] )
  
  EY = EY_num / EY_denom
  
  return(EY)
}

run_iptw_boot = function(d, indices, sval, adj=1){
  d_boot = d[indices, ]
  
  res = numeric(length(sval))
  
  for(j in 1:length(sval) ){
    res[j] = run_iptw_s2(sval[j], d_boot, adj)  
  }
  
  return(res)
}

#------------------------------------------------------------------------------#
## Discrete-Time Methods  ----
#------------------------------------------------------------------------------#

run_iptw_discrete = function(d, adj=1){
  
  ## partition interval
  endpoints = seq(1, ceiling(max(d$Y)), 1) 
  start = c(0,endpoints)[ 1:( (length(endpoints) +1) - 1 )]
  end = c(0,endpoints)[ 2:(length(endpoints) + 1) ]
  n_int = length(end)
  
  d$W_int = findInterval(d$W1, vec = c(0, endpoints))
  d$W_int[!d$next_trt] = Inf
  
  
  ds2 = contToDisc(d, timeColumn = 'Y', intervalLimits = endpoints )
  ds3 = discSurv::dataLong(ds2, 'timeDisc', eventColumn = 'delta')
  
  ds3$start = start[ds3$timeInt]
  ds3$end = end[ds3$timeInt]
  
  ## interval of second treatment start
  ds3$W1_int = ifelse( ds3$start < ds3$W1 & ds3$W1 < ds3$end & ds3$next_trt, ds3$timeInt, NA)
  ## has second treatment decision occurred by interval k?
  ds3$W1k = ifelse(ds3$W1 < ds3$end & ds3$next_trt, 1, 0) 
  ## which second treatment has been administered by interval k?
  ds3$A2k = ifelse(ds3$W1 < ds3$end & ds3$next_trt, ds3$A2, NA) 
  
  
  ###--- censoring model ---###
  ds2$delta_cen = 1 - ds2$delta
  dsc = discSurv::dataLong(ds2, 'timeDisc', eventColumn = 'delta_cen')
  dsc$timeInt = as.factor(dsc$timeInt)
  
  dsc_pre = dsc[as.numeric(dsc$timeInt) < dsc$W_int,]
  
  haz_c_pre = glm(data=dsc_pre, formula = y ~ A1 + L1, family=binomial('logit'))
  dsc_pre$haz_c_pre = 1-predict(haz_c_pre, dsc_pre, 'response')

  dsc_post=dsc[as.numeric(dsc$timeInt) >= dsc$W_int,]
  
  haz_c_post = glm(data=dsc_post, formula = y ~ timeInt + A2 + L2 + I(W_int>15), family=binomial('logit')  )  
  dsc_post$haz_c_post = 1-predict(haz_c_post, dsc_post, 'response')
  
  ###---                 ---###
  
  ### --- Train propensity score model for treatment 1 --- ###
  ds3_bsl = ds3[ds3$timeInt==1,]
  mod = glm(data = ds3_bsl, formula =  A1 ~ L1, family=binomial('logit'))
  
  ### --- Estimate propensity score for treatment 1 --- ###
  ds3_bsl$pscore1 = predict(mod, newdata = ds3_bsl, type='response')
  
  ### --- Train propensity score model for treatment 2 --- ###
  ds3_bsl2 = ds3[!is.na(ds3$W1_int) & ds3$next_trt,]
  if(adj==1){
    mod2 = glm(data=ds3_bsl2, formula = A2 ~ L2 + I(W1_int>15), family=binomial('logit')) 
  }else{
    mod2 = glm(data=ds3_bsl2, formula = A2 ~ L2, family=binomial('logit')) 
  }
  ### --- Estimate propensity score for treatment 2 --- ###
  ds3_bsl2$pscore2 = predict(mod2, newdata = ds3_bsl2, type='response')
  
  # calculuate time-varying treamtent weights and subset to just 
  #   person-time consistent with (1,1) treatment strategy
  dsc_pre$timeInt = as.numeric(dsc_pre$timeInt)
  dsc_post$timeInt = as.numeric(dsc_post$timeInt)
  
  dmsm = ds3 %>%
    left_join( select(ds3_bsl, obj, timeInt,  pscore1), by=c('obj','timeInt') ) %>%
    left_join( select(ds3_bsl2,obj, timeInt, pscore2), by=c('obj','timeInt') ) %>%
    left_join( select(dsc_pre,obj, timeInt, haz_c_pre), by=c('obj','timeInt') ) %>%
    left_join( select(dsc_post,obj, timeInt, haz_c_post), by=c('obj','timeInt') ) %>%
    filter( (next_trt & A1==1 & A2==1) | (!next_trt & A1==1) | (next_trt & A1==1 & A2==0) ) %>%
    ## exclude person-time that is no longer consistent with strategy
    mutate(exclude = ifelse(next_trt & !is.na(A2k) & A2k==0,1,0 )) 
  
  dmsm  = dmsm %>% filter(exclude==0)
  
  dmsm$pscore1[is.na(dmsm$pscore1)] = 1
  dmsm$pscore2[is.na(dmsm$pscore2)] = 1
  dmsm$haz_c_pre[is.na(dmsm$haz_c_pre)] = 1
  dmsm$haz_c_post[is.na(dmsm$haz_c_post)] = 1
  
  dmsm = dmsm %>%
    group_by(obj) %>%
    mutate(joint_pscore = pscore1*pscore2, 
           tv_joint_pscore = cumprod(joint_pscore), 
           joint_cscore = haz_c_pre*haz_c_post,
           tv_joint_cscore = cumprod(joint_cscore), 
           combined_score = tv_joint_pscore*tv_joint_cscore)  %>%
    mutate(weight = 1/combined_score)
  
  dmsm$timeInt = factor(dmsm$timeInt, levels=1:n_int, labels=1:n_int)
  
  ## estimate marginal structural hazard model parameters via weighted 
  ## logistic regression 
  mod_y = glm(data = dmsm, formula =  y ~ timeInt , 
              family=binomial('logit'), weights = weight)
  
  ## back out the survival probability from the hazards
  newdata = data.frame( timeInt = as.factor(1:20) )
  newdata$haz = (1-predict(mod_y, newdata = newdata, type='response'))
  newdata$surv_prob = cumprod(newdata$haz)
  
  return(newdata)
}

run_iptw_discrete_boot = function(d, indices , sval , adj=1){
  d_boot = d[ indices , ]
  res = run_iptw_discrete_s2( d_boot , adj )
  return(res$surv_prob[sval])
}

