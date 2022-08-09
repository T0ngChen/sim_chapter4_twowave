##################
##              ##      
## rare disease ##
##              ##
##################
impute.id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
name = paste0("sim", impute.id, ".rds")




if (impute.id == 1){
  beta = c(-3.553709, 0, 1)
  para.val = c(-3, 3, 1, 1)
} else if (impute.id == 2){
  beta = c(-3.379317, -0.5, 1)
  para.val = c(-3, 3, 1, 1)
} else if (impute.id == 3){
  beta = c(-3.251742, -1, 1)
  para.val = c(-3, 3, 1, 1)
} else if (impute.id == 4){
  beta = c(-3.553709, 0, 1)
  para.val = c(-3, 4, 1, 1)
} else if (impute.id == 5){
  beta = c(-3.379317, -0.5,1)
  para.val = c(-3, 4, 1, 1)
} else if (impute.id == 6){
  beta = c(-3.251742, -1, 1)
  para.val = c(-3, 4, 1, 1)
} else if (impute.id == 7){
  beta = c(-3.553709, 0, 1)
  para.val = c(-3, 5, 1, 1)
} else if (impute.id == 8){
  beta = c(-3.379317, -0.5, 1)
  para.val = c(-3, 5, 1, 1)
} else if (impute.id == 9){
  beta = c(-3.251742, -1, 1)
  para.val = c(-3, 5, 1, 1)
}


library(survey)
options(survey.lonely.psu="remove")


## influence function for logistic regression
inf.fun <- function(fit) {
  dm <- model.matrix(fit)
  Ihat <- (t(dm) %*% (dm * fit$fitted.values * (1 - fit$fitted.values))) / nrow(dm)
  ## influence function
  infl <- (dm * resid(fit, type = "response")) %*% solve(Ihat)
  infl
}

#### Influence function from a weighted model
inf.fun.w <- function(fit) {
  dm <- model.matrix(fit)
  w <- fit$prior.weights
  Ihat <- (t(dm) %*% (dm * w * fit$fitted.values * (1 - fit$fitted.values)))
  ## influence function
  infl <- (dm * w * resid(fit, type = "response")) %*% solve(Ihat)
  infl
}


## Exact Neyman allocation (Wright, 2017) https://doi.org/10.1016/j.spl.2017.04.026
integer.neyman.w1 = function(n.strata, NS, sample.size, upper){
  nc = max(upper+1)
  s = 1:nc
  arr = matrix(rep(NS, each = nc)/sqrt(s*(s+1)), nrow = n.strata, byrow = T)
  arr.list = as.list(as.data.frame(t(arr)))
  for(i in 1:length(arr.list)){
    arr.list[[i]][upper[i]:nc] = 0
  }
  arr = do.call(cbind, arr.list)[-1,]
  rk <- order(arr, na.last=TRUE, decreasing=TRUE)[1:(sample.size - 2 * n.strata)]
  re.om.zero = table(rk%/%(nc-1) + 1)
  re.zero = rep(2, n.strata)
  re.zero[1:n.strata %in% names(re.om.zero)] = 2 + re.om.zero
  re.zero
}

integer.neyman.w2 = function(n.strata, NS, sample.size, upper){
  nc = max(upper + 1)
  s = 1:nc
  arr = matrix(rep(NS, each = nc)/sqrt(s*(s+1)), nrow = n.strata, byrow = T)
  arr.list = as.list(as.data.frame(t(arr)))
  for(i in 1:length(arr.list)){
    arr.list[[i]][upper[i]:nc] = 0
  }
  arr = do.call(cbind, arr.list)
  rk <- order(arr, na.last=TRUE, decreasing=TRUE)[1:(sample.size - n.strata)]
  re.om.zero = table(rk%/%nc + 1) + 1
  re.zero = rep(0, n.strata)
  re.zero[1:n.strata %in% names(re.om.zero)] = re.om.zero
  re.zero[which.max(upper - re.zero)] =  re.zero[which.max(upper - re.zero)] + length(which(re.zero == 0))
  re.zero
}

## perform raking and IPW

impute.raking = function(data){
  twophase.w2 <- twophase(id = list(~1, ~1), strata = list(NULL, ~stra),
                          subset = ~insample, data = data)
  svyest_ipw = svyglm(Y ~ X + Z, family = quasibinomial, design = twophase.w2)
  impmodel <- svyglm(X ~ X_tilde + Z, family = quasibinomial, design = twophase.w2)
  data$imputex <- as.vector(predict(impmodel, newdata=data, type="response", se.fit=FALSE))
  phase1model_imp <- glm(Y ~ imputex + Z, family = binomial, data = data)
  inffun_imp <- inf.fun(phase1model_imp)
  colnames(inffun_imp)<-paste0("if",1:ncol(inffun_imp))
  twophase.w2.imp <- twophase(id=list(~1,~1),  strata = list(NULL, ~stra),
                              subset = ~insample, data = cbind(data, inffun_imp), method="simple")
  calformula <- make.formula(colnames(inffun_imp))
  cal_twophase_imp <- calibrate(twophase.w2.imp, calformula, phase=2, calfun = "raking")
  svyest_rak<-svyglm(Y ~ X + Z, family = quasibinomial, design=cal_twophase_imp)
  list(svyest_ipw, svyest_rak)
}






pps.prob = function(n1 = n1, w1.n = n){
  if(any(n1/sum(n1) * w1.n < 3)){
    prob = n1/sum(n1)
    prob[prob<0.05] = 0.05
    prob[which.max(prob)] = prob[which.max(prob)] - (sum(prob) - 1)
  } else {
    prob = n1/sum(n1)
  }
  prob
}



one.sim<-function(beta, N = 15000, n = 400, para.val){
  ## generate data
  df <- data.frame(Z = rbinom(N, 1, .5))
  df$X <- rbinom(N, size = 1, prob = 0.4)
  modelterm = beta[1] + beta[2]*df$X + beta[3]*df$Z
  p_y <- exp(modelterm) / (1 + exp(modelterm))
  df$Y <- rbinom(N, size = 1, p_y)
  
  dm = cbind(rep(1, nrow(df)), df$X, df$Y, df$Z)
  fitted = exp(rowSums(dm %*% diag(para.val))) / (1+exp(rowSums(dm %*% diag(para.val))))
  
  
  df$X_tilde = rbinom(N,1,fitted)
  
  df$id<-1:nrow(df)
  
  # stratification on X_tilde
  df$stra = as.numeric(interaction(df$X_tilde, df$Y))
  n1 = xtabs(~stra, df)
  
  strata = list()
  for (index in 1:length(n1)){
    strata[[index]] = which(df$stra == index)  
  }
  
  
  
  
  ###########################
  ##                       ##
  ## Optimal design IF-IPW ##
  ##                       ##
  ###########################
  infl = inf.fun(glm(Y~X+Z, data=df, family = binomial))[, 2]
  sd.stra1 = numeric()
  for(i in 1:length(strata)){
    sd.stra1 = c(sd.stra1, sd(infl[strata[[i]]]))
  }
  ney = integer.neyman.w1(n.strata = length(strata), NS = n1 * sd.stra1, sample.size = n, upper = n1)
  
  
  ## sampling and calibration
  df0 = df
  s.ney = list()
  for(i in 1:length(strata)){
    s.ney[[i]] = sample(strata[[i]], ney[i])
  }
  df0$insample = (1:nrow(df0)) %in% unlist(s.ney)
  IFIPW = impute.raking(df0)
  ipw.IFIPW = IFIPW[[1]]
  rak.IFIPW = IFIPW[[2]]
  
  
  
  ##############
  ##          ##
  ##   BSS    ##
  ##          ##
  ##############
  df1 = df
  s.bss = list()
  for(i in 1:length(strata)){
    s.bss[[i]] = sample(strata[[i]], round(n/length(strata)))
  }
  
  df1$insample = (1:nrow(df1)) %in% unlist(s.bss)
  BSS = impute.raking(df1)
  ipw.BSS = BSS[[1]]
  rak.BSS = BSS[[2]]
  
  ##############
  ##          ##
  ##   PSS    ##
  ##          ##
  ##############
  df2 = df
  s.pss = list()
  pps.num = round(n * pps.prob(n1 = n1, w1.n = n))
  for(i in 1:length(strata)){
    s.pss[[i]] = sample(strata[[i]], pps.num[i])
  }
  
  df2$insample = (1:nrow(df2)) %in% unlist(s.pss)
  PSS = impute.raking(df2)
  ipw.PSS = PSS[[1]]
  rak.PSS = PSS[[2]]
  
  
  ##############
  ##          ##
  ##  wave-1  ##
  ##          ##
  ##############
  df3 = df
  wave1 = list()
  for(i in 1:length(strata)){
    wave1[[i]] = sample(strata[[i]], n/2/length(strata))
  }
  df3$w1.insample = (1:nrow(df3)) %in% unlist(wave1)
  
  ################
  ##            ##
  ##  weighted  ##
  ##            ##  
  ################
  df4 = df3
  design <- twophase(id=list(~1,~1), strata=list(NULL,~stra),
                     subset=~w1.insample,
                     data=df4, method="simple")
  fit.wt = svyglm(Y~X+Z, family="quasibinomial", design=design)
  infl.w1.wt<-inf.fun.w(fit.wt)[,2]
  df4$infl = NA
  df4[df4$w1.insample,]$infl = infl.w1.wt
  
  
  sd.w1.wt = numeric()
  for(i in 1:length(strata)){
    sd.w1.wt = c(sd.w1.wt, sd(df4$infl[strata[[i]]][!is.na(df4$infl[strata[[i]]])]))
  }
  wt.w2.size = integer.neyman.w2(n.strata = length(strata), NS = n1 * sd.w1.wt, sample.size = 200, upper = n1)
  # sample wave 2
  wt.wave2 = list()
  for (i in 1:length(strata)) {
    if(length(strata[[i]][-match(wave1[[i]], strata[[i]])]) != 0){
      wt.wave2[[i]] <- sample(strata[[i]][-match(wave1[[i]], strata[[i]])], wt.w2.size[i])
    } else if(length(strata[[i]][-match(wave1[[i]], strata[[i]])]) == 0){
      wt.wave2[[i]] <- NULL
    }
  }
  wt.sample = mapply(c, wave1, wt.wave2, SIMPLIFY=FALSE)
  df4$insample <- (1:nrow(df4)) %in% unlist(wt.sample)
  TWO.wt = impute.raking(df4)
  ipw.TWO.wt = TWO.wt[[1]]
  rak.TWO.wt = TWO.wt[[2]]
  
  ##################
  ##              ##
  ##  imputation  ##
  ##              ##  
  ##################
  df5 = df3
  design1 <- twophase(id = list(~1, ~1), strata = list(NULL, ~stra),
                      subset = ~w1.insample, data = df5)
  imp.model <- svyglm(X ~ X_tilde + Z, family = quasibinomial, design = design1)
  df5$imputex <- as.vector(predict(imp.model, newdata=df5, type="response", se.fit=FALSE))
  fit.imp <- glm(Y ~ imputex + Z, family = binomial, data = df5)
  df5$infl <- inf.fun(fit.imp)[,2]
  
  sd.w1.imp = numeric()
  for(i in 1:length(strata)){
    sd.w1.imp = c(sd.w1.imp, sd(df5$infl[strata[[i]]]))
  }
  imp.w2.size = integer.neyman.w2(n.strata = length(strata), NS = n1 * sd.w1.imp, sample.size = 200, upper = n1)
  
  imp.wave2 = list()
  for (i in 1:length(strata)) {
    if(length(strata[[i]][-match(wave1[[i]], strata[[i]])]) != 0){
      imp.wave2[[i]] <- sample(strata[[i]][-match(wave1[[i]], strata[[i]])], imp.w2.size[i])
    } else if(length(strata[[i]][-match(wave1[[i]], strata[[i]])]) == 0){
      imp.wave2[[i]] <- NULL
    }
  }
  imp.sample = mapply(c, wave1, imp.wave2, SIMPLIFY=FALSE)
  df5$insample <- (1:nrow(df5)) %in% unlist(imp.sample)
  TWO.imp = impute.raking(df5)
  ipw.TWO.imp = TWO.imp[[1]]
  rak.TWO.imp = TWO.imp[[2]]
  
  
  coef.ipw = c((coef(ipw.IFIPW) - beta)[2], (coef(ipw.BSS) - beta)[2],
               (coef(ipw.PSS) - beta)[2], (coef(ipw.TWO.wt) - beta)[2],
               (coef(ipw.TWO.imp) - beta)[2])
  
  coef.rak = c((coef(rak.IFIPW) - beta)[2], (coef(rak.BSS) - beta)[2],
               (coef(rak.PSS) - beta)[2], (coef(rak.TWO.wt) - beta)[2],
               (coef(rak.TWO.imp) - beta)[2])
  
  sd.ipw = c(summary(ipw.IFIPW)$coefficients[2,2],
             summary(ipw.BSS)$coefficients[2,2],
             summary(ipw.PSS)$coefficients[2,2],
             summary(ipw.TWO.wt)$coefficients[2,2],
             summary(ipw.TWO.imp)$coefficients[2,2])
  
  sd.rak = c(summary(rak.IFIPW)$coefficients[2,2],
             summary(rak.BSS)$coefficients[2,2],
             summary(rak.PSS)$coefficients[2,2],
             summary(rak.TWO.wt)$coefficients[2,2],
             summary(rak.TWO.imp)$coefficients[2,2])
  
  
  
  sample.size = cbind(ney, unlist(lapply(wt.sample,length)),
                      unlist(lapply(imp.sample,length)))
  names(coef.ipw) = names(coef.rak) = names(sd.ipw) = names(sd.rak) = c("OPT", "BSS", "PSS", "TWO.wt", "TWO.imp")
  
  
  list(coef.ipw = coef.ipw, coef.rak = coef.rak,
       sd.rak = sd.rak, sd.ipw = sd.ipw,
       sample.size = sample.size,
       cor = cor(df$X,df$X_tilde),
       sen = table(df$X,df$X_tilde+2)[1,1]/sum(table(df$X,df$X_tilde+2)[1,]),
       spe = table(df$X,df$X_tilde+2)[2,2]/sum(table(df$X,df$X_tilde+2)[2,]))
}


sim.fun = function(){
  while(TRUE){
    df <- try(one.sim(beta = beta, para.val = para.val), silent=TRUE)
    if(!is(df, 'try-error')) break
  }
  df
}


set.seed(910111)
out = replicate(2000, sim.fun())

saveRDS(out, file = name)
