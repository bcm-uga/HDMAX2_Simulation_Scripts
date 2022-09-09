library(useFunc)
library(lfmm)
library(corrplot)
library(tidyverse)
library(reshape2)
library(foreach)

source("../function/function.R")

'%!in%' <- function(x,y){
  !('%in%'(x,y))
}


# estimate CT prop 
prop.CT <- readr::read_csv("../data/omega_selprobes.csv")
prop.CT <- prop.CT[, -c(1:2)]
prop.CT <- colMeans(prop.CT)


simu.parallel <- function(n = 400,
                          p = 10000,
                          k = 2,
                          nc = 10,
                          propVX = 0.1,
                          propVY = 0.1,
                          rho = rep(seq(0, 0.3, 0.1), 3),
                          sigma = 1,
                          meanA = c(0.25, 0.5),
                          meanB = c(0.25, 0.5),
                          sdU = 1,
                          sdV = 1,
                          nb.core = 1,
                          prop.CT = prop.CT,
                          log = "test/log.txt") {
  
  
  
  st <- expand.grid(n, p, k, nc, propVX, propVY, rho, sigma, meanA, meanB, sdU, sdV)
  colnames(st) <- c("n", "p", "k", "nc", "propVX", "propVY", 
                    "rho", "sigma", "meanA", "meanB", "sdU", "sdV")
  
  
  nb.sim <- nrow(st)
  
  print(paste0("number of simulation : ", nb.sim))
  
  writeLines(c(""), log)
  
  SIMU <- function(i) {
    
    max2 <- function (pval1, pval2, diagnostic.plot = F) 
    {
      pval <- apply(cbind(pval1, pval2), 1, max)^2
      eta0 <- fdrtool::pval.estimate.eta0(pval, diagnostic.plot = diagnostic.plot)
      qval <- fdrtool::fdrtool(pval, statistic = "pvalue", plot = F, 
                               verbose = F, ...)$qval
      return(list(pval = pval, eta0 = eta0, qval = qval))
    }
    
    source("../function/function.R") # ATTENTION REPERTOIRE
    source("../function/refactor.R") # ATTENTION REPERTOIRE
    # be careful --> ici cY = nc 
    
    # BE CAREFUL this sdu is correct for the chromosome 1 of EDEN
    
    propCX = 2 * (st$nc[i] / st$p[i])
    propCY = 2 * (st$nc[i] / st$p[i])
    
    simu <- r_mediation2(n = st$n[i], p = st$p[i], K = st$k[i], 
                         prop.causal.x = propCX,  prop.causal.y = propCY, 
                         prop.variance.x = st$propVX[i], prop.variance.y = st$propVY[i], 
                         rho = st$rho[i],
                         sigma = st$sigma[i], 
                         mean.A = st$meanA[i], 
                         mean.B = st$meanB[i], 
                         # sd.U = st$sdU[i],
                         sd.U = c(3, 2, 1),
                         sd.V = st$sdV[i], 
                         prop.causal.ylx = 1, 
                         sd.A = 0.02, sd.B = 0.02, 
                         prop.CT = prop.CT,
                         K.ct = length(prop.CT)) # From paulina data
    
    X <- simu$X
    Y <- simu$Y
    M <- simu$M
    causal <- simu$mediators
    K <- simu$K + simu$K.ct
    
    # lfmm
    vande <- Djordjilovic(X = X, Y = Y, M = M, k = K)
    tmp <- max2(pval1 = vande$pValue[, 1], pval2 = vande$pValue[, 2])$qval
    vande <- max2(pval1 = vande$calibrated.pvalue[, 1], pval2 = vande$calibrated.pvalue[, 2])$qval
    
    lfm <- data.frame(cal = vande, nocal = tmp)
    
    # reffreeewas
    ct <- RefFreeEWAS::RefFreeCellMix(t(M), K = K)$Omega
    conf <- NULL
    
    #### VanderWeele
    
    # first regression
    
    # pValeur
    pv1 <- matrix(nrow = ncol(M), ncol = 1)
    
    # score
    sc1 <- matrix(nrow = ncol(M), ncol = 1)
    
    
    for (h in 1:ncol(M)) {
      
      dat <- data.frame(mi = M[, h], X = X, cbind(ct, conf))
      dat <- summary(lm(mi ~ X + ., data = dat))$coeff[2, 3:4]
      
      sc1[h] <- dat[1]
      pv1[h] <- dat[2]
      
    }
    
    # Calibrate
    
    dat <- calibrate_by_gif(as.matrix(sc1))
    
    sc1.cal <- dat$calibrated.score2[, 1]
    pv1.cal <- dat$calibrated.pvalue[, 1]  
    
    
    # pValeur
    pv2 <- matrix(nrow = ncol(M), ncol = 1)
    
    # score
    sc2 <- matrix(nrow = ncol(M), ncol = 1)
    
    
    # Second regression
    for (h in 1:ncol(M)) {
      
      dat <- data.frame(Y = Y, mi = M[, h], X = X, cbind(ct, conf))
      dat <- summary(lm(Y ~ mi + X + ., data = dat))$coeff[2, 3:4]
      
      sc2[h] <- dat[1]
      pv2[h] <- dat[2]
      
    }
    
    # Calibrate
    
    dat <- calibrate_by_gif(as.matrix(sc2))
    
    sc2.cal <- dat$calibrated.score2[, 1]
    pv2.cal <- dat$calibrated.pvalue[, 1]
    
    pv.cal <- max2(pval1 = pv1.cal, pval2 = pv2.cal)$qval
    pv <- max2(pval1 = pv1, pval2 = pv2)$qval
    
    rfe <- data.frame(pv.cal, pv)
    
    
    # refactor
    ct <- refactor2(data = t(M), k = K)$refactor_components
    conf <- NULL
    
    #### VanderWeele
    
    # first regression
    
    # pValeur
    pv1 <- matrix(nrow = ncol(M), ncol = 1)
    
    # score
    sc1 <- matrix(nrow = ncol(M), ncol = 1)
    
    
    for (h in 1:ncol(M)) {
      
      dat <- data.frame(mi = M[, h], X = X, cbind(ct, conf))
      dat <- summary(lm(mi ~ X + ., data = dat))$coeff[2, 3:4]
      
      sc1[h] <- dat[1]
      pv1[h] <- dat[2]
      
    }
    
    # Calibrate
    
    dat <- calibrate_by_gif(as.matrix(sc1))
    
    sc1.cal <- dat$calibrated.score2[, 1]
    pv1.cal <- dat$calibrated.pvalue[, 1]  
    
    
    # pValeur
    pv2 <- matrix(nrow = ncol(M), ncol = 1)
    
    # score
    sc2 <- matrix(nrow = ncol(M), ncol = 1)
    
    
    # Second regression
    for (h in 1:ncol(M)) {
      
      dat <- data.frame(Y = Y, mi = M[, h], X = X, cbind(ct, conf))
      dat <- summary(lm(Y ~ mi + X + ., data = dat))$coeff[2, 3:4]
      
      sc2[h] <- dat[1]
      pv2[h] <- dat[2]
      
    }
    
    # Calibrate
    
    dat <- calibrate_by_gif(as.matrix(sc2))
    
    sc2.cal <- dat$calibrated.score2[, 1]
    pv2.cal <- dat$calibrated.pvalue[, 1]
    
    pv.cal <- max2(pval1 = pv1.cal, pval2 = pv2.cal)$qval
    pv <- max2(pval1 = pv1, pval2 = pv2)$qval
    
    
    ref <- data.frame(pv.cal, pv)
    
    # stat
    
    res <- data.frame(lfmm = lfm$nocal, lfmm_cal = lfm$cal, 
                      RefFreeEWAS = rfe$pv, RefFreeEWAS_cal = rfe$pv.cal,
                      refactor = ref$pv, refactor_cal = ref$pv.cal)
    
    F1 <- matrix(NA, nrow = 1, ncol = 6)
    AUC <- matrix(NA, nrow = 1, ncol = 6)
    PRE <- matrix(NA, nrow = 1, ncol = 6)
    REC <- matrix(NA, nrow = 1, ncol = 6)
    
    for (j in 1:ncol(res)) {
      
      rp <- useFunc::rank.pwer(pval = res[, j], known.mediator = causal, 
                               toplist = length(causal), 
                               ral = 0.05 * st$p[i])
      
      F1[j] <- rp$f1_score 
      AUC[j] <- rp$auc.norm
      PRE[j] <- rp$precision 
      REC[j] <- rp$recall
      
    }
    
    print(F1)
    return(c(i, F1, AUC, PRE, REC))
  }
  
  cl <- parallel::makeCluster(nb.core)
  doParallel::registerDoParallel(cl)
  
  res <- foreach(i = 1:nb.sim, .combine = 'rbind', .packages = c("parallel"), 
                 .errorhandling = 'pass') %dopar% {
                   
                   sink(log, append = T)
                   cat(paste("Simulation number : ", i, "\n"))
                   sink()
                   
                   SIMU(i)
                 }
  
  parallel::stopCluster(cl)
  
  index <- res[,1]
  
  res <- res[, -1]
  
  a <- c("lfmm", "lfmm_cal", "RefFreeEWAS", "RefFreeEWAS_cal", "refactor", "refactor_cal")
  
  res[is.na(res)] <- 0
  
  F1 <- res[, 1:6]
  AUC <- res[, 7:12]
  PRE <- res[, 13:18]
  REC <- res[, 19:24]
  
  colnames(F1) <- a
  colnames(AUC) <- a
  colnames(PRE) <- a
  colnames(REC) <- a
  
  return(list(F1 = data.frame(st[, apply(st, 2, var) != 0], F1),
              AUC = data.frame(st[, apply(st, 2, var) != 0], AUC),
              PRE = data.frame(st[, apply(st, 2, var) != 0], PRE),
              REC = data.frame(st[, apply(st, 2, var) != 0], REC),
              para = st,
              index = index))
}


res <- simu.parallel(n = 500,
                     p = 38000,
                     k = 3,
                     nc = rep(c(8, 16, 32), 300),
                     propVX = 0.1,
                     propVY = 0.1,
                     rho = 0.1,
                     sigma = 1,
                     meanA = c(0.2, 0.4),
                     meanB = c(0.2, 0.4),
                     sdU = 1,
                     sdV = 1,
                     nb.core = 14,
                     log = "test/log.txt", prop.CT = prop.CT)

saveRDS(res, file = paste0("data_res/SIMU-latent_", Sys.Date(), ".RDS"))
