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
                          k = 5,
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
    
    
    mod3 <- Djordjilovic(X = X, Y = Y, M = M, k = K) 
    mod4 <- Tobi(X = X, Y = Y, M = M, k = K)
    mod5 <- hima2(X = X, Y = Y, M = M, k = K, conf = NULL)
    
    # max2 
    m3 <- max2(pval1 = mod3$calibrated.pvalue[, 1], pval2 = mod3$calibrated.pvalue[, 2])
    
    # tobi
    
    m4 <- fdrtool::fdrtool(as.vector(mod4$calibrated.pvalue), statistic = "pvalue", plot = F)$qval
    
    # sreenmin et sbmh
    ss3 <- combi_pv(p1 = mod3$calibrated.pvalue[, 1], p2 = mod3$calibrated.pvalue[, 2], 
                    FDR = 0.1)$pv.fdr[, 1:2]
    
    # tobi 
    tmp <- as.vector(mod4$pValue)
    tmp2 <- stats::qnorm(tmp)
    score <- calibrate_by_gif(as.matrix(tmp2))
    qv <- fdrtool::fdrtool(as.vector(score$calibrated.pvalue), statistic = "pvalue", plot = F)
    
    # hdmt 
    
    input_pvalues <- mod3$calibrated.pvalue
    nullprop <- HDMT::null_estimation(input_pvalues, lambda = 0.5)
    
    pnull1 <- HDMT::adjust_quantile(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                                    nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)
    
    pmax <- apply(input_pvalues, 1, max)
    
    hdmt <- HDMT::fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                          nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)
    
    
    qval <- data.frame(VanderWeele_max2 = m3$qval, VanderWeele_SM = ss3$screeM, 
                       VanderWeele_SBMH = ss3$SBMH, VanderWeele_HDMT = hdmt,
                       Tobi = qv$qval, hima = mod5)
    
    F1 <- matrix(NA, nrow = 1, ncol = ncol(qval))
    AUC <- matrix(NA, nrow = 1, ncol = ncol(qval))
    PRE <- matrix(NA, nrow = 1, ncol = ncol(qval))
    REC <- matrix(NA, nrow = 1, ncol = ncol(qval))
    
    for (j in 1:ncol(qval)) {
      
      rp <- useFunc::rank.pwer(pval = qval[, j], known.mediator = causal, 
                               toplist = length(causal), 
                               ral = 0.05 * ncol(M))
      
      F1[j] <- rp$f1_score 
      AUC[j] <- rp$auc.norm
      PRE[j] <- rp$precision 
      REC[j] <- rp$recall
      
      
    }
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
  
  a <- c("max2", "SM", "SBMH", "HDMT", "Tobi", "hima")
  
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

saveRDS(res, file = paste0("data_res/SIMU-mediation_", Sys.Date(), ".RDS"))
