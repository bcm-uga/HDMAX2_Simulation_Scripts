library(foreach)
simu.parallel <- function(n = 400,
                          p = 10000,
                          k = 5,
                          nc = 10,
                          propVX = 0.1,
                          propVY = 0.1,
                          rho = 0.1,
                          meanA = c(0.25, 0.5),
                          meanB = c(0.25, 0.5),
                          sims = 10,
                          nb.core = 1,
                          log = "test/log.txt") {
  
  
  
  st <- expand.grid(n, p, k, nc, propVX, propVY, rho, meanA, meanB)
  colnames(st) <- c("n", "p", "k", "nc", "propVX", "propVY", "rho", "meanA", "meanB")
  
  # print(st)
  
  nb.sim <- nrow(st)
  
  print(paste0("number of simulation : ", nb.sim))
  
  writeLines(c(""), log)
  
  SIMU <- function(i) {
    
    source("../function/function.R")
    source("../function/refactor.R")
    # be careful --> ici cY = nc 
    
    # BE CAREFUL this sdu is correct for the chromosome 1 of EDEN
    
    propCX = 2 * (st$nc[i] / st$p[i])
    propCY = 2 * (st$nc[i] / st$p[i])
    
    simu <- useFunc::r_mediation(n = st$n[i], p = st$p[i], K = st$k[i], 
                                 prop.causal.x = propCX,  prop.causal.y = propCY, 
                                 prop.variance.x = st$propVX[i], prop.variance.y = st$propVY[i], 
                                 rho = st$rho[i],
                                 sigma = 1, 
                                 mean.A = st$meanA[i], mean.B = st$meanB[i], 
                                 sd.U = 1, 
                                 sd.V = 1, 
                                 prop.causal.ylx = 0.05, sd.A = 0.02, sd.B = 0.02)
    
    
    X <- simu$X
    Y <- simu$Y
    M <- simu$M
    
    causal <- simu$mediators
    colnames(M) <- paste0("M-", 1:ncol(M))
    
    mod <- Djordjilovic(X = X, Y = Y, M = M, k = 6)
    MAX2 <- highmed::max2(pval1 = mod$calibrated.pvalue[, 1], pval2 = mod$calibrated.pvalue[, 2])
    
    # refactor for cardenas
    ct <- refactor2(data = t(M), k = 6)$refactor_components
    
    # m <- M[, 1:100]
    
    pV.cardenas <- matrix(ncol = 2, nrow = ncol(M))
    pV.morales <- matrix(ncol = 1, nrow = ncol(M))
    sc.cardenas <- matrix(ncol = 2, nrow = ncol(M))
    sc.morales <- matrix(ncol = 1, nrow = ncol(M))
    
    
    for (h in 1:ncol(M)) {
      
      # EWAS M ~ X + ct + C
      mod <- summary(lm(M[, h] ~ X + ct))$coef
      pV.cardenas[h, 1] <- mod[2, 4]
      sc.cardenas[h, 1] <- mod[2, 3]
      
      # EWAS M ~ Y + ct + C
      mod <- summary(lm(M[, h] ~ Y + ct))$coef
      pV.cardenas[h, 2] <- mod[2, 4]
      sc.cardenas[h, 2] <- mod[2, 3]
      
      # EWAS M ~ X + C
      mod <- summary(lm(M[, h] ~ X))$coef
      pV.morales[h] <- mod[2, 4]
      sc.morales[h] <- mod[2, 3]
    }
    
    cal_cardenas <- calibrate_by_gif(sc.cardenas)
    cal_morales <- calibrate_by_gif(as.matrix(sc.morales))
    
    qv <- as.data.frame(matrix(1, ncol = 6, nrow = ncol(M)))
    colnames(qv) <- c("Cardenas", "Max2", "Morales", "Bonf-cardenas", "Bonf-Max2", "Bonf-morales")
    
    test <- sum(cal_cardenas$calibrated.pvalue[, 1] <= 0.05 / ncol(M) | cal_cardenas$calibrated.pvalue[, 2] <= 0.05)
    
    if (test > 0) {
      # cardenas
      m <- M[, cal_cardenas$calibrated.pvalue[, 1] <= 0.05 / ncol(M) | cal_cardenas$calibrated.pvalue[, 2] <= 0.05]
      cardenas <- highmed::univariate_mediation(qval = rep(0, ncol(m)), X = X, Y = Y, M = m, 
                                                covar = NULL, U = ct, FDR = 0.05, sims = sims)
      
      cardenas <- cardenas$ACME[, c(5, 4)]
      sign.car <- cardenas[cardenas$pval <= 0.05, ]
      bonf.car <- cardenas[cardenas$pval <= (0.05 / nrow(cardenas)), ]
      sign.car <- tidyr::separate(sign.car, CpG, c("M", "Posi"), "-")[, 2:3]
      bonf.car <- tidyr::separate(bonf.car, CpG, c("M", "Posi"), "-")[, 2:3]
      qv$Cardenas[as.numeric(sign.car$Posi)] <- 0
      qv$`Bonf-cardenas`[as.numeric(bonf.car$Posi)] <- 0
    }
    
    # morales
    # FDR correction
    qv.morales <- fdrtool::fdrtool(as.vector(cal_morales$calibrated.pvalue), 
                                   statistic = "pvalue", verbose = F, plot = F)
    test <- sum(qv.morales$qval <= 0.05)
    
    if (test > 0) {
      m <- M[, (qv.morales$qval <= 0.05)] 
      morales <- useFunc::sob.parallel(X = X, Y = Y, M = m, conf = NULL, nb.core = 1)
      
      morales <- cbind(CpG = colnames(M)[(qv.morales$qval <= 0.05)], morales)
      # data management
      morales <- morales[, 1:2]
      sign.mor <- morales[morales$pv.sobel <= 0.05, ]
      bonf.mor <- morales[morales$pv.sobel <= (0.05 / nrow(morales)), ]
      sign.mor <- tidyr::separate(sign.mor, CpG, c("M", "Posi"), "-")[, 2:3]
      bonf.mor <- tidyr::separate(bonf.mor, CpG, c("M", "Posi"), "-")[, 2:3]
      qv$Morales[as.numeric(sign.mor$Posi)] <- 0
      qv$`Bonf-morales`[as.numeric(bonf.mor$Posi)] <- 0
    }
    
    test <- sum(MAX2$pval <= (0.05 / ncol(M)))
    # test <- sum(MAX2$qval <= 0.05)
    
    if (test > 0) {
      # max2
      # m <- M[, MAX2$qval <= 0.05]
      m <- M[, MAX2$pval <= (0.05 / ncol(M))]
      max2 <- highmed::univariate_mediation(qval = rep(0, ncol(m)), X = X, Y = Y, M = m, 
                                            covar = NULL, U = ct, FDR = 0.05, sims = 10)
      max2 <- max2$ACME[, c(5, 4)]
      # mediation signi ?
      sign.max <- max2[max2$pval <= 0.05, ]
      # mediation signi Bonferroni
      bonf.max <- max2[max2$pval <= (0.05 / nrow(max2)), ]
      # data manage
      sign.max <- tidyr::separate(sign.max, CpG, c("M", "Posi"), "-")[, 2:3]
      bonf.max <- tidyr::separate(bonf.max, CpG, c("M", "Posi"), "-")[, 2:3]
      qv$Max2[as.numeric(sign.max$Posi)] <- 0
      qv$`Bonf-Max2`[as.numeric(bonf.max$Posi)] <- 0
    }
    
    # f1 score
    
    F1 <- matrix(NA, nrow = 1, ncol = 6)
    AUC <- matrix(NA, nrow = 1, ncol = 6)
    PRE <- matrix(NA, nrow = 1, ncol = 6)
    REC <- matrix(NA, nrow = 1, ncol = 6)
    
    for (j in 1:ncol(qv)) {
      
      rp <- useFunc::rank.pwer(pval = qv[, j], known.mediator = causal, toplist = length(causal), 
                               ral = 0.05)
      
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
  
  # print(index)
  
  res <- res[, -1]
  
  a <- c("Cardenas", "Max2", "Morales", "Bonf-cardenas", "Bonf-Max2", "Bonf-morales")
  
  res[is.na(res)] <- 0
  
  F1 <- res[, 1:6]
  AUC <- res[, 7:12]
  PRE <- res[, 13:18]
  REC <- res[, 19:24]
  # 
  colnames(F1) <- a
  colnames(AUC) <- a
  colnames(PRE) <- a
  colnames(REC) <- a
  # 
  return(list(F1 = data.frame(st, F1),
              AUC = data.frame(st, AUC),
              PRE = data.frame(st, PRE),
              REC = data.frame(st, REC),
              para = st,
              index = index))
  
  # return(res)
}

res <- simu.parallel(n = 500,
                     p = 38000,
                     k = 6,
                     nc = rep(c(8, 16, 32), 200),
                     propVX = 0.1,
                     propVY = 0.1,
                     rho = 0.1,
                     meanA = c(0.2, 0.4, 0.8, 1),
                     meanB = c(0.2, 0.4, 0.8, 1),
                     nb.core = 10,
                     log = "test/log.txt", sims = 100)

saveRDS(res, "data_res/simu_cardenas_morales_raw.RDS")

# F1
f1 <- res$F1
meth <- f1[, 10:15]
meth <- as.data.frame(meth)

F1 <- matrix(NA, ncol = 6, nrow = nrow(f1))

for (j in 1:ncol(meth)) {
  tmp <- meth[, j]
  
  temp <- matrix(NA, ncol = 1, nrow = length(tmp))
  
  for (i in 1:length(tmp)) {
    if (is.numeric(tmp[[i]])) {
      temp[i, ] <- tmp[[i]]
    }
  }
  F1[, j] <- temp
}

colnames(F1) <- colnames(f1)[10:15]

f1 <- data.frame(f1[, 1:9], F1)

# pre
pre <- res$PRE
meth <- pre[, 10:15]
meth <- as.data.frame(meth)

PRE <- matrix(NA, ncol = 6, nrow = nrow(f1))

for (j in 1:ncol(meth)) {
  tmp <- meth[, j]
  
  temp <- matrix(NA, ncol = 1, nrow = length(tmp))
  
  for (i in 1:length(tmp)) {
    if (is.numeric(tmp[[i]])) {
      temp[i, ] <- tmp[[i]]
    }
  }
  PRE[, j] <- temp
}

colnames(PRE) <- colnames(pre)[10:15]

pre <- data.frame(pre[, 1:9], PRE)

# F1
rec <- res$REC
meth <- rec[, 10:15]
meth <- as.data.frame(meth)

REC <- matrix(NA, ncol = 6, nrow = nrow(f1))

for (j in 1:ncol(meth)) {
  tmp <- meth[, j]
  
  temp <- matrix(NA, ncol = 1, nrow = length(tmp))
  
  for (i in 1:length(tmp)) {
    if (is.numeric(tmp[[i]])) {
      temp[i, ] <- tmp[[i]]
    }
  }
  REC[, j] <- temp
}

colnames(REC) <- colnames(rec)[10:15]

rec <- data.frame(rec[, 1:9], REC)

res <- list(f1 = f1,
            pre = pre, 
            rec = rec)

# saveRDS(res, "data_res/.............RDS")
