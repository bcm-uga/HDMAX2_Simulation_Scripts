compute_gif <- function(score) {
  score2 <- score ^ 2
  apply(score2, 2, stats::median, na.rm = TRUE) / stats::qchisq(0.5, df = 1)
}

compute_pvalue_from_zscore2 <- function(score2, df = 1) {
  apply(score2, 1:2, function(z) stats::pchisq(z, lower.tail = FALSE, df = df))
}

calibrate_by_gif <- function(score) {
  gif <- compute_gif(score)
  calibrated.score2 <- sweep(score ^ 2, 2, gif, FUN = "/")
  calibrated.pvalue <- compute_pvalue_from_zscore2(calibrated.score2, df = 1)
  
  return(list(calibrated.score2 = calibrated.score2,
              calibrated.pvalue = calibrated.pvalue,
              gif = gif))
}


#### 
VanderWeele <- function(X, Y, M, k, conf = NULL) {
  
  # First regression with lfmm
  dat <- lfmm::lfmm_ridge(Y = M, X = cbind(X, conf), K = k)
  dat <- lfmm::lfmm_test(Y = M, X = cbind(X, conf), lfmm = dat)
  
  pv1 <- dat$pvalue[, 1]
  sc1 <- dat$score[, 1]
  
  sc1.cal <- dat$calibrated.score2[, 1]
  pv1.cal <- dat$calibrated.pvalue[, 1]
  
  gif1 <- dat$gif
  
  
    
  # Estimate latent factors for second regression
  U <- lfmm::lfmm_ridge(Y = M, X = cbind(X, Y, conf), K = k)$U
  
  # pValeur
  pv2 <- matrix(nrow = ncol(M), ncol = 1)
  
  # score
  sc2 <- matrix(nrow = ncol(M), ncol = 1)
  
  
  # Second regression
  for (i in 1:ncol(M)) {
    
    dat <- data.frame(Y = Y, mi = M[, i], X = X, cbind(U, conf))
    dat <- summary(lm(Y ~ mi + X + ., data = dat))$coeff[2, 3:4]
    
    sc2[i] <- dat[1]
    pv2[i] <- dat[2]
    
  }
  
  # Calibrate
  
  dat <- calibrate_by_gif(as.matrix(sc2))
  
  sc2.cal <- dat$calibrated.score2[, 1]
  pv2.cal <- dat$calibrated.pvalue[, 1]
  
  gif2 <- dat$gif
  
  return(list(score = cbind(sc1, sc2),
              pValue = cbind(pv1, pv2),
              calibrated.score2 = cbind(sc1.cal, sc2.cal),
              calibrated.pvalue = cbind(pv1.cal, pv2.cal),
              gif = c(gif1, gif2)))
  
}

EWAS_multivariate <- function(X, Y, M, k, conf = NULL) {
  
  dat <- lfmm::lfmm_ridge(Y = M, X = cbind(X, Y, conf), K = k)
  dat <- lfmm::lfmm_test(Y = M, X = cbind(X, Y, conf), lfmm = dat)
  
  return(list(score = dat$score[, 1:2],
              pValue = dat$pvalue[, 1:2],
              calibrated.score2 = dat$calibrated.score2[, 1:2],
              calibrated.pvalue = dat$calibrated.pvalue[, 1:2],
              gif = dat$gif))
}

Djordjilovic <- function(X, Y, M, k, conf = NULL) {
  
  res <- list()
  
  # First regression 
  dat <- lfmm::lfmm_ridge(Y = M, X = cbind(X, conf), K = k)
  res[[1]] <- dat 
  U1 <- dat$U
  dat <- lfmm::lfmm_test(Y = M, X = cbind(X, conf), lfmm = dat)
  
  pv1 <- dat$pvalue[, 1]
  sc1 <- dat$score[, 1]
  
  sc1.cal <- dat$calibrated.score2[, 1]
  pv1.cal <- dat$calibrated.pvalue[, 1]
  
  gif1 <- dat$gif[1]
  
  # Second regression
  dat <- lfmm::lfmm_ridge(Y = M, X = cbind(X, Y, conf), K = k)
  res[[2]] <- dat
  # ajout
  U2 <- dat$U
  dat <- lfmm::lfmm_test(Y = M, X = cbind(X, Y, conf), lfmm = dat)
  
  pv2 <- dat$pvalue[, 2]
  sc2 <- dat$score[, 2]
  
  sc2.cal <- dat$calibrated.score2[, 2]
  pv2.cal <- dat$calibrated.pvalue[, 2]
  
  gif2 <- dat$gif[2]
  
  names(res) <- c("mod1", "mod2")
  
  return(list(score = cbind(sc1, sc2),
              pValue = cbind(pv1, pv2),
              calibrated.score2 = cbind(sc1.cal, sc2.cal),
              calibrated.pvalue = cbind(pv1.cal, pv2.cal),
              gif = c(gif1, gif2),
              U1 = U1,
              U2 = U2,
              lfmm = res))
}

Tobi <- function(X, Y, M, k, conf = NULL) {
  
  # Estimate latent factors for second regression
  U <- lfmm::lfmm_ridge(Y = M, X = cbind(X, Y, conf), K = k)$U
  
  pValue <- matrix(nrow = ncol(M), ncol = 1)
  score <- matrix(nrow = ncol(M), ncol = 1)
  
  for (i in 1:ncol(M)) {
    
    dat <- data.frame(mi = M[, i], X = X, cbind(U, conf))
    p1 <- lm(mi ~ X + ., data = dat)
    
    dat <- data.frame(mi = M[, i], X = X, Y = Y, cbind(U, conf))
    p2 <- lm(mi ~ X + Y + ., data = dat)
    
    dat <- anova(p1, p2)[2, 5:6]
    
    score[i] <- dat[1, 1]
    pValue[i] <- dat[1, 2]
    
  }
  
  # Calibrate
  
  dat <- calibrate_by_gif(as.matrix(score))
  
  
  # ne pas utiliser la pValeur calibré car elle est calculé sur une F-stat
  return(list(score = score,
              pValue = pValue,
              calibrated.score2 = dat$calibrated.score2,
              calibrated.pvalue = dat$calibrated.pvalue,
              gif = dat$gif))
  
}

EWAS_two <- function(X, Y, M, k, conf = NULL) {
  
  # First regression 
  dat <- lfmm::lfmm_ridge(Y = M, X = cbind(X, conf), K = k)
  dat <- lfmm::lfmm_test(Y = M, X = cbind(X, conf), lfmm = dat)
  
  pv1 <- dat$pvalue[, 1]
  sc1 <- dat$score[, 1]
  
  sc1.cal <- dat$calibrated.score2[, 1]
  pv1.cal <- dat$calibrated.pvalue[, 1]
  
  gif1 <- dat$gif
  
  # Second regression
  dat <- lfmm::lfmm_ridge(Y = M, X = cbind(Y, conf), K = k)
  dat <- lfmm::lfmm_test(Y = M, X = cbind(Y, conf), lfmm = dat)
  
  pv2 <- dat$pvalue[, 1]
  sc2 <- dat$score[, 1]
  
  sc2.cal <- dat$calibrated.score2[, 1]
  pv2.cal <- dat$calibrated.pvalue[, 1]
  
  gif2 <- dat$gif
  
  return(list(score = cbind(sc1, sc2),
              pValue = cbind(pv1, pv2),
              calibrated.score2 = cbind(sc1.cal, sc2.cal),
              calibrated.pvalue = cbind(pv1.cal, pv2.cal),
              gif = c(gif1, gif2)))
  
}

ScreenMinBonf <- function(p1, p2, c){
  # SreenMin procedure
  #
  # Args:
  #  p1... an m-length vector of p-values (study1)
  #  p2... second m-length vector of p-values (study2)
  #  c... threshold for selection
  #
  # Returns:
  # an m-length vector of Bonferroni adjusted p-values
  #
  minp <- pmin(p1, p2)
  maxp <- pmax(p1, p2)
  
  m <- length(p1)
  
  S <- which(minp <=  c)
  SM <- which(maxp <= c / length(S))
  SM <- intersect(S, SM)
  adj.pval <- rep(1, m)
  names(adj.pval) <- 1:m
  adj.pval[SM] <- maxp[SM]*length(SM) 
  return(adj.pval)
}


ScreenMinBonf2 <- function(p1, p2, c){
  # SreenMin procedure
  #
  # Args:
  #  p1... an m-length vector of p-values (study1)
  #  p2... second m-length vector of p-values (study2)
  #  c... threshold for selection
  #
  # Returns:
  # an m-length vector of Bonferroni adjusted p-values
  #
  alpha <- c
  c <- c / length(p1)
  
  minp <- pmin(p1, p2)
  maxp <- pmax(p1, p2)
  
  m <- length(p1)
  
  S <- which(minp <=  c)
  SM <- which(maxp <= alpha / length(S))
  SM <- intersect(S, SM)
  adj.pval <- rep(1, m)
  names(adj.pval) <- 1:m
  adj.pval[SM] <- maxp[SM]*length(SM) 
  return(adj.pval)
}


combi_pv <- function(p1, p2, FDR = 0.05) {
  
  temp <- cbind(p1, p2)
  
  # max2
  pmax2 <- apply(temp, 1, max)^2
  
  # fisher
  fisher <- apply(temp, 1, function(x) metap::sumlog(x)$p)
  
  # SBMH
  pSBMH <- MultiMed::medTest.SBMH(pEM = temp[, 1], pMY = temp[,  2], MCP.type = "FDR")
  
  # ScreenMin
  screeM <- ScreenMinBonf2(p1 = temp[, 1], p2 = temp[,  2], c = FDR)
  
  return(list(pv = data.frame(pmax2 = pmax2, fisher = fisher),
              pv.fdr = data.frame(SBMH = pSBMH, 
                                  screeM = screeM, 
                                  pmax2 = p.adjust(pmax2, method = "fdr"),
                                  fisher = p.adjust(fisher, method = "fdr")),
              p1 = temp[, 1], p2 = temp[, 2]))
}


hima2 <- function(X, Y, M , k, conf) {
  U <- lfmm::lfmm_ridge(Y = M, X = cbind(X, Y, conf), K = k)$U
  colnames(M) <- 1:ncol(M)
  temp <- HIMA::hima(X, Y, M, COV.XM = cbind(U, conf), parallel = F, ncore = 1)
  hima <- rep(1, ncol(M))
  a <- row.names(temp)
  a <- gsub("`","",a)
  hima[as.numeric(a)] <- temp$BH.FDR 
  
  return(hima)
}

conditional_mediation_simu3 <- function(meth, nc = 5, eff_xm = 5, sig_x = 1, K = 5,
                                        eff_my = 5, eff_xy = 0.1, sd_min = 0.1, cY = 5, sig_y = 1) {
  
  eff_xm <- rep(eff_xm, nc)
  eff_my <- rep(eff_my, nc)
  
  n <- nrow(meth)
  p <- ncol(meth)
  
  # resampling
  sple <- sample(1:n, n, replace = F)
  meth_sple <- meth[sple, ]
  
  # causal mediator
  # only CpGs with an sd greater than: sd_min
  sd <- apply(meth_sple, 2, sd)
  list_mediator <- sample(which(sd >= sd_min), nc)
  mediator <- meth_sple[, list_mediator]
  
  # simulate X
  X <- useFunc::conditional_simu(mediator, eff_xm, sig_x)
  
  # pc <- prcomp(meth_sple)
  # factors <- pc$x[, 1:K]
  # sigma <- sqrt(sum(pc$sdev^2) - sum(pc$sdev[1:K]^2))
  
  if (cY >= 1) {
    sd <- apply(meth_sple[, -(list_mediator)], 2, sd)
    list_cY <- sample(which(sd >= sd_min), cY)
    cY <- meth_sple[, list_cY]
    
    mediator <- cbind(mediator, cY)
    eff_my <- rep(eff_my[1], ncol(mediator))
  }
  
  U <- rowSums(RSpectra::svds(meth_sple, K)$u)
  
  Y <- as.vector(eff_xy * X + eff_my %*% t(mediator) + U + rnorm(n, 0, sd = sig_y))
  
  return(list(X = X,
              M = meth_sple,
              Y = Y,
              mediator = list_mediator,
              dat = data.frame(mediator = list_mediator,
                               eff_xm = eff_xm,
                               eff_my = eff_my[1:nc],
                               eff_med = eff_xm * eff_my[1:nc])))
}

r_mediation2 <- function (n = 400, 
                          p = 10000,  
                          K = 2, 
                          K.ct = 5, 
                          freq = NULL, 
                          prop.causal.x = 0.01, 
                          prop.causal.y = 0.01,
                          prop.causal.ylx = 1, 
                          prop.variance.y = 0.1,
                          prop.variance.x = 0.1, 
                          rho = 0.1, 
                          sigma = 1, 
                          sd.A = 0.2, 
                          mean.A = 1, 
                          sd.B = 0.2, 
                          mean.B = 1, 
                          sd.U = 1,
                          sd.V = 1,
                          prop.CT = NULL) 
{
  causal.x <- sample.int(p, prop.causal.x * p)
  causal.ylx <- sample(causal.x, prop.causal.ylx * length(causal.x))
  if (prop.causal.y * p < prop.causal.ylx * length(causal.x)) {
    stop("# causal y < # mediators")
  }
  else {
    causal.y <- c(causal.ylx, sample.int(p, prop.causal.y * 
                                           p - prop.causal.ylx * length(causal.x)))
  }
  x.nb = length(causal.x)
  y.nb = length(causal.y)
  if (is.null(freq)) 
    freq <- runif(n = p, min = 0.2, max = 0.8)
  if (prop.variance.y + rho^2 > 1) 
    stop("prop.variance.y + rho^2 > 1")
  if (prop.variance.x + rho^2 > 1) 
    stop("prop.variance.x + rho^2 > 1")
  cs.y <- runif(K, min = -1, max = 1)
  cs.x <- runif(K, min = -1, max = 1)
  theta.y <- sqrt(prop.variance.y/sum((cs.y/sd.U)^2))
  theta.x <- sqrt(prop.variance.x/sum((cs.x/sd.U)^2))
  Sigma <- diag(x = sd.U^2, nrow = K, ncol = K)
  Sigma <- rbind(Sigma, matrix(cs.y * theta.y, nrow = 1))
  Sigma <- rbind(Sigma, matrix(cs.x * theta.x, nrow = 1))
  Sigma <- cbind(Sigma, matrix(c(cs.y * theta.y, 1, rho), ncol = 1))
  Sigma <- cbind(Sigma, matrix(c(cs.x * theta.x, rho, 1), ncol = 1))
  UYX <- MASS::mvrnorm(n, mu = rep(0, K + 2), Sigma = Sigma)
  U <- UYX[, 1:K, drop = FALSE]
  Y <- UYX[, K + 1, drop = FALSE]
  X <- UYX[, K + 2, drop = FALSE]
  V <- MASS::mvrnorm(p, mu = rep(0, K), Sigma = sd.V^2 * diag(K))
  A <- matrix(0, p, 1)
  A[causal.x, 1] <- rnorm(x.nb, mean.A, sd.A)
  B <- matrix(0, p, 1)
  B[causal.y, 1] <- rnorm(y.nb, mean.B, sd.B)
  Epsilon <- apply(matrix(rep(0, p), nrow = 1), 2, function(x) rnorm(n, 
                                                                     x, sigma))
  
  if(is.null(prop.CT)) {
    prop.CT <- runif(K.ct, min = 0.05, max = 0.95) 
    prop.CT <- prop.CT / sum(prop.CT) 
  }
  
  CT <- gtools::rdirichlet(n = n, alpha = prop.CT)
  TC <- gtools::rdirichlet(n = p, alpha = prop.CT)
  
  Z = CT %*% t(TC) + U %*% t(V) + X %*% t(A) + Y %*% t(B) + Epsilon
  M = matrix(rep(qnorm(freq), n), nrow = n, byrow = T) + Z
  M = pnorm(M)
  return(list(M = M, 
              Y = Y, 
              B = B, 
              X = X, 
              A = A, 
              mediators = sort(causal.ylx), 
              causal.x = sort(causal.x), 
              causal.y = sort(causal.y), 
              U = U, 
              V = V, 
              CT = CT,
              TC = TC, 
              K.ct = K.ct,
              prop.CT = prop.CT,
              K = K,
              freq = freq, 
              Sigma = Sigma, 
              controls = !(1:p %in% unique(sort(c(sort(causal.x), sort(causal.y)))))))
}


r_mediation3 <- function (n = 400, 
                          p = 10000,  
                          K = 2, 
                          K.ct = 5, 
                          freq = NULL, 
                          prop.causal.x = 0.01, 
                          prop.causal.y = 0.01,
                          prop.causal.ylx = 1, 
                          prop.variance.y = 0.1,
                          prop.variance.x = 0.1, 
                          rho = 0.1, 
                          sigma = 1, 
                          sd.A = 0.2, 
                          mean.A = 1, 
                          sd.B = 0.2, 
                          mean.B = 1, 
                          sd.U = 1,
                          sd.V = 1,
                          prop.CT = NULL, 
                          only.CT = T) 
{
  causal.x <- sample.int(p, prop.causal.x * p)
  causal.ylx <- sample(causal.x, prop.causal.ylx * length(causal.x))
  if (prop.causal.y * p < prop.causal.ylx * length(causal.x)) {
    stop("# causal y < # mediators")
  }
  else {
    causal.y <- c(causal.ylx, sample.int(p, prop.causal.y * 
                                           p - prop.causal.ylx * length(causal.x)))
  }
  
  x.nb = length(causal.x)
  y.nb = length(causal.y)
  
  if (only.CT) {
    K <- 1
    prop.variance.y <- 0
    prop.variance.x <- 0
  }
  
  if (is.null(freq)) 
    freq <- runif(n = p, min = 0.2, max = 0.8)
  
  if (prop.variance.y + rho^2 > 1) 
    stop("prop.variance.y + rho^2 > 1")
  
  if (prop.variance.x + rho^2 > 1) 
    stop("prop.variance.x + rho^2 > 1")
  
  cs.y <- runif(K, min = -1, max = 1)
  cs.x <- runif(K, min = -1, max = 1)
  
  theta.y <- sqrt(prop.variance.y/sum((cs.y/sd.U)^2))
  theta.x <- sqrt(prop.variance.x/sum((cs.x/sd.U)^2))
  
  Sigma <- diag(x = sd.U^2, nrow = K, ncol = K)
  Sigma <- rbind(Sigma, matrix(cs.y * theta.y, nrow = 1))
  Sigma <- rbind(Sigma, matrix(cs.x * theta.x, nrow = 1))
  Sigma <- cbind(Sigma, matrix(c(cs.y * theta.y, 1, rho), ncol = 1))
  Sigma <- cbind(Sigma, matrix(c(cs.x * theta.x, rho, 1), ncol = 1))
  
  UYX <- MASS::mvrnorm(n, mu = rep(0, K + 2), Sigma = Sigma)
  
  U <- UYX[, 1:K, drop = FALSE]
  Y <- UYX[, K + 1, drop = FALSE]
  X <- UYX[, K + 2, drop = FALSE]
  
  V <- MASS::mvrnorm(p, mu = rep(0, K), Sigma = sd.V^2 * diag(K))
  
  A <- matrix(0, p, 1)
  A[causal.x, 1] <- rnorm(x.nb, mean.A, sd.A)
  
  B <- matrix(0, p, 1)
  B[causal.y, 1] <- rnorm(y.nb, mean.B, sd.B)
  
  Epsilon <- apply(matrix(rep(0, p), nrow = 1), 2, function(x) rnorm(n, 
                                                                     x, sigma))
  
  if(is.null(prop.CT)) {
    prop.CT <- runif(K.ct, min = 0.05, max = 0.95) 
    prop.CT <- prop.CT / sum(prop.CT) 
  }
  
  CT <- gtools::rdirichlet(n = n, alpha = prop.CT)
  TC <- gtools::rdirichlet(n = p, alpha = prop.CT)
  
  if (only.CT) {
    
    Z <- CT %*% t(TC) + X %*% t(A) + Y %*% t(B) + Epsilon
    M <- matrix(rep(qnorm(freq), n), nrow = n, byrow = T) + Z
    M <- pnorm(M)
    U <- NULL
    V <- NULL
    
  } else {
    
    Z <- CT %*% t(TC) + U %*% t(V) + X %*% t(A) + Y %*% t(B) + Epsilon
    M <- matrix(rep(qnorm(freq), n), nrow = n, byrow = T) + Z
    M <- pnorm(M)
    
  }
  
  return(list(M = M, 
              Y = Y, 
              B = B, 
              X = X, 
              A = A, 
              mediators = sort(causal.ylx), 
              causal.x = sort(causal.x), 
              causal.y = sort(causal.y), 
              U = U, 
              V = V, 
              CT = CT,
              TC = TC, 
              K.ct = K.ct,
              prop.CT = prop.CT,
              K = K,
              freq = freq, 
              Sigma = Sigma, 
              controls = !(1:p %in% unique(sort(c(sort(causal.x), sort(causal.y)))))))
}

wrap_multimediate <- function(qval, X , Y, M, covar = NULL, U = NULL, FDR = 0.1, sims = 10) {
  
  m <- M[, qval <= FDR]
  
  CpG <- c("Joint", colnames(m))
  colnames(m) <- NULL
  
  if (is.null(covar)) {
    data1 <- data.frame(M = m, Treatment = X, U = U)
  } else {
    data1 <- data.frame(M = m, Treatment = X, U = U, C = covar)
  }
  
  llm <- list()
  
  for (i in 1:ncol(m)) {
    regression <- paste0(paste0("M.", i), " ~ Treatment + .")
    llm[[i]] <- lm(as.formula(regression), data = data1)
  }
  
  # names(llm) <- paste0("mod.", 1:ncol(m))
  
  if (is.null(covar)) {
    data1 <- data.frame(M = m, Treatment = X, U = U, Outcome = Y)
  } else {
    data1 <- data.frame(M = m, Treatment = X, U = U, C = covar, Outcome = Y)
  }
  
  Yreg <- lm(Outcome ~ Treatment + ., data = data1)
  
  # summary(Yreg)
  
  med.analysis <- multimediate::multimediate(lmodel.m = llm, 
                                             correlated = T, model.y = Yreg, 
                                             treat = "Treatment", J = sims, conf.level = 0.95)
  
  tmp <- summary(med.analysis, opt = "avg")
  
  acme <- tmp[-c((nrow(tmp) - 1), nrow(tmp)), ] 
  acme <- acme[seq(1, nrow(acme), 2), ]
  
  pm <- tmp[-c((nrow(tmp) - 1), nrow(tmp)), ]
  pm <- pm[seq(2, nrow(pm), 2), ]
  
  acme <- data.frame(CpG, acme)
  pm <- data.frame(CpG, pm)
  
  return(list(ACME = acme[, -2], PM = pm[, -2]))
  
}