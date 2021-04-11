# Parameter calculation program for simulations of Maruo et al. (2017). 
# Stat. Med. 36, 15, 2420-2434
# Setting: power normal distributions (PN)

#install.packages("Matrix")
#install.packages("doParallel")
#install.packages("nlme")
#install.packages("bcmixed)
library(doParallel)

#Reparametrization of the power normal distribution (PND)
PND_reparm <- function(lambda, xi, tau){
  xik_ms <- function(lambda, xi, K){
    z5 <- zps(lambda, K, 0.5)
    mu <- (K * (xi ^ lambda - 1) - z5) / (lambda * (K + z5))
    sigma <- (1 + lambda * mu) / (lambda * K)
    return(list(mu, sigma))
  }
  zps <- function(lambda, K, p){
    AK <- pnorm(sign(lambda) * K)
    if (lambda < 0) ps <- AK * p
    if (lambda > 0) ps <- 1 - AK * (1 - p)
    zps <- qnorm(ps)
    return(zps)
  }
  if (abs(lambda) > 0.01){
    K1 <- sign(lambda) * seq(0.1, 100, 0.1)
    kp <- length(K1)
    ms1 <- xik_ms(lambda, xi, K1)
    mu1 <- ms1[[1]]
    sigma1 <- ms1[[2]]
    d1 <- pnd_quantile(lambda, mu1, sigma1, 0.75) - 
      pnd_quantile(lambda, mu1, sigma1, 0.25) - xi * tau
    ld1 <- which.min(abs(d1))
    if (sign(d1[ld1 - 1]) == sign(d1[ld1])) {
      K2 <- seq(K1[ld1], K1[ld1 + 1], sign(lambda) * 1e-4)
    } else {
      K2 <- seq(K1[ld1 - 1], K1[ld1], sign(lambda) * 1e-4)
    }
    kp <- length(K2)
    ms2 <- xik_ms(lambda, xi, K2)
    mu2 <- ms2[[1]]
    sigma2 <- ms2[[2]]
    d2 <- pnd_quantile(lambda, mu2, sigma2, 0.75) - 
      pnd_quantile(lambda, mu2, sigma2, 0.25) - xi * tau
    ld2 <- which.min(abs(d2))
    mu <- mu2[ld2]
    sigma <- sigma2[ld2]
  }
  else{
    mu <- log(xi)
    sigma <- log((tau + sqrt(tau ^ 2 + 4)) / 2) / qnorm(0.75)
  }
  return(c(mu, sigma))
}

# 100p percentile of PND(lambda, mu, sigma)
pnd_quantile <- function(lambda, mu, sigma, p){
  zps <- function(lambda, K, p){
    AK <- pnorm(sign(lambda) * K)
    if (lambda < 0) ps <- AK * p
    if (lambda > 0) ps <- 1 - AK * (1 - p)
    zps <- qnorm(ps)
    return(zps)
  }
  if(lambda != 0){
    K <- (1 + lambda * mu) / (lambda * sigma)
    zps <- zps(lambda, K, p)
    q <- (lambda * (mu + sigma * zps) + 1) ^ (1 / lambda)
  }
  else q <- exp(mu + sigma * qnorm(p))
  return(q)
}

parmset <- function(lmd, tau, npg, cov, sgn, drop){
  
  n <- 200000
  tn <- 3
  
  ms1 <- PND_reparm(lmd, 100, tau * 0.8)
  ms2 <- PND_reparm(lmd, 100-sgn * 10, tau * 0.9)
  ms3 <- PND_reparm(lmd, 100-sgn * 20, tau)
  m0 <- c(ms1[1], ms1[1], ms1[2], ms1[3])
  rho <- 0.7
  R <- matrix(0, tn + 1, tn + 1)
  for (i in 1:(tn + 1)){
    for (j in 1:(tn + 1)){
      R[i, j] <- rho^abs(i-j)
    }
  }
  y0 <- mvrnorm(n, numeric(tn + 1), R)
  
  y0[, 1:2] <- y0[, 1:2] * ms1[2] + ms1[1] 
  y0[, 3] <- y0[, 3] * ms2[2] + ms2[1] 
  y0[, 4] <- y0[, 4] * ms3[2] + ms3[1] 
  ms_s <- c(ms1[2], ms2[2], ms3[2]) 
  xi1 <- c(100, 100 - sgn * 10, 100 - sgn * 20)
  if (lmd < 0){
    zup0 <- t(((c(numeric(4) + xi1[1]) * 50) ^ lmd-1) / lmd) %x% 
      (numeric(nrow(y0)) + 1)
    t000 <- (y0 < -1 / lmd) & (y0 < zup0)
  }
  if (lmd > 0){
    t000 <- y0 > -1 / lmd
  }
  if (lmd  ==  0){
    t000 <- y0 > -Inf
  }
  
  y0 <- y0[t000[, 1] & t000[, 2] & t000[, 3] & t000[, 4], ]
  
  #drop <- 0.5
  
  drop_f <- function(Ilgt, y0, runi, ms_s, tn, sgn, drop){
    lgt <- c()
    for (i in 2:tn){
      lgt <- cbind(lgt, Ilgt - sgn * y0[, i] / ms_s[i - 1])
    }
    elgt <- exp(lgt)
    pdp <- elgt / (1 + elgt)
    drp <- c()
    for (i in 1:(tn - 1)){
      drp <- cbind(drp, 1 * (runif(nrow(y0)) < pdp[, i]))
    }
    y01 <- y0
    
    for (i in 3:(tn + 1)){
      y01[is.na(y01[, i - 1])|drp[, i - 2] == 1, i] <- NA
    }
    return(mean(is.na(y01[, tn + 1])) - drop)  
  }
  
  if (drop>0) {
    Ilgt <- uniroot(drop_f, y0 = y0, runi = runi, ms_s = ms_s, tn = tn, 
                    sgn = sgn, drop = drop, interval = c(-100, 100))$root
  } else {
    Ilgt <- -Inf
  }
  
  dlt <- sqrt(2 / (npg * (1 - drop) * 1.1)) * ms3[2] * 
    (qnorm(0.975) + qnorm(0.8))
  return(list(lmd = lmd, npg = npg, cov = cov, ms1 = ms1, ms2 = ms2, ms3 = ms3, 
              Ilgt = Ilgt, dlt = sgn * dlt, drop = drop))
}


# Data generation
rand_parm_bc_sim <- function(lmd, npg, cov, m1, m2, m3, s1, s2, s3, Ilgt, dlt, h){
  ms1 <- c(m1, s1)
  ms2 <- c(m2, s2)
  ms3 <- c(m3, s3)
  sgn <- sign(dlt)
  xi1 <- c(100, 100 - sgn * 10, -sgn * 20)
  tn <- 3
  rho <- 0.7
  ms_s <- c(ms1[2], ms2[2], ms3[2])
  R <- matrix(0, tn + 1, tn + 1)
  for (i in 1:(tn + 1)){
    for (j in 1:(tn + 1)){
      R[i, j] <- rho ^ abs(i - j)
    }
  }
  y00 <- mvrnorm(npg * 4, numeric(tn + 1), R)
  xg0 <- 1 * ((1:npg * 4) > npg * 2)
  
  y00[, 1] <- y00[, 1] * ms1[2] + ms1[1]
  if (h == 0){
    y00[, 2] <- y00[, 2] * ms1[2] + ms1[1]
    y00[, 3] <- y00[, 3] * ms2[2] + ms2[1] 
    y00[, 4] <- y00[, 4] * ms3[2] + ms3[1]
    
  }
  if (h == 1){
    y00[, 2] <- y00[, 2] * ms1[2] + ms1[1] + xg0 * dlt / 2
    y00[, 3] <- y00[, 3] * ms2[2] + ms2[1] + xg0 * dlt / 2
    y00[, 4] <- y00[, 4] * ms3[2] + ms3[1]
  }
  if (h == 2){
    y00[, 2] <- y00[, 2] * ms1[2] + ms1[1] + xg0 * dlt / 2
    y00[, 3] <- y00[, 3] * ms2[2] + ms2[1] + xg0 * dlt 
    y00[, 4] <- y00[, 4] * ms3[2] + ms3[1] + xg0 * dlt
  }
  
  y000 <- y00[xg0 == 0, ]
  y001 <- y00[xg0 == 1, ]
  
  if (lmd < 0){
    zup0 <- t(((c(numeric(4) + xi1[1]) * 50) ^ lmd - 1) / lmd) %x% 
      (numeric(nrow(y000)) + 1)
    zup1 <- t(((c(numeric(4) + xi1[1]) * 50) ^ lmd - 1) / lmd) %x% 
      (numeric(nrow(y001)) + 1)
    t000 <- (y000 < -1 / lmd) & (y000 < zup0)
    t001 <- (y001 < -1 / lmd) & (y001 < zup1)
  }
  if (lmd > 0){
    t000 <- y000 > -1 / lmd
    t001 <- y001 > -1 / lmd
  }
  if (lmd == 0){
    t000 <- y000 > -Inf
    t001 <- y001 > -Inf
  }
  
  y000 <- y000[t000[, 1] & t000[, 2] & t000[, 3] & t000[, 4], ]
  y001 <- y001[t001[, 1] & t001[, 2] & t001[, 3] & t001[, 4], ]
  y000 <- y000[1:npg, ]
  y001 <- y001[1:npg, ]
  y00 <- rbind(y000, y001)
  lgt <- c()
  for (i in 2:tn){
    lgt <- cbind(lgt, Ilgt - sgn * y00[, i] / ms_s[i - 1])
  }
  elgt <- exp(lgt)
  pdp <- elgt / (1 + elgt)
  
  drp <- c()
  for (i in 1:(tn - 1)){
    drp <- cbind(drp,  1 * (runif(npg * 2) < pdp[, i]))
  }
  
  for (i in 3:(tn + 1)){
    y00[is.na(y00[, i - 1])|drp[, i - 2] == 1, i] <- NA
  }
  
  z <- c(y00[1:npg, 2:(tn + 1)], y00[(npg + 1):(2 * npg), 2:(tn + 1)])
  if (lmd == 0){  y <- exp(z)} else {  y <- (lmd * z + 1) ^ (1 / lmd)}
  zbl <- c((numeric(tn) + 1) %x% y00[1:npg, 1], 
           (numeric(tn) + 1) %x% y00[(npg + 1):(2 * npg), 1])
  if (lmd == 0){bl <- exp(zbl)}else {bl <- (lmd * zbl + 1) ^ (1 / lmd)}
  time <- c(1, 1) %x% ((1:tn) %x% (numeric(npg) + 1))
  group <- c(0, 1) %x% (numeric(npg * tn) + 1)
  id <- c(1, 1) %x% (numeric(tn) + 1) %x% (1:npg)
  id <- paste(group, "-", id)
  dat <- data.frame(y = y, bl = bl, time = time, id = id, group = group, zt = z, 
                    lblt = zbl)
  
  
  dat0 <- dat[order(dat$id, dat$time), ]
  
  l_by <- bcreg(bl ~ 1, data = dat0)$lambda
  
  dat0$lbl <- (dat0$bl ^ l_by - 1) / l_by
  return(dat0)
}

# MMRM analysis
MMRM <- function(outcome, group, time, id, timepoint, data, covv = NULL, 
                 cfactor = NULL){
  tr <- function(X) sum(diag(X)) 
  corandcov <- function(glsob, nt){
    flg <- 0
    j <- 0
    while (flg == 0){
      j <- j + 1
      corm <- corMatrix(glsob$modelStruct$corStruct)[[j]]
      if (nrow(corm) == nt) flg <- 1
    }
    varstruct <- glsob$modelStruct$varStruct  
    varests <- coef(varstruct, uncons = F, allCoef = T)
    covm <- corm * glsob$sigma ^ 2 * t(t(varests)) %*% t(varests)
    return(covm)
  }
  y <- deparse(substitute(outcome))
  group <- deparse(substitute(group))
  time <- deparse(substitute(time))
  cov <- deparse(substitute(cov))
  id <- deparse(substitute(id))
  data$y <- data[, y]
  data$group <- data[, group]
  data$time0 <- data[, time]
  data$id <- data[, id]
  datatable <- sort(unique(data$time0))
  for (i in 1:length(datatable)){
    data[data$time0 == datatable[i], "time"] <- i
    if (datatable[i] == timepoint) tp0 <- i
  }
  covform <- ""
  cct <- 0
  cfactor <- c(0, 1, 0)
  if (!is.null(covv)){
    cct <- numeric(length(covv)) + 1
    for (i in 1:length(covv)){
      if (cfactor[i] == 0) covform <- paste(covform, "+", covv[i])
      else {
        covform <- paste(covform, "+as.factor(", covv[i], ")")
        cct[i] <- length(unique(data[, covv[i]])) - 1
      }
    }
  }
  form <- paste("y~ as.factor(group) + as.factor(time) + as.factor(group):as.factor(time)", covform)
  ng <- length(unique(data$group))
  nt <- length(unique(data$time))
  nc <- sum(cct)
  UN <- gls(as.formula(form), data = data, corr = corSymm(form =~ time|id),  
            weights = varIdent(form =~ 1|time), method = "REML",  
            control = lmeControl(msMaxIter = 100, msVerbose = F), 
            na.action = na.omit)
  V <- corandcov(UN, nt)
  beta <- UN$coefficients
  ntm <- nrow(V)
  alp <- c()
  for (i in 1:ntm){
    for (j in i:ntm){
      alp <- c(alp, V[i, j])
    }
  }
  nv <- length(alp)
  Vi <- ginv(V)
  N <- length(table(data$id))
  
  nna <- !is.na(data$y)
  nn <- sum(nna)
  y <- data$y[nna]
  options(na.action = "na.pass")
  X <- model.matrix(as.formula(form), data = data)
  dimnames(X) <- NULL
  Va <- diag(N) %x% V
  Va <- Va[nna, nna]
  nti <- table(data$id[!is.na(data$y)])
  ntic <- cumsum(nti)
  ntic <- c(0, ntic)
  Via <- list()
  for (j in 1:N){
    jj <- (ntic[j] + 1):ntic[j + 1]
    Via[[j]] <- ginv(Va[jj, jj])
  }
  xcm0 <- 0
  xcm <- 0
  if (!is.null(covv)){
    if (nc == 1)  {
      xcm0 <- mean(X[data$time == unique(data$time)[1], ng + nt])
    } else {
      xcm0 <- c(colMeans(X[data$time == unique(data$time)[1], 
                           (ng + nt):(ng + nt + nc - 1)]))
    }
    xcm <- c()
    for (i in 1:length(covv)){
      if (cfactor[i] == 0) {
        xcm <- c(xcm, xcm0[i])
      } else {
        nci <- length(unique(data[, covv[i]]))
        xcm <- c(xcm, numeric(nci - 1) + 1 / nci)
      }
    }
  } 
  X <- X[nna, ]
  Vl <- list()
  for (i in 1:nv){
    vl0 <- diag(N) %x% (1 * (V == alp[i]))
    Vl[[i]] <- vl0[nna, nna]
  } 
  AA <- list()
  for (i in 1:nv){
    BB <- list()
    for (j in 1:N){
      jj <- (ntic[j] + 1):ntic[j + 1]
      BB[[j]] <- Via[[j]] %*% Vl[[i]][jj, jj] %*% Via[[j]]
    }
    AA[[i]] <- BB
  } 
  re <- y - X %*% beta
  l_bb <- as.matrix(-t(X) %*% bdiag(Via) %*% X)
  C <- ginv(-l_bb)
  XC <- X %*% C
  Ci <<- C
  Xd <- X %*% t(chol(C + diag(nrow(C)) * diag(diag(C)) * 1e-8))
  H2i <- list()
  H3i <- list()
  Viaa <- as.matrix(bdiag(Via))
  Viaa <- bdiag(Via)
  ViXd <- Viaa %*% Xd
  Vire <- Viaa %*% re
  VlVire <- list()
  VlViXd <- list()
  for (i in 1:nv){
    H2i[[i]] <- t(ViXd) %*% Vl[[i]] %*% Vire
    VlVire[[i]] <- Vl[[i]] %*% Vire
    H3i[[i]] <- t(ViXd) %*% Vl[[i]] %*% ViXd
    VlViXd[[i]] <- Vl[[i]] %*% ViXd
  }  
  Hv <- matrix(0, nv, nv)
  for (i in 1:nv){
    for (j in 1:nv){
      if (i <= j){
        BB <- list()
        for (k in 1:N){
          jj <- (ntic[k] + 1):ntic[k + 1]
          BB[[k]] <- AA[[i]][[k]] %*% Vl[[j]][jj, jj]
        }
        H2ij <- 2 * t(VlVire[[i]]) %*% Viaa %*% VlVire[[j]]
        H3ij <- 2 * t(VlViXd[[i]]) %*% Viaa %*% VlViXd[[j]]
        H1 <- tr(-bdiag(BB))
        H2 <- as.numeric(H2ij - 2 * t(H2i[[i]]) %*% H2i[[j]])
        H3 <- tr(as.matrix(H3ij - H3i[[i]] %*% H3i[[j]]))
        Hv[i, j] <- -1 / 2 * (H1 + H2 + H3)
      } else{
        Hv[i, j] <- Hv[j, i]
      } 
    }
  }
  #print("e")
  lsm <- numeric(ng)
  dbt <- numeric(nt - 1)
  iIIv <- ginv(-Hv)
  
  if (tp0 != 1){
    dbt[tp0 - 1] <- 1
  } 
  SElsm <- numeric(ng)
  lsm <- SElsm
  ell <- matrix(0, length(beta), ng)
  for (i in 1:ng){
    dbg <- numeric(ng - 1)
    if (i != 1) {
      dbg[i - 1] <- 1
    }
    dbgt <- numeric((ng - 1) * (nt - 1))
    if (i != 1 & tp0 != 1){
      bgti <- (tp0 - 2) * (ng - 1) + i - 1
      dbgt[bgti] <- 1 
    } 
    if (!is.null(covv)) ell[, i] <- c(1, dbg, dbt, xcm, dbgt)
    else ell[, i] <- c(1, dbg, dbt, dbgt)
    lsm[i] <- t(ell[, i]) %*% beta
    Vlsm <- c(t(ell[, i]) %*% C %*% ell[, i])
    SElsm[i] <- sqrt(Vlsm)
  }
  
  cbn <- choose(ng, 2)
  comb <- matrix(0, cbn, 2)
  count <- 1
  for (ii in 1:ng){
    for (jj in 1:ng){
      if (ii<jj) {
        comb[count, 1] <- ii
        comb[count, 2] <- jj
        count <- count + 1
      }
    }
  }
  FF <- list()
  for (i in 1:nv){
    FF0 <- 0
    for (k in 1:N){
      jj <- (ntic[k] + 1):ntic[k + 1]
      Xjj <- X[jj, ]
      if (length(jj) == 1) Xjj <- t(Xjj)
      FF00 <- Via[[k]] %*% Xjj %*% C
      FF0 <- FF0 + t(FF00) %*% Vl[[i]][jj, jj] %*% FF00
    }
    FF[[i]] <- FF0
  }
  
  delta <- numeric(cbn)
  SEd <- numeric(cbn)
  df <- numeric(cbn)
  for (i in 1:cbn){
    delta[i] <- lsm[comb[i, 2]] - lsm[comb[i, 1]]
    elld <- ell[, comb[i, 2]] - ell[, comb[i, 1]]
    Vd <- c(t(elld) %*% C %*% elld)
    SEd[i] <- sqrt(Vd)
    gg <- numeric(nv)
    for (j in 1:nv){
      gg[j] <- as.numeric(t(elld) %*% FF[[j]] %*% elld)
    }
    df[i] <- c(2 * Vd ^ 2 / t(gg) %*% iIIv %*% gg)
  }
  tvalue <- delta / SEd
  alpha <- 0.05
  tt <- qt(1 - alpha / 2, df)
  pvalue <- (1 - pt(abs(tvalue), df)) * 2
  lower <- delta - SEd * tt
  upper <- delta + SEd * tt
  lsmeans <- data.frame(group = 1:ng, estimate = lsm, SE = SElsm) 
  lsmeans_diff <- data.frame(group1 = comb[, 2], group0 = comb[, 1], estimate = delta, 
                             SE = SEd, lower_CL = lower, upper_CL = upper, df = df, 
                             tvalue = tvalue, pvalue = pvalue)
  result <- list(lsmeans = lsmeans, lsmeans_diff = lsmeans_diff)
  return(result)
}

simset <- function(cov, sgn, lmd, drop){
  tau <- 1
  npg <- 50
  
  parm <- parmset(lmd, tau, npg, cov, sgn, drop)

  
  dlt <- parm$dlt
  sgm <- parm$ms3[2]
  
  dlta <- sign(dlt) * seq(sign(dlt) * dlt * 0.5, sign(dlt) * dlt * 1.8, 
                        by = sgm / 18)
  
  resm <- c()
  ms1 <- parm$ms1
  ms2 <- parm$ms2
  ms3 <- parm$ms3
  for (dlt0 in dlta){
    resa <- c()
    for (i in 1:200){
      dat0 <- rand_parm_bc_sim(lmd, npg, cov, ms1[1], ms2[1], ms3[1], 
                               ms1[2], ms2[2], ms3[2], parm$Ilgt, dlt0, 2)
      datt <- dat0
      datt$y <- dat0$zt
      datt$bl <- dat0$lblt
      covv <- NULL
      cfactor <- NULL
      if (cov == 0) {
        res <- c(MMRM(outcome = y, group = group, time = time, id = id, 
                      timepoint = 3, data = datt)$lsmeans_diff$pvalue < 0.05)
      } else {
        res <- c(MMRM(outcome = y, group = group, time = time, id = id, 
                      timepoint = 3, data = datt, covv = "bl", 
                      cfactor = 0)$lsmeans_diff$pvalue < 0.05)
      }
      resa <- c(resa, res)
    }
    resm <- rbind(resm, cbind(dlt0, resa))
  }
  resm <- as.data.frame(resm)
  pcf <- glm(resa ~ dlt0, family = binomial(link = "probit"), 
             data = resm)$coefficients
  dltt <- (qnorm(0.8) - pcf[1]) / pcf[2]
  Ilgt <- parm$Ilgt
  xi2 <- c(pnd_quantile(lmd, ms1[1] + dltt / 2, ms1[2], .5),
           pnd_quantile(lmd, ms2[1] + dltt, ms2[2], .5),
           pnd_quantile(lmd, ms3[1] + dltt, ms3[2], .5))
  datl <- rand_parm_bc_sim(lmd, 500000, cov, ms1[1], ms2[1], ms3[1], 
                           ms1[2], ms2[2], ms3[2], Ilgt, dltt, 2)
  rr1 <- c(1 - mean(is.na(datl$y[datl$time == 2 & datl$group == 0])), 
           1 - mean(is.na(datl$y[datl$time == 3 & datl$group == 0])))
  rr2 <- c(1 - mean(is.na(datl$y[datl$time == 2 & datl$group == 1])), 
           1 - mean(is.na(datl$y[datl$time == 3 & datl$group == 1])))
  return(c(cov, sgn, lmd, drop, tau, npg, ms1, ms2, ms3, Ilgt, dltt, 
           xi2, rr1, rr2))
}

conda <- c()
for (cov in c(0, 1)){
  for (sgn in c(-1, 1)){
    for (lmd in c(-0.5, 0, 0.5)){
      for (drop in c(0, 0.3, 0.5)){
        conda <- rbind(conda, c(cov, sgn, lmd, drop))
      }
    }
  }
}
conda <- as.data.frame(conda)
names(conda) <- c("cov", "sgn", "lmd", "drop")

gc()
t<-proc.time()

cluster = makeCluster(40)
clusterSetRNGStream(cluster, 12345)
registerDoParallel(cluster)
simr <- foreach (cov = conda$cov, sgn = conda$sgn, lmd = conda$lmd, 
                 drop = conda$drop, .combine = rbind, 
                 .packages = c("nlme", "Matrix", "MASS","bcmixed")) %dopar%  
  {
    simset(cov, sgn, lmd, drop)
  }
stopCluster(cluster)
simr


as.data.frame(simr)

colnames(simr) <- c("cov", "sgn", "lmd", "drop", "tau", "npg", "m1", "s1", "m2", 
                    "s2", "m3", "s3", "Ilgt", "dltt", "xi2_1", "xi2_2", "xi2_3", 
                    "drrate1_2", "drrate1_3", "drratte2_2", "drrate2_3")
proc.time()-t
write.csv(simr,"Simulation_bcmixed_PN_parm.csv",row.names=F)
