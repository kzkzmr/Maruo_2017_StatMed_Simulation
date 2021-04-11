# Parameter calculation program for simulations of Maruo et al. (2017). 
# Stat. Med. 36, 15, 2420-2434
# Setting: Weibull and gamma normal distributions (WG)

#install.packages("Matrix")
#install.packages("doParallel")
#install.packages("nlme")
#install.packages("bcmixed)
library(doParallel)

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



parmset_wg <- function(k, sgn, wg, drop){
  parmset_w <- function(k, sgn, drop){
    
    n <- 200000
    tn <- 3
    
    th1_1 <- 100 / (log(2) ^ (1 / k))
    th1_2 <- (100 - sgn*10) / (log(2) ^ (1 / k))
    th1_3 <- (100 - sgn*20) / (log(2) ^ (1 / k))
    rho <- 0.7
    R <- matrix(0, tn + 1, tn + 1)
    for (i in 1:(tn + 1)){
      for (j in 1:(tn + 1)){
        R[i, j] <- rho ^ abs(i - j)
      }
    }
    c0 <- pnorm(mvrnorm(n, numeric(tn + 1), R))
    y0 <- c0
    y0[, 1:2] <- qweibull(c0[, 1:2], k, scale = th1_1) 
    y0[, 3] <- qweibull(c0[, 3], k, scale = th1_2) 
    y0[, 4] <- qweibull(c0[, 4], k, scale = th1_3) 
    sd0 <- sqrt(gamma(1 + 2 / k) - gamma(1 + 1 / k) ^ 2) 
    sd <- c(th1_1*sd0, th1_2*sd0, th1_3*sd0) 
    xi1 <- c(100, 100 - sgn*10, 100 - sgn*20)
    
    drop_f <- function(Ilgt, y0, runi, sd, tn, sgn, drop){
      lgt <- c()
      for (i in 2:tn){
        lgt <- cbind(lgt, Ilgt - sgn*y0[, i] / sd[i - 1])
      }
      elgt <- exp(lgt)
      pdp <- elgt / (1 + elgt)
      drp <- c()
      for (i in 1:(tn - 1)){
        drp <- cbind(drp,  1*(runif(nrow(y0)) < pdp[, i]))
      }
      y01 <- y0
      
      for (i in 3:(tn + 1)){
        y01[is.na(y01[, i - 1])|drp[, i - 2]==1, i] <- NA
      }
      return(mean(is.na(y01[, tn + 1])) - drop)  
    }
    th1 <- c(th1_1, th1_2, th1_3)
    if (drop==0){
      Ilgt <- -Inf
    } else {
      Ilgt <- uniroot(drop_f, y0 = y0, runi = runi, sd = sd, tn = tn, 
                      sgn = sgn, drop = drop, interval = c(-100, 100))$root
    }
    return(list(th1 = th1, Ilgt = Ilgt, sd = sd))
  }
  
  dg_theta <- function(x, p, k){
    dg <- function(x, p, k, t){
      return(p - pgamma(x, k, scale = t))
    }
    return(uniroot(dg, c(1e-5, 1000), x = x, p = p, k = k)$root)
  }
  
  parmset_g <- function(k, sgn, drop){
    
    n <- 200000
    tn <- 3
    
    th1_1 <- dg_theta(100, 0.5, k)
    th1_2 <- dg_theta(100 - sgn*10, 0.5, k)
    th1_3 <- dg_theta(100 - sgn*20, 0.5, k)
    rho <- 0.7
    R <- matrix(0, tn + 1, tn + 1)
    for (i in 1:(tn + 1)){
      for (j in 1:(tn + 1)){
        R[i, j] <- rho ^ abs(i - j)
      }
    }
    c0 <- pnorm(mvrnorm(n, numeric(tn + 1), R))
    y0 <- c0
    y0[, 1:2] <- qgamma(c0[, 1:2], k, scale = th1_1) 
    y0[, 3] <- qgamma(c0[, 3], k, scale = th1_2) 
    y0[, 4] <- qgamma(c0[, 4], k, scale = th1_3) 
    sd <- c(sqrt(k)*th1_1, sqrt(k)*th1_2, sqrt(k)*th1_3) 
    xi1 <- c(100, 100 - sgn*10, 100 - sgn*20)
    

    drop_f <- function(Ilgt, y0, runi, sd, tn, sgn, drop){
      lgt <- c()
      for (i in 2:tn){
        lgt <- cbind(lgt, Ilgt - sgn*y0[, i] / sd[i - 1])
      }
      elgt <- exp(lgt)
      pdp <- elgt / (1 + elgt)
      drp <- c()
      for (i in 1:(tn - 1)){
        drp <- cbind(drp,  1*(runif(nrow(y0)) < pdp[, i]))
      }
      y01 <- y0
      
      for (i in 3:(tn + 1)){
        y01[is.na(y01[, i - 1])|drp[, i - 2]==1, i] <- NA
      }
      return(mean(is.na(y01[, tn + 1])) - drop)  
    }
    th1 <- c(th1_1, th1_2, th1_3)
    if (drop==0) {
      Ilgt <- -Inf
    } else {
      Ilgt <- uniroot(drop_f, y0 = y0, runi = runi, sd = sd, tn = tn, 
                      sgn = sgn, drop = drop, interval = c(-100, 100))$root
    }
    return(list(th1 = th1, Ilgt = Ilgt, sd = sd))
  }
  if (wg==0) return(parmset_w(k, sgn, drop))
  else return(parmset_g(k, sgn, drop))
}

rand_parm_wg <- function(cov, sgn, wg, k, npg, th1, Ilgt, sd, th2_3, h){
  q_wg <- function(p, k, th, wg){
    if (wg == 0) return(qweibull(p, k, scale = th))
    else return(qgamma(p, k, scale = th))
  }
  xi1 <- c(100, 100 - sgn * 10, 100 - sgn * 20)
  tn <- 3
  
  rho <- 0.7
  R <- matrix(0, tn + 1, tn + 1)
  for (i in 1:(tn + 1)){
    for (j in 1:(tn + 1)){
      R[i, j] <- rho ^ abs(i - j)
    }
  }
  y00 <- pnorm(mvrnorm(npg * 2, numeric(tn + 1), R))
  xg0 <- 1 * ((1:(npg * 2)) > npg)
  y000 <- y00[xg0 == 0, ]
  y001 <- y00[xg0 == 1, ]
  
  y00[, 1] <- q_wg(y00[, 1], k, th1[1], wg)
  if (h == 0){
    y00[, 2] <- q_wg(y00[, 2], k, th1[1], wg)
    y00[, 3] <- q_wg(y00[, 3], k, th1[2], wg)
    y00[, 4] <- q_wg(y00[, 4], k, th1[3], wg)
  }
  if (h == 1){
    th2 <- c(0, 0, 0)
    d3 <- th2_3 - th1[3]
    th2[1] <- th1[1] + d3/2
    th2[2] <- th1[2] + d3/2
    th2[3] <- th1[3]
    y00[xg0 == 0, 2] <- q_wg(y00[xg0 == 0, 2], k, th1[1], wg)
    y00[xg0 == 0, 3] <- q_wg(y00[xg0 == 0, 3], k, th1[2], wg)
    y00[xg0 == 0, 4] <- q_wg(y00[xg0 == 0, 4], k, th1[3], wg)
    y00[xg0 == 1, 2] <- q_wg(y00[xg0 == 1, 2], k, th2[1], wg)
    y00[xg0 == 1, 3] <- q_wg(y00[xg0 == 1, 3], k, th2[2], wg)
    y00[xg0 == 1, 4] <- q_wg(y00[xg0 == 1, 4], k, th2[3], wg)
  }
  if (h == 2){
    th2 <- c(0, 0, th2_3)
    d3 <- th2_3 - th1[3]
    th2[1] <- th1[1] + d3/2
    th2[2] <- th1[2] + d3
    y00[xg0 == 0, 2] <- q_wg(y00[xg0 == 0, 2], k, th1[1], wg)
    y00[xg0 == 0, 3] <- q_wg(y00[xg0 == 0, 3], k, th1[2], wg)
    y00[xg0 == 0, 4] <- q_wg(y00[xg0 == 0, 4], k, th1[3], wg)
    y00[xg0 == 1, 2] <- q_wg(y00[xg0 == 1, 2], k, th2[1], wg)
    y00[xg0 == 1, 3] <- q_wg(y00[xg0 == 1, 3], k, th2[2], wg)
    y00[xg0 == 1, 4] <- q_wg(y00[xg0 == 1, 4], k, th2[3], wg)
  }
  lgt <- c()
  for (i in 2:tn){
    lgt <- cbind(lgt, Ilgt - sgn * y00[, i] / sd[i - 1])
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
  
  y <- c(y00[1:npg, 2:(tn + 1)], y00[(npg + 1):(2 * npg), 2:(tn + 1)])
  
  bl <- c((numeric(tn) + 1) %x% y00[1:npg, 1], (numeric(tn) + 1) %x% y00[(npg + 1):(2 * npg), 1])
  time <- c(1, 1) %x% ((1:tn) %x% (numeric(npg) + 1))
  group <- c(0, 1) %x% (numeric(npg * tn) + 1)
  id <- c(1, 1) %x% (numeric(tn) + 1) %x% (1:npg)
  id <- paste(group,  '-',  id)
  dat <- data.frame(y = y, bl = bl, time = time, id = id, group = group)
  
  dat0 <- dat[order(dat$id, dat$time), ]
  
  l_by <- bcreg(bl ~ 1, data = dat0)$lambda
  dat0$lbl <- (dat0$bl ^ l_by - 1) / l_by
  return(dat0)
}


simset_wg <- function(cov, sgn, wg, k, npg, drop){
  parm <- parmset_wg(k, sgn, wg, drop)
  th1 <- parm$th1
  Ilgt <- parm$Ilgt
  sd <- parm$sd

  sgm <- sd[3] * (1 + sgn * 0.2)
  dlt <- sqrt(2 / npg / 0.8) * (qnorm(0.8) + qnorm(0.975)) * sgm
  dlta <- seq(dlt * 0.6, dlt * 2.5, by = sgm / 15)
  if (wg == 1) {
    dlta <- dlta * 0.6
  } else {
    dlta <- dlta * 1.1
  }
  resm <- c()
  resmm <- c()
  
  for (dlt0 in dlta){
    resa <- c()
    for (i in 1:1000){
      th2_3 <- th1[3] + sgn * dlt0
      dat0 <- rand_parm_wg(cov, sgn, wg, k, npg, th1, Ilgt, sd, th2_3, h = 2)
      if (cov == 0) {
        res <- c(MMRM(outcome = y, group = group, time = time, id = id, 
                      timepoint = 3, data = dat0)$lsmeans_diff$pvalue < 0.05)
      } else{
        res <- c(MMRM(outcome = y, group = group, time = time, id = id, 
                      timepoint = 3, data = dat0, covv = 'bl', 
                      cfactor = 0)$lsmeans_diff$pvalue < 0.05)
      }
      resa <- c(resa, res)
    }
    resmm <- c(resmm, mean(resa))
    resm <- rbind(resm, cbind(dlt0, resa))
  }
  resm <- as.data.frame(resm)
  dam <- 1:length(resmm)
  dlta1 <- dlta[min(dam[resmm > 0.5]):(min(dam[resmm > 0.9]) + 1)]
  cho <- numeric(nrow(resm)) == 1
  for (i in 1:length(dlta1)){
    cho <- cho|resm$dlt0 == dlta1[i]
  }
  resm <- resm[cho, ]
  pcf <- glm(resa ~ dlt0,  family = binomial(link = "probit"), 
             data = resm)$coefficients
  dltt <- (qnorm(0.8) - pcf[1]) / pcf[2]
  th2_3 <- th1[3] + sgn * dltt
  th2 <- c(0, 0, th2_3)
  d3 <- th2_3 - th1[3]
  th2[1] <- th1[1] + d3 / 2
  th2[2] <- th1[2] + d3
  q_wg <- function(p, k, th, wg){
    if (wg == 0) return(qweibull(p, k, scale = th))
    else return(qgamma(p, k, scale = th))
  }
  xi2 <- c(q_wg(0.5, k, th2[1], wg), q_wg(0.5, k, th2[2], wg), 
           q_wg(0.5, k, th2[3], wg))
  datl <- rand_parm_wg(cov, sgn, wg, k, 500000, th1, Ilgt, sd, th2_3, h = 2)
  rr1 <- c(1 - mean(is.na(datl$y[datl$time == 2 & datl$group == 0])), 
           1 - mean(is.na(datl$y[datl$time == 3 & datl$group == 0])))
  rr2 <- c(1 - mean(is.na(datl$y[datl$time == 2 & datl$group == 1])), 
           1 - mean(is.na(datl$y[datl$time == 3 & datl$group == 1])))
  return(c(cov, sgn, wg, k, drop, npg, th1, th2_3, sd, Ilgt, xi2, rr1, rr2))
}


conda <- c()
for (cov in c(1)){
  for (sgn in c(-1)){
    for (wg in c(0, 1)){
      for (k in c(1.5)){
        for (drop in c(0, 0.3, 0.5)){
          conda <- rbind(conda, c(cov, sgn, wg, k, drop))
        }
      }
    }
  }
}

conda <- as.data.frame(conda)
names(conda) <- c("cov", "sgn", "wg", "k", "drop")


gc()
t<-proc.time()

cluster = makeCluster(40)
clusterSetRNGStream(cluster, 12345)
registerDoParallel(cluster)
simr <- foreach (cov = conda$cov, sgn = conda$sgn, wg= conda$wg, 
                 drop = conda$drop, .combine = rbind, 
                 .packages = c("nlme", "Matrix", "MASS","bcmixed")) %dopar%  
  {
    simset_wg(cov, sgn, wg, 1.5, 50, drop)
  }
stopCluster(cluster)
simr



simr <- as.data.frame(simr)


names(simr) <- c('cov', 'sgn', 'wg', 'k', 'drop', 'npg', 
                       'th1_1', 'th1_2', 'th1_3', 'th2_3', 'sd1', 'sd2', 'sd3', 
                       'Ilgt', 'xi2_1', 'xi2_2', 'xi2_3', 'drrate1_2', 
                       'drrate1_3', 'drratte2_2', 'drrate2_3')
proc.time()-t
write.csv(simr, 'Simulation_bcmixed_WG_parm.csv', row.names = F)
