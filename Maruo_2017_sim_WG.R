# Simulation program for Maruo et al. (2017). Stat. Med. 36, 15, 2420-2434
# Setting: Weibull and gamma distributions


#install.packages("doParallel")
#install.packages("nlme")
#install.packages("Matrix")
#install.packages("contrast")
#install.packages("e1071")
library(doParallel)

# Data generation
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
  
  bl <- c((numeric(tn) + 1) %x% y00[1:npg, 1], (numeric(tn) + 1) %x% 
            y00[(npg + 1):(2 * npg), 1])
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

# Calculate LOCF variable from data
LOCF_dat <- function(dat){
  N <- nrow(dat) / 3
  dat1 <- dat
  j <- 1
  for (i in 1:N){
    ns <- (i - 1) * 3 + 1
    ndat <- dat1$y[ns:(ns + 2)]
    miss <- sum(is.na(ndat))
    if (miss == 1) ndat[3] <- ndat[2]
    if (miss == 2) ndat[2:3] <- ndat[1]
    dat1$y[ns:(ns + 2)] <- ndat
  }
  return(dat1)
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
  form <- paste("y~ as.factor(group) + as.factor(time) + as.factor(group):as.factor(time)", 
                covform)
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
  lsmeans_diff <- data.frame(group1 = comb[, 2], group0 = comb[, 1], 
                             estimate = delta, SE = SEd, lower_CL = lower, 
                             upper_CL = upper, df = df, 
                             tvalue = tvalue, pvalue = pvalue)
  result <- list(lsmeans = lsmeans, lsmeans_diff = lsmeans_diff)
  return(result)
}

# bcmixed part of simulation
MMRM_bc_0 <- function(data, cov, deltat){
  if (cov == 1){
    res00 <- try(bcmmrm(outcome = y, group = group, data = data, time = time, 
                        id = id, covv = "lbl", cfactor = 0))
    res10 <- try(bcmmrm(outcome = y, group = group, data = data, time = time, 
                        id = id, covv = "lbl", cfactor = 0, 
                        structure = "CS"))
  } else {
    res00 <- try(bcmmrm(outcome = y, group = group, data = data, time = time, 
                        id = id))
    res10 <- try(bcmmrm(outcome = y, group = group, data = data, time = time, 
                        id = id, structure = "CS"))
  }
  if (class(res00) == "try-error"|res00 == "Not converged") {
    res01 <- c(numeric(14) * NA, 0)
  } else {
    cn0 <- 100 * (res00$meddif.mod[[3]]$lower.CL < deltat &
                    deltat < res00$meddif.mod[[3]]$upper.CL)
    ca0 <- 100 * (res00$meddif.mod.adj[[3]]$lower.CL < deltat &
                    deltat < res00$meddif.mod.adj[[3]]$upper.CL)
    crn0 <- 100 * (res00$meddif.rob[[3]]$lower.CL < deltat &
                     deltat < res00$meddif.rob[[3]]$upper.CL)
    cra0 <- 100 * (res00$meddif.rob.adj[[3]]$lower.CL < deltat &
                     deltat < res00$meddif.rob.adj[[3]]$upper.CL)
    sn0 <- 100 * (res00$meddif.mod[[3]]$p.value < 0.05)
    sa0 <- 100 * (res00$meddif.mod.adj[[3]]$p.value < 0.05)
    srn0 <- 100 * (res00$meddif.rob[[3]]$p.value < 0.05)
    sra0 <- 100 * (res00$meddif.rob.adj[[3]]$p.value < 0.05)
    res01 <- c(res00$lambda, as.numeric(res00$meddif.mod[[3]][c(3, 4)]), 
               res00$meddif.rob[[3]]$SE, res00$meddif.mod.adj[[3]]$SE, 
               res00$meddif.rob.adj[[3]]$SE, 
               cn0, crn0, ca0, cra0, sn0, srn0, sa0, sra0, 100)
  }
  if (class(res10) == "try-error"|res10 == "Not converged") {
    res11 <- c(numeric(14) * NA, 0)
  } else {
    cn1 <- 100 * (res10$meddif.mod[[3]]$lower.CL < deltat &
                    deltat < res10$meddif.mod[[3]]$upper.CL)
    ca1 <- 100 * (res10$meddif.mod.adj[[3]]$lower.CL < deltat &
                    deltat < res10$meddif.mod.adj[[3]]$upper.CL)
    crn1 <- 100 * (res10$meddif.rob[[3]]$lower.CL < deltat &
                     deltat < res10$meddif.rob[[3]]$upper.CL)
    cra1 <- 100 * (res10$meddif.rob.adj[[3]]$lower.CL < deltat &
                     deltat < res10$meddif.rob.adj[[3]]$upper.CL)
    sn1 <- 100 * (res10$meddif.mod[[3]]$p.value < 0.05)
    sa1 <- 100 * (res10$meddif.mod.adj[[3]]$p.value < 0.05)
    srn1 <- 100 * (res10$meddif.rob[[3]]$p.value < 0.05)
    sra1 <- 100 * (res10$meddif.rob.adj[[3]]$p.value < 0.05)
    res11 <- c(res10$lambda, as.numeric(res10$meddif.mod[[3]][c(3, 4)]), 
               res10$meddif.rob[[3]]$SE, res10$meddif.mod.adj[[3]]$SE, 
               res10$meddif.rob.adj[[3]]$SE, 
               cn1, crn1, ca1, cra1, sn1, srn1, sa1, sra1, 100)
  }
  return(list(res01, res11))  
}

# Simulation conducting function
bcmixed_sim_wg <- function(cov, sgn, wg, k, drop, npg, th1_1, th1_2, th1_3, 
                           th2_3, sd1, sd2, sd3, Ilgt, xi2_3, h, simn){
  th1 <- c(th1_1, th1_2, th1_3)
  sd <- c(sd1, sd2, sd3)
  deltat <- xi2_3 - (100 - sgn * 20)
  if (h!=2) deltat <- 0
  res00a <- c()
  res01a <- c()
  res1a <- c()
  res2a <- c()
  res3a <- c()
  res4a <- c()
  res5a <- c()
  for (i in 1:simn){
    dat0 <- rand_parm_wg(cov, sgn, wg, k, npg, th1, Ilgt, sd, th2_3, h)
    dat0$lny <- log(dat0$y)
    dat0$lnbl <- log(dat0$bl) 
    datl <- LOCF_dat(dat0)
    res0 <- MMRM_bc_0(dat0, cov, deltat)
    if (cov == 1){
      res10 <- try(bcmmrm(outcome = y, group = group, 
                          data = subset(dat0, time == 3),
                          covv = "lbl", cfactor = 0))
      res20 <- try(bcmmrm(outcome = y, group = group, data = datl,
                          covv = "lbl", cfactor = 0))
    }
    if (cov == 0){
      res10 <- try(bcmmrm(outcome = y, group = group, 
                          data = subset(dat0, time == 3)))
      res20 <- try(bcmmrm(outcome = y, group = group, data = datl))
    }
    if (class(res10) == "try-error") {
      res1 <- NA
    } else {
      res11 <- res10$meddif.mod.adj[[1]]
      res1 <- c(res11$delta, 100 * (res11$lower.CL < deltat & 
                                      deltat < res11$upper.CL),
                100 * (res11$p.value < 0.05))
    }
    if (class(res20) == "try-error") {
      res2 <- NA
    } else {
      res21 <- res20$meddif.mod.adj[[1]]
      res2 <- c(res21$delta, 100 * (res21$lower.CL < deltat & 
                                      deltat < res21$upper.CL),
                100 * (res21$p.value < 0.05))
    }
    if (cov == 1){
      res3 <- try(100 * c(MMRM(outcome = y, group = group, time = time, id = id, 
                               timepoint = 3, data = dat0, covv = "bl", 
                               cfactor = 0)$lsmeans_diff$pvalue < 0.05))
      res4 <- try(100 * c(MMRM(outcome = lny, group = group, time = time, id = id, 
                               timepoint = 3, data = dat0, covv = "lnbl", 
                               cfactor = 0)$lsmeans_diff$pvalue < 0.05))
    } else {
      res3 <- try(100 * c(MMRM(outcome = y, group = group, time = time, id = id, 
                               timepoint = 3, 
                               data = dat0)$lsmeans_diff$pvalue < 0.05))
      res4 <- try(100 * c(MMRM(outcome = lny, group = group, time = time, 
                               id = id, timepoint = 3, 
                               data = dat0)$lsmeans_diff$pvalue < 0.05))
    }
    
    if (class(res2) == "try-error") res2 <- NA
    if (class(res3) == "try-error") res3 <- NA
    if (class(res4) == "try-error") res4 <- NA

    res00a <- rbind(res00a, res0[[1]])
    res01a <- rbind(res01a, res0[[2]])
    res1a <- rbind(res1a, res1)
    res2a <- rbind(res2a, res2)
    res3a <- c(res3a, res3)
    res4a <- c(res4a, res4)
  }
  res00m <- colMeans(res00a, na.rm = T)
  res01m <- colMeans(res01a, na.rm = T)
  return(as.numeric(c(cov, sgn, wg, k, drop, npg, h, deltat, res00m[1:6], 
                      sd(res00a[, 2], na.rm = T), 
                      res00m[7:15], res01m[1:6], sd(res01a[, 2], na.rm = T), 
                      res01m[7:15], colMeans(res1a, na.rm = T), 
                      colMeans(res2a, na.rm = T), 
                      mean(res3a, na.rm = T), mean(res4a, na.rm = T))))
}


gc()
t<-proc.time()


parm <- read.csv('Simulation_bcmixed_WG_parm.csv')
parm <- rbind(cbind(parm, h = 0), cbind(parm, h = 1), cbind(parm, h = 2))
parma <- parm[parm$cov == 1 & parm$sgn == -1,]

cluster = makeCluster(40)
clusterSetRNGStream(cluster, 12345)
registerDoParallel(cluster)
simr <- foreach (wg = parma$wg, drop = parma$drop, th1_1 = parma$th1_1, 
                 th1_2 = parma$th1_2, th1_3 = parma$th1_3, th2_3 = parma$th2_3, 
                 sd1 = parma$sd1, sd2 = parma$sd2, sd3 = parma$sd3, 
                 Ilgt = parma$Ilgt, xi2_3 = parma$xi2_3,  h = parma$h,  
                 .combine = rbind,
                 .packages = c("nlme", "Matrix", "MASS",
                               "e1071", "bcmixed")) %dopar%  
  {
    bcmixed_sim_wg(1, -1, wg, 1.5, drop, 50, th1_1, th1_2, th1_3, 
                   th2_3, sd1, sd2, sd3, Ilgt, xi2_3, h, 10000)
  }
stopCluster(cluster)

proc.time()-t

simr <- as.data.frame(simr)
names(simr) <- c('cov', 'sgn', 'wg', 'k', 'drop', 'npg' ,'H', 'delta_t',
                  "lambda_e_un", "delta_un", "SE_mod_un", "SE_rob_un", 
                  "SE_mod_adj_un", "SE_rob_adj_un", "SD_delta_un", "cp_mod_un", 
                  "cp_rob_un", "cp_mod_adj_un", "cp_rob_adj_un", 
                  "power_mod_un", "power_rob_un", "power_mod_adj_un", 
                  "power_rob_adj_un", "conv_un", "lambda_e_cs", 
                  "delta_cs", "SE_mod_cs", "SE_rob_cs", "SE_mod_adj_cs", 
                  "SE_rob_adj_cs", "SD_delta_cs", "cp_mod_cs", "cp_rob_cs", 
                  "cp_mod_adj_cs", "cp_rob_adj_cs", "power_mod_cs", 
                  "power_rob_cs", "power_mod_adj_cs", "power_rob_adj_cs", 
                  "conv_cs", "delta_cca", "cp_cca", "power_cca", "delta_locf", 
                  "cp_locf", "power_locf", "power_mmrm", "power_mmrm_l")


write.csv(simr,"Simulation_bcmixed_WG.csv",row.names=F)

