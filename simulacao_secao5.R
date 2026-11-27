# =====================================================================
# Simulation study of Section 5
#   "Discrete Data Modeling with Optimal Variable Selection via Cross-Validation"
#
# Design (as described in Section 5, now with independent replications):
#   - 8 covariates X1..X8; true model uses X1, X5, X6, X8 (categorical X7, X8
#     enter as blocks of indicators); n = 2000 training observations and an
#     independent test sample of n_test = 5000 observations per replicate
#   - true distributions: P, NB (phi = 5), GP (phi = 0.7)
#   - fitted distributions: P, NB, GP
#   - best subset selection with q* = 4 and q* = 5 (all 219 subsets of size
#     <= 5 are fitted once; q* = 4 uses the 163 subsets of size <= 4)
#   - selection procedures (Sections 4 and 4.4): CV with k = 5 and k = 10,
#     training error, AIC and BIC on the full training sample
#   - distribution chosen by BIC (Section 4.2); test error (Section 4.3);
#     corrected Friedman test (Section 4.1)
#
# Estimation:
#   - Poisson: glm.fit (IRLS)
#   - NB and GP: nloptr (NLOPT_LD_LBFGS) with analytical gradient, parameter
#     vector (beta, log phi), log phi bounded in [-10, 10]; starting values for
#     beta = Poisson glm estimates, log phi from moment estimates
#
# Output: one .rds file per replicate in OUT_DIR, then summary tables.
# =====================================================================


rm(list=ls(all=TRUE))
closeAllConnections()          # closes clusters left open by previous runs
setwd("C:\\Marinho\\Trabalho\\Orientacoes\\Domingos-Katiane\\Programa")

library(nloptr)
library(parallel)

## ---------------- settings ------------------------------------------
R_REP    <- 100                 # replications per true distribution
R_TARGET <- 100                 # used only for the time estimate
FRESH    <- TRUE                # TRUE: delete OUT_DIR and start from scratch
                                # FALSE: resume (keep replicates already saved)
N_TRAIN  <- 2000
N_TEST   <- 5000
K_FOLDS  <- c(5, 10)
Q_MAX    <- 5                   # q* = 5 (q* = 4 obtained by restriction)
PHI_TRUE <- c(NB = 5, GP = 0.7)
DISTS    <- c("P", "NB", "GP")
TRUE_BLOCKS <- c(1, 5, 6, 8)
LPHI_LB  <- -10; LPHI_UB <- 10
N_CORES  <- max(1, detectCores() - 1)
OUT_DIR  <- "sim_out"
if (FRESH) unlink(OUT_DIR, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE)
if (FRESH) unlink("sim_summary.rds")

## ---------------- data generation -----------------------------------
gen_x <- function(n) {
  data.frame(x1 = rnorm(n, 45, 15), x2 = rnorm(n, 2000, 800),
             x3 = rbinom(n, 1, .25), x4 = rbinom(n, 1, .50), x5 = rbinom(n, 1, .75),
             x6 = rpois(n, 1.5), x7 = rbinom(n, 2, .4), x8 = rbinom(n, 3, .7))
}

# design matrix: intercept + blocks 1..8; continuous columns are scaled by
# fixed constants (numerical stability only; fitted means are unchanged)
design <- function(X) {
  M <- cbind(1, (X$x1 - 45) / 15, (X$x2 - 2000) / 800, X$x3, X$x4, X$x5,
             (X$x6 - 1.5) / sqrt(1.5),
             X$x7 == 1, X$x7 == 2, X$x8 == 1, X$x8 == 2, X$x8 == 3)
  storage.mode(M) <- "double"
  attr(M, "block") <- c(0, 1, 2, 3, 4, 5, 6, 7, 7, 8, 8, 8)
  M
}

true_mu <- function(X) {
  exp(0.300 + 0.015 * X$x1 + 0.350 * X$x5 + 0.080 * X$x6 +
        0.450 * (X$x8 == 1) + 0.750 * (X$x8 == 2) + 1.050 * (X$x8 == 3))
}

# GP with mean mu and Var = mu (1 + phi mu)^2, generated as the total progeny
# of a Poisson branching process (theta = mu/(1+phi mu), lambda = phi mu/(1+phi mu))
rgpois <- function(mu, phi) {
  theta <- mu / (1 + phi * mu); lam <- phi * mu / (1 + phi * mu)
  cur <- rpois(length(mu), theta); tot <- cur
  while (any(cur > 0)) { cur <- rpois(length(mu), lam * cur); tot <- tot + cur }
  tot
}

gen_y <- function(dist, mu) {
  switch(dist,
         P  = rpois(length(mu), mu),
         NB = rnbinom(length(mu), size = PHI_TRUE["NB"], mu = mu),
         GP = rgpois(mu, PHI_TRUE["GP"]))
}

## ---------------- estimation ----------------------------------------
fit_pois <- function(X, y) {
  g <- glm.fit(X, y, family = poisson())
  b <- g$coefficients; mu <- drop(exp(X %*% b))
  list(beta = b, phi = NA, ll = sum(dpois(y, mu, log = TRUE)), ok = g$converged)
}

# negative log-likelihood and gradient, parameter = c(beta, log(phi))
nll_nb <- function(par, X, y, lgy1) {
  p <- length(par); b <- par[-p]; phi <- exp(par[p])
  mu <- drop(exp(X %*% b))
  ll <- lgamma(y + phi) - lgamma(phi) - lgy1 + phi * log(phi / (phi + mu)) +
    y * log(mu / (phi + mu))
  s_eta <- (y - mu) * phi / (phi + mu)
  s_phi <- digamma(y + phi) - digamma(phi) + log(phi / (phi + mu)) + (mu - y) / (phi + mu)
  list(objective = -sum(ll), gradient = -c(crossprod(X, s_eta), sum(s_phi) * phi))
}

nll_gp <- function(par, X, y, lgy1) {
  p <- length(par); b <- par[-p]; phi <- exp(par[p])
  mu <- drop(exp(X %*% b)); a <- 1 + phi * mu; cc <- 1 + phi * y
  ll <- y * log(mu) - y * log(a) + (y - 1) * log(cc) - lgy1 - mu * cc / a
  s_eta <- (y - mu) / a^2
  s_phi <- -y * mu / a + y * (y - 1) / cc - mu * (y - mu) / a^2
  list(objective = -sum(ll), gradient = -c(crossprod(X, s_eta), sum(s_phi) * phi))
}

NLOPT_OPTS <- list(algorithm = "NLOPT_LD_LBFGS", xtol_rel = 1e-8, maxeval = 1000)

fit_disp <- function(dist, X, y) {
  p0 <- fit_pois(X, y); mu <- drop(exp(X %*% p0$beta))
  if (dist == "NB") {
    den <- sum((y - mu)^2 - mu)
    phi0 <- if (den > 0) sum(mu^2) / den else 100
    f <- nll_nb
  } else {
    r <- mean((y - mu)^2 / mu)
    phi0 <- if (r > 1) (sqrt(r) - 1) / mean(mu) else 1e-3
    f <- nll_gp
  }
  lphi0 <- min(max(log(min(max(phi0, 1e-4), 1e4)), LPHI_LB), LPHI_UB)
  lgy1 <- lgamma(y + 1)
  pp <- ncol(X) + 1
  lb <- c(rep(-Inf, pp - 1), LPHI_LB); ub <- c(rep(Inf, pp - 1), LPHI_UB)
  run <- function(x0) nloptr(x0, eval_f = f, lb = lb, ub = ub, opts = NLOPT_OPTS,
                             X = X, y = y, lgy1 = lgy1)
  res <- run(c(p0$beta, lphi0))
  if (res$status < 0 || !is.finite(res$objective))       # restart once
    res <- run(c(p0$beta, 0))
  list(beta = res$solution[-pp], phi = exp(res$solution[pp]),
       ll = -res$objective, ok = res$status > 0)
}

fit_model <- function(dist, X, y) if (dist == "P") fit_pois(X, y) else fit_disp(dist, X, y)

## ---------------- subsets and folds ---------------------------------
SUBSETS <- unlist(lapply(0:Q_MAX, function(k)
  if (k == 0) list(integer(0)) else combn(8, k, simplify = FALSE)), recursive = FALSE)
SIZE <- lengths(SUBSETS)                               # 219 subsets
cols_of <- function(s, block) which(block %in% c(0, s))
make_folds <- function(n, k) sample(rep_len(1:k, n))

## ---------------- one replicate -------------------------------------
run_replicate <- function(true_dist, seed) {
  set.seed(seed)
  Xtr <- gen_x(N_TRAIN); Xte <- gen_x(N_TEST)
  ytr <- gen_y(true_dist, true_mu(Xtr)); yte <- gen_y(true_dist, true_mu(Xte))
  Dtr <- design(Xtr); Dte <- design(Xte); block <- attr(Dtr, "block")
  fl <- lapply(K_FOLDS, function(k) make_folds(N_TRAIN, k)); names(fl) <- K_FOLDS
  nS <- length(SUBSETS)
  A <- function() matrix(NA_real_, 3, nS, dimnames = list(DISTS, NULL))
  out <- list(mse = A(), ll = A(), p = A(), test = A(), conv = A(),
              cv = lapply(K_FOLDS, function(k) array(NA_real_, c(3, nS, k))))
  names(out$cv) <- K_FOLDS
  for (d in seq_along(DISTS)) for (j in seq_len(nS)) {
    dist <- DISTS[d]; cl <- cols_of(SUBSETS[[j]], block); X <- Dtr[, cl, drop = FALSE]
    ft <- fit_model(dist, X, ytr)
    out$mse[d, j]  <- mean((ytr - exp(X %*% ft$beta))^2)
    out$ll[d, j]   <- ft$ll
    out$p[d, j]    <- length(cl) + (dist != "P")
    out$conv[d, j] <- ft$ok
    out$test[d, j] <- mean((yte - exp(Dte[, cl, drop = FALSE] %*% ft$beta))^2)
    for (kk in names(fl)) {
      f <- fl[[kk]]
      for (i in seq_len(as.integer(kk))) {
        tr <- f != i
        fi <- fit_model(dist, X[tr, , drop = FALSE], ytr[tr])
        out$cv[[kk]][d, j, i] <- mean((ytr[!tr] - exp(X[!tr, , drop = FALSE] %*% fi$beta))^2)
      }
    }
  }
  out$true <- true_dist; out$seed <- seed
  out
}



## ---------------- run (parallel, resumable) -------------------------
jobs <- expand.grid(rep = seq_len(R_REP), true = DISTS, stringsAsFactors = FALSE)
seed_base <- c(P = 10000, NB = 20000, GP = 30000)

run_job <- function(i) {
  fn <- file.path(OUT_DIR, sprintf("%s_%03d.rds", jobs$true[i], jobs$rep[i]))
  if (!file.exists(fn))
    saveRDS(run_replicate(jobs$true[i], seed_base[jobs$true[i]] + jobs$rep[i]), fn)
  NULL
}

n_todo <- sum(!file.exists(file.path(OUT_DIR,
                sprintf("%s_%03d.rds", jobs$true, jobs$rep))))
cat("Replicates to run:", n_todo, "of", nrow(jobs), "using", N_CORES, "cores\n")

tempo <- system.time({
  if (.Platform$OS.type == "windows") {
    cl <- makeCluster(N_CORES)                     # PSOCK cluster (works on Windows)
    tryCatch({
      clusterEvalQ(cl, library(nloptr))
      clusterExport(cl, ls(globalenv()))           # functions, settings, jobs, ...
      invisible(parLapplyLB(cl, seq_len(nrow(jobs)), run_job))
    }, finally = stopCluster(cl))                  # closes the cluster even on error
  } else {
    invisible(mclapply(seq_len(nrow(jobs)), run_job, mc.cores = N_CORES))
  }
})

if (n_todo > 0) {
  # each job = one replicate of one true distribution; jobs run N_CORES at a time
  t_job   <- tempo["elapsed"] / ceiling(n_todo / N_CORES)          # seconds per round
  t_total <- t_job * ceiling(length(DISTS) * R_TARGET / N_CORES)
  cat(sprintf("Elapsed: %.1f min for %d jobs (%.1f min per job, approx.)\n",
              tempo["elapsed"] / 60, n_todo, t_job / 60))
  cat(sprintf("Estimated time for %d replicates per distribution: %.1f hours (%d cores)\n",
              R_TARGET, t_total / 3600, N_CORES))
} else cat("Nothing to run: all replicates were already in", OUT_DIR, "\n")

## ---------------- corrected Friedman test (Section 4.1) -------------
friedman_corrected <- function(E) {
  k <- nrow(E); B <- ncol(E)
  Rk <- t(apply(E, 1, rank)); Rs <- colSums(Rk)
  Q <- 12 / (k * B * (B + 1)) * sum((Rs - k * (B + 1) / 2)^2)
  ck <- (2 * k - 1) / k; ke <- k / ck; Qc <- Q / ck
  Fc <- if (ke * (B - 1) - Qc > 0) (ke - 1) * Qc / (ke * (B - 1) - Qc) else Inf
  c(p_unc = pchisq(Q, B - 1, lower.tail = FALSE),
    p_cor = pf(Fc, B - 1, (B - 1) * (ke - 1), lower.tail = FALSE))
}

## ---------------- summaries -----------------------------------------
JTRUE <- which(sapply(SUBSETS, function(s) identical(as.integer(s), as.integer(TRUE_BLOCKS))))
PROCS <- c("CV5", "CV10", "TRAIN", "AIC", "BIC")

summarise_rep <- function(o) {
  ti <- match(o$true, DISTS)
  aic <- -2 * o$ll + 2 * o$p; bic <- -2 * o$ll + log(N_TRAIN) * o$p
  cvm <- lapply(o$cv, function(a) apply(a, c(1, 2), mean))
  crit <- list(CV5 = cvm[["5"]], CV10 = cvm[["10"]], TRAIN = o$mse, AIC = aic, BIC = bic)
  res <- list()
  for (qs in c(4, 5)) for (proc in PROCS) {
    C <- crit[[proc]]; C[, SIZE > qs] <- Inf
    sel <- apply(C, 1, which.min)
    dstar <- which.min(sapply(1:3, function(d) bic[d, sel[d]]))
    st <- SUBSETS[[sel[ti]]]
    r <- data.frame(true = o$true, q = qs, proc = proc,
                    exact = setequal(st, TRUE_BLOCKS) * 1,
                    fp = length(setdiff(st, TRUE_BLOCKS)),
                    fn = length(setdiff(TRUE_BLOCKS, st)), size = length(st),
                    dist_ok = (dstar == ti) * 1,
                    test_final = sqrt(o$test[dstar, sel[dstar]]),
                    test_true = sqrt(o$test[ti, sel[ti]]),
                    oracle = sqrt(o$test[ti, JTRUE]),
                    opt = NA, fu = NA, fc = NA)
    for (v in 1:8) r[[paste0("X", v)]] <- (v %in% st) * 1
    if (proc %in% c("CV5", "CV10")) {
      kk <- if (proc == "CV5") "5" else "10"
      r$opt <- sqrt(crit[[proc]][ti, sel[ti]]) - sqrt(o$test[ti, sel[ti]])
      E <- sapply(1:3, function(d) o$cv[[kk]][d, sel[d], ])
      fr <- friedman_corrected(E); r$fu <- fr["p_unc"]; r$fc <- fr["p_cor"]
    }
    res[[length(res) + 1]] <- r
  }
  do.call(rbind, res)
}

files <- list.files(OUT_DIR, pattern = "\\.rds$", full.names = TRUE)
SUM <- do.call(rbind, lapply(files, function(f) summarise_rep(readRDS(f))))
conv_rate <- mean(unlist(lapply(files, function(f) readRDS(f)$conv)))
cat("Replicates:", length(files), " convergence rate:", round(conv_rate, 4), "\n")

# Table: variable selection under the true distribution
print(aggregate(cbind(exact, fp, fn, size) ~ true + q + proc, data = SUM, FUN = mean))
# Table: inclusion frequencies (q* = 5)
print(aggregate(cbind(X1, X2, X3, X4, X5, X6, X7, X8) ~ true + proc,
                data = subset(SUM, q == 5), FUN = mean))
# Table: complete procedure
print(aggregate(cbind(dist_ok, test_final, oracle) ~ true + q + proc, data = SUM, FUN = mean))
print(aggregate(opt ~ true + q + proc, data = subset(SUM, !is.na(opt)), FUN = mean))
# Table: Friedman rejection rates at 5%
SUM$ru <- SUM$fu < 0.05; SUM$rc <- SUM$fc < 0.05
print(aggregate(cbind(ru, rc) ~ true + q + proc, data = subset(SUM, !is.na(fu)), FUN = mean))
saveRDS(SUM, "sim_summary.rds")


setwd("C:/Marinho/Trabalho/Orientacoes/Domingos-Katiane/Programa")

# 1) resultados por réplica (todas as tabelas da Seção 5)
SUM <- readRDS("sim_summary.rds")
write.csv(SUM, "sim_summary.csv", row.names = FALSE)

# 2) números citados no texto: convergência e diferenças de BIC
files <- list.files("sim_out", pattern = "\\.rds$", full.names = TRUE)
SUBSETS <- unlist(lapply(0:5, function(k)
  if (k == 0) list(integer(0)) else combn(8, k, simplify = FALSE)), recursive = FALSE)
JTRUE <- which(sapply(SUBSETS, function(s) identical(as.integer(s), c(1L, 5L, 6L, 8L))))
extra <- do.call(rbind, lapply(files, function(f) {
  o <- readRDS(f)
  bic <- -2 * o$ll + log(2000) * o$p
  data.frame(true = o$true, seed = o$seed, conv = mean(o$conv),
             bic_P = bic[1, JTRUE], bic_NB = bic[2, JTRUE], bic_GP = bic[3, JTRUE])
}))
write.csv(extra, "sim_extra.csv", row.names = FALSE)
cat("Arquivos:", length(files), "\n")
