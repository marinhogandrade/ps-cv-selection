# =====================================================================
# Simulation study of Section 5 -- STUDY 2
#   "Discrete Data Modeling with Optimal Variable Selection via Cross-Validation"
#
# Same design as Study 1 (simulacao_secao5.R), with two changes:
#   (1) different training sample sizes (N_SET, e.g. 500, 2000 and 10000);
#   (2) correlated covariates: X1 (relevant) and X2 (irrelevant) are bivariate
#       normal with correlation RHO (marginals N(45, 15^2) and N(2000, 800^2),
#       as in Study 1).
# All other covariates, the true model (X1, X5, X6, X8), the coefficients,
# the distributions (P, NB phi = 5, GP phi = 0.7), the selection procedures
# (CV k = 5 and k = 10, training error, AIC, BIC), q* = 4 and 5, and the
# independent test sample (N_TEST) are as in Study 1.
#
# Results are saved in separate files so that Study 1 is not overwritten:
#   folder  sim2_out/            one .rds per replicate: n<N>_<dist>_<rep>.rds
#   files   sim2_summary.rds/.csv, sim2_extra.csv
# =====================================================================

rm(list = ls(all = TRUE))
closeAllConnections()          # closes clusters left open by previous runs
setwd("C:\\Marinho\\Trabalho\\Orientacoes\\Domingos-Katiane\\Programa")

library(nloptr)
library(parallel)

## ---------------- settings ------------------------------------------
N_SET    <- c(500, 2000, 10000) # training sample sizes (scenarios)
RHO      <- 0.7                 # correlation between X1 and X2
R_REP    <- 100                 # replications per scenario and true distribution
N_TEST   <- 5000                # independent test sample (same for all n)
K_FOLDS  <- c(5, 10)
Q_MAX    <- 5                   # q* = 5 (q* = 4 obtained by restriction)
PHI_TRUE <- c(NB = 5, GP = 0.7)
DISTS    <- c("P", "NB", "GP")
TRUE_BLOCKS <- c(1, 5, 6, 8)
LPHI_LB  <- -10; LPHI_UB <- 10
N_CORES  <- max(1, detectCores() - 1)
FRESH    <- FALSE               # TRUE: delete OUT_DIR and start from scratch
                                # FALSE: resume (keep replicates already saved)
OUT_DIR  <- "sim2_out"
SUMMARY_RDS <- "sim2_summary.rds"
SUMMARY_CSV <- "sim2_summary.csv"
EXTRA_CSV   <- "sim2_extra.csv"
if (FRESH) { unlink(OUT_DIR, recursive = TRUE); unlink(c(SUMMARY_RDS, SUMMARY_CSV, EXTRA_CSV)) }
dir.create(OUT_DIR, showWarnings = FALSE)

## ---------------- data generation -----------------------------------
# X1 and X2 bivariate normal with correlation rho; other covariates as in Study 1
gen_x <- function(n, rho) {
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  data.frame(x1 = 45 + 15 * z1, x2 = 2000 + 800 * z2,
             x3 = rbinom(n, 1, .25), x4 = rbinom(n, 1, .50), x5 = rbinom(n, 1, .75),
             x6 = rpois(n, 1.5), x7 = rbinom(n, 2, .4), x8 = rbinom(n, 3, .7))
}

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

## ---------------- estimation (as in Study 1) ------------------------
fit_pois <- function(X, y) {
  g <- glm.fit(X, y, family = poisson())
  b <- g$coefficients; mu <- drop(exp(X %*% b))
  list(beta = b, phi = NA, ll = sum(dpois(y, mu, log = TRUE)), ok = g$converged)
}

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
    den <- sum((y - mu)^2 - mu); phi0 <- if (den > 0) sum(mu^2) / den else 100; f <- nll_nb
  } else {
    r <- mean((y - mu)^2 / mu); phi0 <- if (r > 1) (sqrt(r) - 1) / mean(mu) else 1e-3; f <- nll_gp
  }
  lphi0 <- min(max(log(min(max(phi0, 1e-4), 1e4)), LPHI_LB), LPHI_UB)
  lgy1 <- lgamma(y + 1); pp <- ncol(X) + 1
  lb <- c(rep(-Inf, pp - 1), LPHI_LB); ub <- c(rep(Inf, pp - 1), LPHI_UB)
  run <- function(x0) nloptr(x0, eval_f = f, lb = lb, ub = ub, opts = NLOPT_OPTS,
                             X = X, y = y, lgy1 = lgy1)
  res <- run(c(p0$beta, lphi0))
  if (res$status < 0 || !is.finite(res$objective)) res <- run(c(p0$beta, 0))
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
run_replicate <- function(true_dist, seed, n, rho) {
  set.seed(seed)
  Xtr <- gen_x(n, rho); Xte <- gen_x(N_TEST, rho)
  ytr <- gen_y(true_dist, true_mu(Xtr)); yte <- gen_y(true_dist, true_mu(Xte))
  Dtr <- design(Xtr); Dte <- design(Xte); block <- attr(Dtr, "block")
  fl <- lapply(K_FOLDS, function(k) make_folds(n, k)); names(fl) <- K_FOLDS
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
  out$true <- true_dist; out$seed <- seed; out$n <- n; out$rho <- rho
  out$cor_x1x2 <- cor(Xtr$x1, Xtr$x2)
  out
}

## ---------------- run (parallel, resumable) -------------------------
# jobs ordered by sample size, so that the faster scenarios finish first
jobs <- expand.grid(rep = seq_len(R_REP), true = DISTS, n = sort(N_SET),
                    stringsAsFactors = FALSE)
jobs <- jobs[order(jobs$n, jobs$rep), ]
# seeds differ from Study 1 (10000/20000/30000) and between sample sizes
seed_base <- c(P = 100000, NB = 200000, GP = 300000)
job_file <- function(i) file.path(OUT_DIR, sprintf("n%05d_%s_%03d.rds",
                                                    jobs$n[i], jobs$true[i], jobs$rep[i]))
run_job <- function(i) {
  fn <- job_file(i)
  if (!file.exists(fn))
    saveRDS(run_replicate(jobs$true[i],
                          seed_base[jobs$true[i]] + 1000 * match(jobs$n[i], sort(N_SET)) + jobs$rep[i],
                          jobs$n[i], RHO), fn)
  NULL
}

todo <- !file.exists(sapply(seq_len(nrow(jobs)), job_file))
cat("Jobs to run:", sum(todo), "of", nrow(jobs), "using", N_CORES, "cores\n")
print(table(n = jobs$n[todo], dist = jobs$true[todo]))

tempo <- system.time({
  if (.Platform$OS.type == "windows") {
    cl <- makeCluster(N_CORES)
    tryCatch({
      clusterEvalQ(cl, library(nloptr))
      clusterExport(cl, ls(globalenv()))
      invisible(parLapplyLB(cl, which(todo), run_job))
    }, finally = stopCluster(cl))
  } else {
    invisible(mclapply(which(todo), run_job, mc.cores = N_CORES))
  }
})
cat(sprintf("Elapsed: %.1f hours for %d jobs\n", tempo["elapsed"] / 3600, sum(todo)))

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
  aic <- -2 * o$ll + 2 * o$p; bic <- -2 * o$ll + log(o$n) * o$p   # BIC uses the scenario n
  cvm <- lapply(o$cv, function(a) apply(a, c(1, 2), mean))
  crit <- list(CV5 = cvm[["5"]], CV10 = cvm[["10"]], TRAIN = o$mse, AIC = aic, BIC = bic)
  res <- list()
  for (qs in c(4, 5)) for (proc in PROCS) {
    C <- crit[[proc]]; C[, SIZE > qs] <- Inf
    sel <- apply(C, 1, which.min)
    dstar <- which.min(sapply(1:3, function(d) bic[d, sel[d]]))
    st <- SUBSETS[[sel[ti]]]
    r <- data.frame(n = o$n, rho = o$rho, true = o$true, q = qs, proc = proc,
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
SUM$ru <- SUM$fu < 0.05; SUM$rc <- SUM$fc < 0.05
saveRDS(SUM, SUMMARY_RDS)
write.csv(SUM, SUMMARY_CSV, row.names = FALSE)

# numbers cited in the text: convergence, sample correlation, BIC at the true subset
extra <- do.call(rbind, lapply(files, function(f) {
  o <- readRDS(f)
  bic <- -2 * o$ll + log(o$n) * o$p
  data.frame(n = o$n, rho = o$rho, true = o$true, seed = o$seed,
             conv = mean(o$conv), cor_x1x2 = o$cor_x1x2,
             bic_P = bic[1, JTRUE], bic_NB = bic[2, JTRUE], bic_GP = bic[3, JTRUE])
}))
write.csv(extra, EXTRA_CSV, row.names = FALSE)
cat("Replicate files:", length(files), "\n")
print(table(n = extra$n, dist = extra$true))

## ---------------- quick look at the results -------------------------
print(aggregate(cbind(exact, fp, fn, size) ~ n + true + q + proc, data = SUM, FUN = mean))
print(aggregate(cbind(X1, X2, X3, X4, X5, X6, X7, X8) ~ n + true + proc,
                data = subset(SUM, q == 5), FUN = mean))
print(aggregate(cbind(dist_ok, test_final, oracle) ~ n + true + q + proc, data = SUM, FUN = mean))
print(aggregate(cbind(ru, rc) ~ n + true + q + proc, data = subset(SUM, !is.na(fu)), FUN = mean))



