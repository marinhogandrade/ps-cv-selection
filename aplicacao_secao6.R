# =====================================================================
# Real data application of Section 6
#   "Discrete Data Modeling with Optimal Variable Selection via Cross-Validation"
#
# Data: hh.csv (25,000 hospital encounters of patients with diabetes)
# Response: n_inpatient (inpatient visits in the year preceding the encounter)
# Candidate predictors (7 blocks):
#   X1 age (6 levels, ref. [40-50)), X2 log(1 + n_emergency),
#   X3 log(1 + n_outpatient),
#   X4 n_medications, X5 n_lab_procedures, X6 n_procedures,
#   X7 diag_1 (7 groups, ref. Circulatory; the 4 "Missing" values -> "Other")
#
# Steps (Sections 3 and 4 of the paper):
#   1. descriptive statistics of Y
#   2. Algorithm 1 for P, NB and GP, k = 5 and k = 10, q* = 5 (120 subsets),
#      same folds for all subsets and distributions
#   3. corrected Friedman test (Section 4.1) and, if needed, pairwise tests
#   4. BIC of the selected models (Section 4.2); AIC/BIC benchmarks (Section 4.4)
#   5. observed vs expected frequencies (Table app_freq)
#   6. estimates of the final model with standard errors (Table app_nb)
#   7. Figure: fold errors and ranks
#
# Estimation as in Section 3: Poisson by glm.fit; NB and GP by nloptr
# (NLOPT_LD_LBFGS, analytical gradient, gamma = log(phi) in [-10, 10],
# starting values from the Poisson fit).
# =====================================================================

rm(list = ls(all = TRUE))
setwd("C:/Marinho/Trabalho/Orientacoes/Domingos-Katiane/Programa")
library(nloptr)

SEED    <- 2026
Q_MAX   <- 5
K_FOLDS <- c(5, 10)
DISTS   <- c("P", "NB", "GP")
LPHI_LB <- -10; LPHI_UB <- 10

## ---------------- data ----------------------------------------------
dat <- read.csv("hh.csv", stringsAsFactors = FALSE)
dat$diag_1[dat$diag_1 == "Missing"] <- "Other"          # 4 encounters
# unordered factors -> treatment contrasts (an *ordered* factor would give
# polynomial contrasts, whose coefficients are not level-vs-reference effects)
dat$age    <- factor(dat$age, levels = c("[40-50)", "[50-60)", "[60-70)",
                                         "[70-80)", "[80-90)", "[90-100)"))
dat$diag_1 <- factor(dat$diag_1, levels = c("Circulatory", "Diabetes", "Digestive",
                                            "Injury", "Musculoskeletal", "Other",
                                            "Respiratory"))
y <- dat$n_inpatient
n <- length(y)

# design matrix with blocks 1..7 (block 0 = intercept)
# the prior-visit counts X2 and X3 are very skewed (maxima 64 and 33, most values
# zero) and enter as log(1 + x) to avoid extreme fitted means under the log link
blocks <- list(age = ~ age, n_emergency = ~ log1p(n_emergency),
               n_outpatient = ~ log1p(n_outpatient),
               n_medications = ~ n_medications, n_lab_procedures = ~ n_lab_procedures,
               n_procedures = ~ n_procedures, diag_1 = ~ diag_1)
Xb <- lapply(blocks, function(f) model.matrix(f, dat)[, -1, drop = FALSE])
X  <- cbind("(Intercept)" = 1, do.call(cbind, Xb))
BLOCK <- c(0, rep(seq_along(Xb), sapply(Xb, ncol)))
VARNAME <- paste0("X", 1:7)

## ---------------- 1. descriptive statistics ------------------------
desc <- c(n = n, max = max(y), mean = mean(y), var = var(y),
          pct_zero = 100 * mean(y == 0), var_over_mean = var(y) / mean(y))
print(round(desc, 3))

## ---------------- estimation (Section 3) ----------------------------
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
  list(beta = res$solution[-pp], phi = exp(res$solution[pp]), gamma = res$solution[pp],
       ll = -res$objective, ok = res$status > 0)
}
fit_model <- function(dist, X, y) if (dist == "P") fit_pois(X, y) else fit_disp(dist, X, y)

## ---------------- 2. Algorithm 1 ------------------------------------
SUBSETS <- unlist(lapply(0:Q_MAX, function(k)
  if (k == 0) list(integer(0)) else combn(7, k, simplify = FALSE)), recursive = FALSE)
NM <- length(SUBSETS)                                   # 120
cols_of <- function(s) which(BLOCK %in% c(0, s))

set.seed(SEED)
FOLDS <- lapply(K_FOLDS, function(k) sample(rep_len(1:k, n))); names(FOLDS) <- K_FOLDS

res <- list()
t0 <- Sys.time()
for (dist in DISTS) for (j in seq_len(NM)) {
  cl <- cols_of(SUBSETS[[j]]); Xs <- X[, cl, drop = FALSE]
  ft <- fit_model(dist, Xs, y)
  r <- list(ll = ft$ll, p = length(cl) + (dist != "P"), ok = ft$ok,
            mse = mean((y - exp(Xs %*% ft$beta))^2))
  for (kk in names(FOLDS)) {
    f <- FOLDS[[kk]]
    r[[paste0("E", kk)]] <- sapply(seq_len(as.integer(kk)), function(i) {
      tr <- f != i
      fi <- fit_model(dist, Xs[tr, , drop = FALSE], y[tr])
      mean((y[!tr] - exp(Xs[!tr, , drop = FALSE] %*% fi$beta))^2)
    })
  }
  res[[paste(dist, j)]] <- r
  if (j %% 20 == 0) cat(dist, j, "of", NM, format(Sys.time() - t0), "\n")
}
saveRDS(list(res = res, SUBSETS = SUBSETS, FOLDS = FOLDS), "app_results.rds")

get <- function(dist, j, what) res[[paste(dist, j)]][[what]]
sub_name <- function(s) if (length(s)) paste(VARNAME[s], collapse = " + ") else "(none)"

## ---------------- 3. selected models, BIC and Friedman --------------
friedman_corrected <- function(E) {
  k <- nrow(E); B <- ncol(E)
  Rk <- t(apply(E, 1, rank)); Rs <- colSums(Rk)
  Q <- 12 / (k * B * (B + 1)) * sum((Rs - k * (B + 1) / 2)^2)
  ck <- (2 * k - 1) / k; ke <- k / ck; Qc <- Q / ck
  Fc <- if (ke * (B - 1) - Qc > 0) (ke - 1) * Qc / (ke * (B - 1) - Qc) else Inf
  c(Q = Q, p_unc = pchisq(Q, B - 1, lower.tail = FALSE),
    Fc = Fc, p_cor = pf(Fc, B - 1, (B - 1) * (ke - 1), lower.tail = FALSE))
}
pairwise_corrected <- function(E) {
  k <- nrow(E); B <- ncol(E); Rs <- colSums(t(apply(E, 1, rank)))
  ck <- (2 * k - 1) / k; pr <- combn(B, 2)
  z <- apply(pr, 2, function(ab) (Rs[ab[1]] - Rs[ab[2]]) / sqrt(ck * k * B * (B + 1) / 6))
  tt <- apply(pr, 2, function(ab) { dd <- E[, ab[1]] - E[, ab[2]]
    mean(dd) / sqrt((1 / k + 1 / (k - 1)) * var(dd)) })
  data.frame(pair = paste(colnames(E)[pr[1, ]], colnames(E)[pr[2, ]], sep = "-"),
             p_z_holm = p.adjust(2 * pnorm(-abs(z)), "holm"),
             p_t_holm = p.adjust(2 * pt(-abs(tt), k - 1), "holm"))
}

SEL <- list()
for (kk in names(FOLDS)) {
  tab <- NULL; E <- NULL
  for (dist in DISTS) {
    cv <- sapply(seq_len(NM), function(j) mean(get(dist, j, paste0("E", kk))))
    j <- which.min(cv)
    bic <- -2 * get(dist, j, "ll") + log(n) * get(dist, j, "p")
    tab <- rbind(tab, data.frame(k = kk, model = dist, subset = sub_name(SUBSETS[[j]]),
                                 sqrt_eps_opt = sqrt(cv[j]), sqrt_mse = sqrt(get(dist, j, "mse")),
                                 BIC = bic, j = j))
    E <- cbind(E, get(dist, j, paste0("E", kk)))
  }
  colnames(E) <- DISTS
  fr <- friedman_corrected(E)
  cat("\n==== k =", kk, "====\n"); print(tab[, -7], digits = 6); print(round(fr, 4))
  if (fr["p_cor"] < 0.05) print(pairwise_corrected(E))
  SEL[[kk]] <- list(tab = tab, E = E, fr = fr)
}

## ---------------- 4. AIC / BIC benchmarks (full data) --------------
cat("\n==== Full-data best subset (Section 4.4) ====\n")
for (dist in DISTS) {
  ll <- sapply(seq_len(NM), function(j) get(dist, j, "ll"))
  p  <- sapply(seq_len(NM), function(j) get(dist, j, "p"))
  cat(dist, " AIC:", sub_name(SUBSETS[[which.min(-2 * ll + 2 * p)]]),
      "  BIC:", sub_name(SUBSETS[[which.min(-2 * ll + log(n) * p)]]), "\n")
}
cat("Convergence rate (full-data fits):",
    mean(sapply(res, function(r) r$ok)), "\n")


## ---------------- 4b. checks cited in the text ----------------------
# CV error of the Poisson model with the subset selected for NB (k = 10)
j_nb <- SEL[["10"]]$tab$j[2]
cat("P with NB subset (", sub_name(SUBSETS[[j_nb]]), "): sqrt CV10 =",
    round(sqrt(mean(get("P", j_nb, "E10"))), 4), "\n")
# CV error and BIC of the full-data BIC subset for each distribution
for (dist in DISTS) {
  ll <- sapply(seq_len(NM), function(j) get(dist, j, "ll"))
  p  <- sapply(seq_len(NM), function(j) get(dist, j, "p"))
  jb <- which.min(-2 * ll + log(n) * p)
  cat(dist, "BIC subset", sub_name(SUBSETS[[jb]]), ": BIC =", round(-2 * ll[jb] + log(n) * p[jb], 2),
      " sqrt CV10 =", round(sqrt(mean(get(dist, jb, "E10"))), 4), "\n")
}

## ---------------- 5. observed vs expected frequencies ---------------
dgpois <- function(y, mu, phi)
  exp(y * log(mu / (1 + phi * mu)) + (y - 1) * log1p(phi * y) - lgamma(y + 1) -
        mu * (1 + phi * y) / (1 + phi * mu))
jsel <- SEL[["10"]]$tab$j; names(jsel) <- DISTS
final <- lapply(DISTS, function(dist) {
  cl <- cols_of(SUBSETS[[jsel[dist]]]); Xs <- X[, cl, drop = FALSE]
  ft <- fit_model(dist, Xs, y); ft$X <- Xs; ft$mu <- drop(exp(Xs %*% ft$beta)); ft
}); names(final) <- DISTS
grid <- 0:3
expf <- function(pmf) { e <- sapply(grid, function(v) sum(pmf(v))); c(e, n - sum(e)) }
freq <- data.frame(y = c(grid, ">=4"),
                   Observed = c(sapply(grid, function(v) sum(y == v)), sum(y > 3)),
                   P  = expf(function(v) dpois(v, final$P$mu)),
                   NB = expf(function(v) dnbinom(v, size = final$NB$phi, mu = final$NB$mu)),
                   GP = expf(function(v) dgpois(v, final$GP$mu, final$GP$phi)))
print(freq, digits = 6)
cat("Largest fitted mean:", round(sapply(final, function(f) max(f$mu)), 2), "\n")
cat("Pearson discrepancy:", round(sapply(DISTS, function(d) sum((freq$Observed - freq[[d]])^2 / freq[[d]])), 1), "\n")

## ---------------- 6. estimates of the final model -------------------
# final distribution = smallest BIC among the selected models
bic_sel <- SEL[["10"]]$tab$BIC; names(bic_sel) <- DISTS
best <- names(which.min(bic_sel)); cat("\nFinal model by BIC:", best, "\n")
est_table <- function(dist) {
  ft <- final[[dist]]; Xs <- ft$X; lgy1 <- lgamma(y + 1)
  if (dist == "P") {
    V <- solve(crossprod(Xs, Xs * ft$mu)); est <- ft$beta; se <- sqrt(diag(V)); nm <- colnames(Xs)
  } else {
    f <- if (dist == "NB") nll_nb else nll_gp
    par <- c(ft$beta, ft$gamma)
    H <- optimHess(par, fn = function(p) f(p, Xs, y, lgy1)$objective,
                   gr = function(p) f(p, Xs, y, lgy1)$gradient)
    V <- solve(H); pp <- length(par)
    est <- c(ft$beta, ft$phi)
    se  <- c(sqrt(diag(V))[-pp], ft$phi * sqrt(V[pp, pp]))     # delta method for phi
    nm  <- c(colnames(Xs), "phi")
  }
  z <- est / se
  data.frame(term = nm, estimate = est, se = se, exp_est = exp(est),
             p_value = ifelse(nm == "phi", NA, 2 * pnorm(-abs(z))))
}
for (dist in DISTS) { cat("\n----", dist, "----\n"); print(est_table(dist), digits = 4) }

## ---------------- 7. figure: fold errors and ranks ------------------
setEPS()
for (kk in names(FOLDS)) {
  E <- SEL[[kk]]$E
  postscript(sprintf("boxplot_realdata_k%s.eps", kk), width = 5, height = 4.5)
  boxplot(E, names = DISTS, ylab = expression(paste("Error (", E[iota], ")")))
  dev.off()
  postscript(sprintf("matplot_realdata_k%s.eps", kk), width = 5, height = 4.5)
  matplot(1:nrow(E), t(apply(E, 1, rank)), type = "b", pch = 1:3, lty = 1:3, col = 1,
          xlab = "Fold", ylab = "Rank", ylim = c(1, 3))
  legend("topright", legend = DISTS, pch = 1:3, lty = 1:3, bty = "n")
  dev.off()
}

source("aplicacao_secao6.R", echo = TRUE)
source("salvar_resultados_secao6.R")
