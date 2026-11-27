f <- list.files("C:/Marinho/Trabalho/Orientacoes/Domingos-Katiane/Programa/sim2_out",
                pattern = "\\.rds$")
table(n = substr(f, 2, 6), dist = sub("^n\\d+_([A-Z]+)_.*$", "\\1", f))

#summarise_rep

setwd("C:/Marinho/Trabalho/Orientacoes/Domingos-Katiane/Programa")
N_DONE <- c(500, 2000)          # scenarios already complete

# definitions needed by summarise_rep (same as in the main script)
DISTS <- c("P", "NB", "GP"); TRUE_BLOCKS <- c(1, 5, 6, 8)
SUBSETS <- unlist(lapply(0:5, function(k)
  if (k == 0) list(integer(0)) else combn(8, k, simplify = FALSE)), recursive = FALSE)
SIZE  <- lengths(SUBSETS)
JTRUE <- which(sapply(SUBSETS, function(s) identical(as.integer(s), as.integer(TRUE_BLOCKS))))
PROCS <- c("CV5", "CV10", "TRAIN", "AIC", "BIC")

friedman_corrected <- function(E) {
  k <- nrow(E); B <- ncol(E)
  Rk <- t(apply(E, 1, rank)); Rs <- colSums(Rk)
  Q <- 12 / (k * B * (B + 1)) * sum((Rs - k * (B + 1) / 2)^2)
  ck <- (2 * k - 1) / k; ke <- k / ck; Qc <- Q / ck
  Fc <- if (ke * (B - 1) - Qc > 0) (ke - 1) * Qc / (ke * (B - 1) - Qc) else Inf
  c(p_unc = pchisq(Q, B - 1, lower.tail = FALSE),
    p_cor = pf(Fc, B - 1, (B - 1) * (ke - 1), lower.tail = FALSE))
}

summarise_rep <- function(o) {
  ti <- match(o$true, DISTS)
  aic <- -2 * o$ll + 2 * o$p; bic <- -2 * o$ll + log(o$n) * o$p
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

# only the files of the completed scenarios
pat   <- paste0("^n(", paste(sprintf("%05d", N_DONE), collapse = "|"), ")_.*\\.rds$")
files <- list.files("sim2_out", pattern = pat, full.names = TRUE)
cat("Files read:", length(files), "(expected", 300 * length(N_DONE), ")\n")

SUM <- do.call(rbind, lapply(files, function(f) summarise_rep(readRDS(f))))
SUM$ru <- SUM$fu < 0.05; SUM$rc <- SUM$fc < 0.05
write.csv(SUM, "sim2_summary_parcial.csv", row.names = FALSE)

extra <- do.call(rbind, lapply(files, function(f) {
  o <- readRDS(f); bic <- -2 * o$ll + log(o$n) * o$p
  data.frame(n = o$n, rho = o$rho, true = o$true, seed = o$seed,
             conv = mean(o$conv), cor_x1x2 = o$cor_x1x2,
             bic_P = bic[1, JTRUE], bic_NB = bic[2, JTRUE], bic_GP = bic[3, JTRUE])
}))
write.csv(extra, "sim2_extra_parcial.csv", row.names = FALSE)
print(table(n = extra$n, dist = extra$true))