# =====================================================================
# Save the results of aplicacao_secao6.R to CSV files
# Run in the SAME R session, right after aplicacao_secao6.R
# (it uses the objects dat, y, n, X, SUBSETS, NM, FOLDS, res, SEL, final,
#  freq, best, desc, DISTS and the functions defined in that script).
# Output: folder "resultados_secao6" with one CSV per table.
# =====================================================================

OUT <- "resultados_secao6"
dir.create(OUT, showWarnings = FALSE)
wcsv <- function(x, name) {
  write.csv(x, file.path(OUT, name), row.names = FALSE)
  cat("saved:", file.path(OUT, name), "\n")
}

## 1. Descriptive statistics of Y and of the skewed covariates ---------
wcsv(data.frame(statistic = names(desc), value = as.numeric(desc)),
     "01_descriptive_y.csv")
wcsv(data.frame(
  variable = c("n_emergency", "n_outpatient"),
  max      = c(max(dat$n_emergency), max(dat$n_outpatient)),
  pct_zero = 100 * c(mean(dat$n_emergency == 0), mean(dat$n_outpatient == 0)),
  skewness = sapply(list(dat$n_emergency, dat$n_outpatient), function(v)
    mean((v - mean(v))^3) / sd(v)^3)),
  "02_descriptive_skewed_covariates.csv")

## 2. All candidate subsets: CV errors, MSE, AIC, BIC, convergence ----
all_sub <- do.call(rbind, lapply(DISTS, function(dist)
  do.call(rbind, lapply(seq_len(NM), function(j) {
    r <- res[[paste(dist, j)]]
    data.frame(model = dist, j = j, subset = sub_name(SUBSETS[[j]]),
               size = length(SUBSETS[[j]]), n_par = r$p, loglik = r$ll,
               AIC = -2 * r$ll + 2 * r$p, BIC = -2 * r$ll + log(n) * r$p,
               sqrt_mse = sqrt(r$mse),
               sqrt_cv5 = sqrt(mean(r$E5)), sqrt_cv10 = sqrt(mean(r$E10)),
               converged = r$ok)
  }))))
wcsv(all_sub, "03_all_subsets.csv")

## 3. Selected models (Table app_results) -----------------------------
sel_tab <- do.call(rbind, lapply(names(SEL), function(kk) {
  tb <- SEL[[kk]]$tab
  data.frame(k = kk, model = tb$model, subset = tb$subset,
             sqrt_eps_opt = tb$sqrt_eps_opt, sqrt_mse = tb$sqrt_mse, BIC = tb$BIC)
}))
wcsv(sel_tab, "04_selected_models.csv")

## 4. Friedman test (uncorrected and corrected) -----------------------
fr_tab <- do.call(rbind, lapply(names(SEL), function(kk) {
  fr <- SEL[[kk]]$fr; k <- as.integer(kk)
  data.frame(k = k, Q = fr["Q"], p_uncorrected = fr["p_unc"],
             Fc = fr["Fc"], p_corrected = fr["p_cor"],
             c_k = (2 * k - 1) / k, k_eff = k^2 / (2 * k - 1),
             rank_sum_P  = colSums(t(apply(SEL[[kk]]$E, 1, rank)))[1],
             rank_sum_NB = colSums(t(apply(SEL[[kk]]$E, 1, rank)))[2],
             rank_sum_GP = colSums(t(apply(SEL[[kk]]$E, 1, rank)))[3])
}))
wcsv(fr_tab, "05_friedman.csv")

## 5. Pairwise comparisons (Table app_pairwise), always computed ------
pw_tab <- do.call(rbind, lapply(names(SEL), function(kk) {
  E <- SEL[[kk]]$E; pw <- pairwise_corrected(E)
  pr <- combn(ncol(E), 2)
  pw$mean_diff <- apply(pr, 2, function(ab) mean(E[, ab[1]] - E[, ab[2]]))
  cbind(k = kk, pw)
}))
wcsv(pw_tab, "06_pairwise.csv")

## 6. Fold errors of the selected models (Figure figura2) -------------
fold_tab <- do.call(rbind, lapply(names(SEL), function(kk) {
  E <- SEL[[kk]]$E; Rk <- t(apply(E, 1, rank))
  data.frame(k = kk, fold = seq_len(nrow(E)),
             E_P = E[, 1], E_NB = E[, 2], E_GP = E[, 3],
             rank_P = Rk[, 1], rank_NB = Rk[, 2], rank_GP = Rk[, 3])
}))
wcsv(fold_tab, "07_fold_errors.csv")

## 7. Checks cited in the text ----------------------------------------
j_nb <- SEL[["10"]]$tab$j[2]
chk <- data.frame(
  check = "P with the subset selected for NB",
  model = "P", subset = sub_name(SUBSETS[[j_nb]]),
  BIC = -2 * get("P", j_nb, "ll") + log(n) * get("P", j_nb, "p"),
  sqrt_cv5  = sqrt(mean(get("P", j_nb, "E5"))),
  sqrt_cv10 = sqrt(mean(get("P", j_nb, "E10"))))
for (dist in DISTS) {
  ll <- sapply(seq_len(NM), function(j) get(dist, j, "ll"))
  p  <- sapply(seq_len(NM), function(j) get(dist, j, "p"))
  for (crit in c("AIC", "BIC")) {
    pen <- if (crit == "AIC") 2 else log(n)
    jb <- which.min(-2 * ll + pen * p)
    chk <- rbind(chk, data.frame(
      check = paste("full-data", crit, "selection"), model = dist,
      subset = sub_name(SUBSETS[[jb]]), BIC = -2 * ll[jb] + log(n) * p[jb],
      sqrt_cv5 = sqrt(mean(get(dist, jb, "E5"))),
      sqrt_cv10 = sqrt(mean(get(dist, jb, "E10")))))
  }
}
wcsv(chk, "08_checks_benchmarks.csv")

## 8. Observed vs expected frequencies (Table app_freq) ---------------
freq_out <- freq
freq_out$y <- as.character(freq_out$y)
pearson <- sapply(DISTS, function(d) sum((freq$Observed - freq[[d]])^2 / freq[[d]]))
freq_out <- rbind(freq_out,
                  data.frame(y = "Pearson", Observed = NA,
                             P = pearson["P"], NB = pearson["NB"], GP = pearson["GP"]))
wcsv(freq_out, "09_frequencies.csv")
wcsv(data.frame(model = DISTS,
                max_fitted_mean = sapply(final, function(f) max(f$mu)),
                phi = sapply(final, function(f) f$phi)),
     "10_fitted_means_phi.csv")

## 9. Estimates of the selected models (Table app_nb) -----------------
est_all <- do.call(rbind, lapply(DISTS, function(dist)
  cbind(model = dist, final_model = (dist == best), est_table(dist))))
wcsv(est_all, "11_estimates.csv")

## 10. Session information --------------------------------------------
wcsv(data.frame(item = c("R version", "nloptr version", "seed", "n", "final model"),
                value = c(R.version.string, as.character(packageVersion("nloptr")),
                          SEED, n, best)),
     "12_session_info.csv")

cat("\nAll results saved in folder", OUT, "\n")
