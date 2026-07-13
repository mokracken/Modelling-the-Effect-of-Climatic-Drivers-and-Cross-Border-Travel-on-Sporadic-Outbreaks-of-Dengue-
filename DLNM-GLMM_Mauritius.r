############################################################
# Title: Modelling the Effect of Climatic Drivers and Cross-Border Travel on Sporadic Outbreaks of Dengue in Mauritius

# Aim: To quantify the association between climatic variables (temperature, rainfall and vapor pressure deficit),
# international arrivals and dengue incidence in Mauritius


# Author: Mohabeer Teeluck, Oyelola Adegboye*, Benn Sartorius, Eloise Skinner, Stephan Karl, Alec Henderson, Dr. M.F. Khodaboccus, Diana P. Iyaloo, Emma McBryde*
# *Lead supervisors: Oyelola Adegboye and Emma McBryde
# Institution: University of Queensland, James Cook University and Charles Darwin University
# Date: 27/04/2026
# R version: 4.4.1

# ==============================================================================
# 1) Libraries
# ==============================================================================
suppressPackageStartupMessages({
  library(glmmTMB)
  library(dlnm)
  library(DHARMa)
  library(splines)
  library(dplyr)
  library(ggplot2)
  library(car)         
  library(corrplot)    
  library(performance) 
  library(mvtnorm)
  library(tidyr)
})

set.seed(202) 

# ==============================================================================
# 2) Helper functions
# ==============================================================================

compute_vif <- function(data, vars) {
  formula_str <- paste("local_cases ~", paste(vars, collapse = " + "))
  lm_model <- lm(as.formula(formula_str), data = data)
  print(round(car::vif(lm_model), 2))
}

evaluate_model <- function(model) {
  rp       <- residuals(model, type = "pearson")
  phi_raw  <- sum(rp^2, na.rm = TRUE) / df.residual(model)
  if (is.finite(phi_raw) && phi_raw > 2) return(NULL)
  
  df_ll <- attr(logLik(model), "df")
  qaic  <- AIC(model) + 2 * (phi_raw - 1) * df_ll
  qbic  <- BIC(model) + log(nrow(model$frame)) * (phi_raw - 1) * df_ll
  
  lag1_acf <- tryCatch(acf(rp, lag.max = 1, plot = FALSE)$acf[2], error = function(e) NA_real_)
  
  dh_p_disp <- dh_p_uni <- NA_real_
  try({
    set.seed(123)
    sim <- DHARMa::simulateResiduals(
      fittedModel = model, n = 1000, refit = FALSE, integerResponse = TRUE,
      simulateFunction = simulate, simulateFunctionArgs = list(cond = TRUE)
    )
    dh_p_disp <- as.numeric(DHARMa::testDispersion(sim)$p.value)
    dh_p_uni  <- as.numeric(DHARMa::testUniformity(sim)$p.value)
  }, silent = TRUE)
  
  r2_pred <- tryCatch({
    mf <- model$frame
    cor(mf$local_cases, fitted(model, type = "response"), use = "complete.obs")^2
  }, error = function(e) NA_real_)
  
  r2_marg <- r2_cond <- NA_real_
  r2obj <- tryCatch(performance::r2_nakagawa(model, tolerance = 0), error = function(e) NULL)
  if (!is.null(r2obj)) {
    if (!is.null(r2obj$R2_marginal))    r2_marg <- as.numeric(r2obj$R2_marginal)
    if (!is.null(r2obj$R2_conditional)) r2_cond <- as.numeric(r2obj$R2_conditional)
    if (is.na(r2_marg) || is.na(r2_cond)) {
      df <- tryCatch(as.data.frame(r2obj), error = function(e) NULL)
      if (!is.null(df)) {
        if ("R2_marginal" %in% names(df))    r2_marg <- as.numeric(df$R2_marginal[1])
        if ("R2_conditional" %in% names(df)) r2_cond <- as.numeric(df$R2_conditional[1])
      }
    }
  }
  
  list(
    QAIC = round(qaic, 3),
    QBIC = round(qbic, 3),
    Overdispersion = round(phi_raw, 4),
    DHARMa_p_dispersion = dh_p_disp,
    DHARMa_p_uniformity = dh_p_uni,
    Lag1_ACF = round(lag1_acf, 4),
    R2_pred = r2_pred,
    R2_Nakagawa_marginal = ifelse(is.na(r2_marg), NA_real_, round(r2_marg, 4)),
    R2_Nakagawa_conditional = ifelse(is.na(r2_cond), NA_real_, round(r2_cond, 4))
  )
}

plot_obs_vs_pred <- function(model, trim_q = 0.99, save_path = NULL, show = TRUE) {
  mf_used <- model$frame
  nat_pred_in <- data.frame(
    observed  = mf_used$local_cases,
    predicted = fitted(model, type = "response")
  )
  rmse_in <- sqrt(mean((nat_pred_in$observed - nat_pred_in$predicted)^2, na.rm = TRUE))
  r2_in   <- cor(nat_pred_in$observed, nat_pred_in$predicted, use = "complete.obs")^2
  ann     <- sprintf("RMSE = %.1f | R² = %.2f", rmse_in, r2_in)
  
  x_cut <- quantile(nat_pred_in$predicted, trim_q, na.rm = TRUE)
  y_cut <- quantile(nat_pred_in$observed,  trim_q, na.rm = TRUE)
  lim   <- c(0, max(x_cut, y_cut, na.rm = TRUE))
  
  idx_trim <- nat_pred_in$predicted <= x_cut & nat_pred_in$observed <= y_cut
  cal_fit  <- lm(observed ~ predicted, data = nat_pred_in[idx_trim, ])
  label_tag <- if (trim_q < 1) "Calibration (trimmed)" else "Calibration (untrimmed)"
  cal_lab  <- sprintf("%s: y = %.2f + %.2f·x",
                      label_tag, unname(coef(cal_fit)[1]), unname(coef(cal_fit)[2]))
  
  p <- ggplot(nat_pred_in, aes(x = predicted, y = observed)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(data = subset(nat_pred_in, idx_trim), method = "lm", se = FALSE, linewidth = 0.8) +
    labs(title = "Observed vs Predicted (National, in-sample)",
         x = "Predicted cases", y = "Observed cases") +
    annotate("text", x = lim[2], y = lim[2], label = ann, hjust = 1.05, vjust = 1.6, size = 4) +
    annotate("text", x = lim[2], y = lim[2], label = cal_lab, hjust = 1.05, vjust = 3.0, size = 3.7) +
    coord_fixed(xlim = lim, ylim = lim, ratio = 1) +
    theme_minimal(base_size = 13)
  
  if (!is.null(save_path)) {
    if (!dir.exists(dirname(save_path))) dir.create(dirname(save_path), recursive = TRUE)
    ggsave(save_path, p, width = 6, height = 5.4, dpi = 300, bg = "white")
  }
  if (show) print(p)
  invisible(NULL)
}

mrv_point <- function(cb, coef_sub, vcov_sub, xvec, within = c(0.05, 0.95), by = NULL) {
  if (is.null(cb)) return(NA_real_)
  qlo <- as.numeric(quantile(xvec, within[1], na.rm = TRUE))
  qhi <- as.numeric(quantile(xvec, within[2], na.rm = TRUE))
  if (!is.finite(qlo) || !is.finite(qhi) || qlo >= qhi) return(NA_real_)
  if (is.null(by)) { rng <- qhi - qlo; by <- max(rng / 200, rng / 1000) }
  at_grid <- seq(qlo, qhi, by = by)
  cen0 <- median(xvec[xvec >= qlo & xvec <= qhi], na.rm = TRUE); if (!is.finite(cen0)) cen0 <- (qlo + qhi) / 2
  cp <- tryCatch(crosspred(cb, coef = coef_sub, vcov = vcov_sub, model.link = "log", at = at_grid, cen = cen0), error = function(e) NULL)
  if (is.null(cp)) return(NA_real_)
  rr <- as.numeric(cp$allRRfit[as.character(at_grid)])
  if (all(!is.finite(rr))) return(NA_real_)
  at_grid[which.min(rr)]
}

mrv_ci_se <- function(cb, coef_sub, vcov_sub, xvec,
                      within = c(0.05, 0.95), by = NULL,
                      nsim = 1000, conf = c(0.025, 0.975)) {
  if (is.null(cb)) return(list(mrv = NA_real_, ci = c(NA_real_, NA_real_), se = NA_real_))
  mrv_hat <- mrv_point(cb, coef_sub, vcov_sub, xvec, within = within, by = by)
  sims <- tryCatch(mvtnorm::rmvnorm(nsim, mean = as.numeric(coef_sub), sigma = as.matrix(vcov_sub)), error = function(e) NULL)
  if (is.null(sims)) return(list(mrv = mrv_hat, ci = c(NA_real_, NA_real_), se = NA_real_))
  qlo <- as.numeric(quantile(xvec, within[1], na.rm = TRUE))
  qhi <- as.numeric(quantile(xvec, within[2], na.rm = TRUE))
  if (is.null(by)) { rng <- qhi - qlo; by <- max(rng / 200, rng / 1000) }
  at_grid <- seq(qlo, qhi, by = by)
  cen0 <- median(xvec[xvec >= qlo & xvec <= qhi], na.rm = TRUE); if (!is.finite(cen0)) cen0 <- (qlo + qhi) / 2
  mrv_s <- vapply(seq_len(nsim), function(i) {
    cp_i <- tryCatch(crosspred(cb, coef = sims[i, ], vcov = vcov_sub, model.link = "log", at = at_grid, cen = cen0), error = function(e) NULL)
    if (!is.null(cp_i)) {
      rr_i <- as.numeric(cp_i$allRRfit[as.character(at_grid)])
      if (any(is.finite(rr_i))) return(at_grid[which.min(rr_i)])
    }
    NA_real_
  }, FUN.VALUE = numeric(1))
  mrv_s <- mrv_s[is.finite(mrv_s)]
  ci <- if (length(mrv_s) > 10) quantile(mrv_s, conf, na.rm = TRUE) else c(NA_real_, NA_real_)
  se <- if (length(mrv_s) > 1) sd(mrv_s, na.rm = TRUE) else NA_real_
  list(mrv = mrv_hat, ci = unname(ci), se = se)
}

summarize_rr <- function(var_label, vec, cb, coef_all, vcov_all,
                         center_fn,
                         prob_pairs = list(c(0.05, 0.95), c(0.04, 0.96), c(0.06, 0.94))) {
  ind <- grep(paste0("^", var_label), names(coef_all)); stopifnot(length(ind) > 0)
  cen <- if (is.function(center_fn)) center_fn(vec, na.rm = TRUE) else as.numeric(center_fn)
  is_bad_row <- function(rr, low, high) {
    any(!is.finite(c(rr, low, high))) || (abs(rr - 1) < 1e-12 && abs(low - 1) < 1e-12 && abs(high - 1) < 1e-12)
  }
  used <- NULL; pred_ok <- NULL
  for (pair in prob_pairs) {
    p_lo <- unname(as.numeric(quantile(vec, pair[1], na.rm = TRUE)))
    p_hi <- unname(as.numeric(quantile(vec, pair[2], na.rm = TRUE)))
    if (is.finite(cen) && (abs(p_lo - cen) < 1e-12 || abs(p_hi - cen) < 1e-12)) next
    pred_try <- crosspred(cb, coef = coef_all[ind], vcov = vcov_all[ind, ind], model.link = "log", cen = cen, at = c(p_lo, p_hi))
    rr_lo <- unname(pred_try$allRRfit[as.character(p_lo)]); lcl_lo <- unname(pred_try$allRRlow[as.character(p_lo)]); ucl_lo <- unname(pred_try$allRRhigh[as.character(p_lo)])
    rr_hi <- unname(pred_try$allRRfit[as.character(p_hi)]); lcl_hi <- unname(pred_try$allRRlow[as.character(p_hi)]); ucl_hi <- unname(pred_try$allRRhigh[as.character(p_hi)])
    if (!is_bad_row(rr_lo, lcl_lo, ucl_lo) && !is_bad_row(rr_hi, lcl_hi, ucl_hi)) { used <- pair; pred_ok <- pred_try; break }
  }
  if (is.null(used)) return(tibble::tibble(variable = var_label, p_lo_prob = NA, p_hi_prob = NA))
  pull_row <- function(x, prob_used, side_lab) {
    rr <- unname(pred_ok$allRRfit[as.character(x)]); low <- unname(pred_ok$allRRlow[as.character(x)]); high <- unname(pred_ok$allRRhigh[as.character(x)])
    logRR <- log(rr); SE <- (log(high) - log(low)) / (2 * 1.96)
    z <- if (is.finite(SE) && SE > 0) logRR / SE else NA_real_
    p <- if (is.finite(z)) 2 * pnorm(abs(z), lower.tail = FALSE) else NA_real_
    tibble::tibble(
      !!paste0("p_", side_lab, "_prob") := prob_used, !!paste0("p_", side_lab, "_value") := x,
      !!paste0("RR_", side_lab) := rr, !!paste0("LCL_", side_lab) := low, !!paste0("UCL_", side_lab) := high, !!paste0("p_", side_lab) := p
    )
  }
  p_lo <- unname(as.numeric(quantile(vec, used[1], na.rm = TRUE)))
  p_hi <- unname(as.numeric(quantile(vec, used[2], na.rm = TRUE)))
  dplyr::bind_cols(
    tibble::tibble(variable = var_label, center_value = cen),
    pull_row(p_lo, used[1], "lo"),
    pull_row(p_hi, used[2], "hi")
  )
}

lagwise_rr <- function(var_label, vec, cb, coef_all, vcov_all,
                       center_fn, probs_used = NULL,
                       prob_pairs = list(c(0.05, 0.95), c(0.04, 0.96), c(0.06, 0.94)),
                       bylag = 0.25) {
  ind <- grep(paste0("^", var_label), names(coef_all)); stopifnot(length(ind) > 0)
  cen <- if (is.function(center_fn)) center_fn(vec, na.rm = TRUE) else as.numeric(center_fn)
  if (is.null(probs_used)) {
    for (pair in prob_pairs) {
      p_lo <- unname(as.numeric(quantile(vec, pair[1], na.rm = TRUE)))
      p_hi <- unname(as.numeric(quantile(vec, pair[2], na.rm = TRUE)))
      if (abs(p_lo - cen) > 1e-12 && abs(p_hi - cen) > 1e-12) { probs_used <- pair; break }
    }
  }
  if (is.null(probs_used)) stop("Could not find valid percentiles for lagwise_rr().")
  p_lo <- unname(as.numeric(quantile(vec, probs_used[1], na.rm = TRUE)))
  p_hi <- unname(as.numeric(quantile(vec, probs_used[2], na.rm = TRUE)))
  pred <- crosspred(cb, coef = coef_all[ind], vcov = vcov_all[ind, ind], model.link = "log", cen = cen, at = c(p_lo, p_hi), bylag = bylag)
  z_cols <- ncol(pred$matRRfit); raw_lag <- pred$lag
  lag_grid <- if (length(raw_lag) == z_cols) raw_lag else seq(min(raw_lag, na.rm = TRUE), max(raw_lag, na.rm = TRUE), length.out = z_cols)
  nearest_row <- function(x, label) {
    i <- which.min(abs(pred$predvar - x))
    tibble::tibble(
      variable = var_label, exposure_quantile = label,
      at_value = pred$predvar[i], lag = lag_grid,
      RR = as.numeric(pred$matRRfit[i, ]), LCL = as.numeric(pred$matRRlow[i, ]), UCL = as.numeric(pred$matRRhigh[i, ])
    ) |>
      dplyr::mutate(
        logRR = log(RR), SE = (log(UCL) - log(LCL)) / (2 * 1.96),
        z = dplyr::if_else(is.finite(SE) & SE > 0, logRR / SE, NA_real_),
        p = dplyr::if_else(is.finite(z), 2 * pnorm(abs(z), lower.tail = FALSE), NA_real_),
        p_used = ifelse(label == "p05", probs_used[1], probs_used[2])
      )
  }
  dplyr::bind_rows(nearest_row(p_lo, "p05"), nearest_row(p_hi, "p95"))
}

xlab_for <- function(var)
  switch(var,
         tmin          = "Minimum temperature (°C)",
         rain          = "Rainfall (mm)",
         VPD           = "VPD (kPa)",
         res_arr_scaled = "International resident arrivals (×1000)",
         var
  )

title_for <- function(var) paste0("Relative Risk (RR) vs ", xlab_for(var))

plot_dlnm_outputs <- function(coef_all, vcov_all, glmm_model_final, statistic,
                              data, cross_basis_var, var, cb, ylim_upper = NULL) {
  indt <- grep(cross_basis_var, names(coef_all))
  if (length(indt) == 0) stop(paste0("No coefficients matched for ", cross_basis_var, "."))
  cen <- if (is.function(statistic)) statistic(data[[var]], na.rm = TRUE) else as.numeric(statistic)
  if (!is.finite(cen)) stop(paste("Statistic (cen) for", var, "is NA."))
  pred <- crosspred(cb, coef = coef_all[indt], vcov = vcov_all[indt, indt, drop = FALSE], model.link = "log", bylag = 0.25, cen = cen)
  
  z_matrix <- pred$matRRfit; x_exposure <- pred$predvar; raw_lag <- pred$lag
  y_lag <- if (length(raw_lag) == ncol(z_matrix)) raw_lag else seq(min(raw_lag, na.rm = TRUE), max(raw_lag, na.rm = TRUE), length.out = ncol(z_matrix))
  pal <- c("#006400", "gold", "red"); levels <- pretty(z_matrix, 20)
  col1 <- colorRampPalette(pal[1:2]); col2 <- colorRampPalette(pal[2:3])
  cols <- c(col1(sum(levels < 1)), col2(sum(levels > 1)))
  
  xlab_text <- xlab_for(var)
  main_str <- paste0("Relative Risk (RR)\nfor ", xlab_text)
  
  filled.contour(
    x = x_exposure, y = y_lag, z = z_matrix,
    col = cols, levels = levels,
    plot.title = title(
      main = main_str,
      xlab = xlab_text,
      ylab = "Lag (weeks)",
      cex.main = par("cex.lab"),
      font.main = 2
    ),
    key.title = title(main = "RR", cex.main = par("cex.lab"), font.main = 2)
  )
  
  rr_df <- data.frame(value = pred$predvar, rr = pred$allRRfit, low = pred$allRRlow, high = pred$allRRhigh) |>
    dplyr::filter(is.finite(rr) & rr > 0) |>
    dplyr::bind_rows(data.frame(value = cen, rr = 1, low = 1, high = 1)) |>
    dplyr::arrange(value)
  
  q_lo_prob <- 0.05; q_hi_prob <- 0.95
  q05_raw <- as.numeric(quantile(data[[var]], q_lo_prob, na.rm = TRUE))
  q95_raw <- as.numeric(quantile(data[[var]], q_hi_prob, na.rm = TRUE))
  if (identical(var, "VPD") && is.finite(cen) && abs(q05_raw - cen) < 1e-12) {
    q_lo_prob <- 0.04; q_hi_prob <- 0.96
    q05_raw <- as.numeric(quantile(data[[var]], q_lo_prob, na.rm = TRUE))
    q95_raw <- as.numeric(quantile(data[[var]], q_hi_prob, na.rm = TRUE))
  }
  
  pred_q <- crosspred(cb, coef = coef_all[indt], vcov = vcov_all[indt, indt, drop = FALSE], model.link = "log", cen = cen, at = c(q05_raw, q95_raw))
  p05_df <- tibble::tibble(value = q05_raw,
                           RR = unname(pred_q$allRRfit[as.character(q05_raw)]),
                           LCL = unname(pred_q$allRRlow[as.character(q05_raw)]),
                           UCL = unname(pred_q$allRRhigh[as.character(q05_raw)]))
  p95_df <- tibble::tibble(value = q95_raw,
                           RR = unname(pred_q$allRRfit[as.character(q95_raw)]),
                           LCL = unname(pred_q$allRRlow[as.character(q95_raw)]),
                           UCL = unname(pred_q$allRRhigh[as.character(q95_raw)]))
  pt_df <- dplyr::bind_rows(p05_df, p95_df)
  mrv_marker <- tibble::tibble(value = cen, RR = 1)
  
  cap <- if (is.null(ylim_upper)) Inf else ylim_upper
  rr_df <- rr_df %>% dplyr::mutate(low = pmax(low, 0), high = pmin(high, cap))
  
  p_overall <- ggplot2::ggplot(rr_df, ggplot2::aes(x = value)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = low, ymax = high), alpha = 0.25, fill = "grey80") +
    ggplot2::geom_line(ggplot2::aes(y = rr), linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = cen, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = c(q05_raw, q95_raw), linetype = "dotted") +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::geom_point(data = mrv_marker, ggplot2::aes(x = value, y = RR), shape = 21, size = 3.3, stroke = 0.9, fill = NA, color = "black") +
    ggplot2::geom_point(data = mrv_marker, ggplot2::aes(x = value, y = RR), shape = 4, size = 3.3, stroke = 0.9, color = "black") +
    ggplot2::geom_errorbar(data = pt_df, ggplot2::aes(x = value, ymin = LCL, ymax = UCL), width = 0) +
    ggplot2::geom_point(data = p05_df, ggplot2::aes(x = value, y = RR), color = "blue", size = 2.5) +
    ggplot2::geom_point(data = p95_df, ggplot2::aes(x = value, y = RR), color = "red", size = 2.5) +
    ggplot2::labs(title = title_for(var), x = xlab_text, y = "Relative Risk (RR)") +
    ggplot2::coord_cartesian(ylim = c(0, cap)) +
    ggplot2::theme_minimal(base_size = 13)
  
  print(p_overall)
  invisible(pred)
}

plot_lagcurve_publication <- function(d, var_label) {
  p <- ggplot(d, aes(x = lag, y = RR, color = exposure_quantile, fill = exposure_quantile)) +
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.20, colour = NA) +
    geom_line(size = 1) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_color_manual(values = c("p05" = "blue", "p95" = "red")) +
    scale_fill_manual(values = c("p05" = "blue", "p95" = "red")) +
    labs(
      title = title_for(var_label),
      x = "Lag (weeks)",
      y = "Relative Risk (RR)",
      color = "Quantile",
      fill  = "Quantile"
    ) +
    theme_minimal(base_size = 13)
  print(p)
  invisible(p)
}

# ==============================================================================
# 3) Load & preprocess data
# ==============================================================================
csv_path <- Sys.getenv("DLNM_DATA", unset = NA)
if (is.na(csv_path) || !nzchar(csv_path)) {
  if (interactive()) { csv_path <- file.choose() } else {
    stop("No CSV specified. Set env var DLNM_DATA to your CSV path or run interactively to pick a file.")
  }
}
data <- read.csv(csv_path, header = TRUE, check.names = TRUE)
cat("Loaded:", csv_path, "\n")

data <- data %>%
  dplyr::filter(
    !is.na(population) & population > 0, !is.na(local_cases), !is.na(tmin),
    !is.na(rain), !is.na(VPD), !is.na(imported_cases), !is.na(res_arr),
    !is.na(state_name), year >= 2014, year <= 2024
  ) %>%
  dplyr::arrange(year, week) %>%
  dplyr::mutate(week = as.integer(week), year = as.integer(year)) %>%
  dplyr::filter((year < 2020) | (year == 2020 & week < 14) | (year == 2021 & week > 40) | (year > 2021)) %>%
  dplyr::group_by(state_name) %>%
  dplyr::mutate(
    month = as.integer(month),
    local_cases = as.integer(local_cases),
    imported_cases = as.numeric(imported_cases),
    res_arr = as.numeric(res_arr),
    res_arr_scaled = res_arr / 1000,
    lag_imported = dplyr::lag(imported_cases, 1),
    lagged_local_cases = dplyr::lag(local_cases, 1),
    tmin = as.numeric(tmin),
    VPD = as.numeric(VPD),
    state_name = as.factor(state_name),
    time_index = dplyr::row_number(),
    outbreak_level = dplyr::case_when(
      year %in% c(2020, 2023, 2024) ~ "high",
      year %in% c(2014, 2015, 2019) ~ "moderate",
      TRUE ~ "none"
    ),
    outbreak_level = factor(outbreak_level, levels = c("none", "moderate", "high"))
  ) %>%
  dplyr::ungroup() %>%
  tidyr::drop_na() %>%
  dplyr::mutate(year = factor(year))
cat("Rows after cleaning:", nrow(data), "\n")

# ==============================================================================
# 4) Quick diagnostics
# ==============================================================================
vars <- data %>% dplyr::select(tmin, rain, VPD, res_arr_scaled, lag_imported, lagged_local_cases, time_index)
cor_matrix <- cor(vars, use = "complete.obs"); print(round(cor_matrix, 2))
corrplot(cor_matrix, method = "color", type = "lower", addCoef.col = "black", tl.col = "black", number.cex = 0.8, diag = FALSE)
compute_vif(data, vars = c("tmin", "rain", "VPD", "res_arr_scaled", "lagged_local_cases", "lag_imported", "time_index"))

# ==============================================================================
# 5) Cross basis
# ==============================================================================
cb_specs <- list(
  tmin   = list(var = "tmin",             lab = "cb_tmin",     lag =  6, df_var = 3, df_lag = 3, within = c(0.05, 0.95), by = 0.05,  ylim = 40),
  rain   = list(var = "rain",             lab = "cb_rain",     lag =  8, df_var = 4, df_lag = 4, within = c(0.05, 0.95), by = 0.5,   ylim = 55),
  VPD    = list(var = "VPD",              lab = "cb_vpd",      lag = 10, df_var = 4, df_lag = 3, within = c(0.05, 0.95), by = 0.2,   ylim = 15),
  res_arr= list(var = "res_arr_scaled",   lab = "cb_res_arr",  lag =  2, df_var = 2, df_lag = 2, within = c(0.10, 0.90), by = 0.005, ylim = 10)
)
make_cb <- function(x, lag, df_var, df_lag) {
  crossbasis(x, lag = lag, argvar = list(fun = "ns", df = df_var), arglag = list(fun = "ns", df = df_lag))
}
cb_list <- lapply(cb_specs, function(s) make_cb(data[[s$var]], s$lag, s$df_var, s$df_lag))
names(cb_list) <- vapply(cb_specs, `[[`, "", "lab")
list2env(cb_list, envir = environment())  # creates cb_tmin, cb_rain, cb_vpd, cb_res_arr

# ==============================================================================
# 6) Fit DLNM-GLMM model
# ==============================================================================
glmm_model_final <- glmmTMB(
  local_cases ~ cb_tmin + cb_rain + cb_vpd + cb_res_arr +
    ns(time_index, df = 11) +
    outbreak_level + lagged_local_cases + lag_imported +
    (1 | state_name) + (1 | year),
  family = nbinom1(),
  offset = log(population),
  data = data,
  control = glmmTMBControl(
    optimizer = optim,
    optArgs   = list(method = "BFGS"),
    optCtrl   = list(maxit = 2000, reltol = 1e-8)
  ))
print(summary(glmm_model_final))
cat("Convergence code:", glmm_model_final$fit$convergence, "\n")
if (!is.null(glmm_model_final$fit$message)) cat("Optimizer message:", glmm_model_final$fit$message, "\n")
cat("pdHess (Hessian positive-definite):", isTRUE(glmm_model_final$sdr$pdHess), "\n")
print(evaluate_model(glmm_model_final))

plot_obs_vs_pred(glmm_model_final, trim_q = 0.99,
                 save_path = "outputs/obs_vs_pred_in_sample.png", show = TRUE)

plot_obs_vs_pred(glmm_model_final, trim_q = 1.00,
                 save_path = "outputs/obs_vs_pred_in_sample_UNTRIMMED.png", show = TRUE)

# ==============================================================================
# 7) MRV per exposure 
# ==============================================================================
coef_all <- fixef(glmm_model_final)$cond
vcov_all <- vcov(glmm_model_final)$cond

get_mrv <- function(s) {
  ind <- grep(paste0("^", s$lab), names(coef_all))
  out <- mrv_ci_se(get(s$lab), coef_all[ind], vcov_all[ind, ind, drop = FALSE],
                   data[[s$var]], within = s$within, by = s$by, nsim = 1000)
  data.frame(variable = s$var, MRV = out$mrv, LCL = out$ci[1], UCL = out$ci[2], SE = out$se)
}
mrv_table <- do.call(rbind, lapply(cb_specs, get_mrv))
print(mrv_table)
write.csv(mrv_table, file = file.path("outputs", "MRV_point_CI_SE.csv"), row.names = FALSE)

# ==============================================================================
# 8) Plots for each exposure (centered at MRV)
# ==============================================================================
for (nm in names(cb_specs)) {
  s <- cb_specs[[nm]]
  cen <- mrv_table$MRV[mrv_table$variable == s$var]
  plot_dlnm_outputs(coef_all, vcov_all, glmm_model_final,
                    statistic = cen, data, s$lab, s$var, get(s$lab), ylim_upper = s$ylim)
}

# ==============================================================================
# 9) Overall RR summaries at p05 & p95 (center = MRV)
# ==============================================================================
summ_one <- function(s) {
  cen <- mrv_table$MRV[mrv_table$variable == s$var]
  summarize_rr(s$lab, data[[s$var]], get(s$lab), coef_all, vcov_all, center_fn = cen)
}
rr_list <- lapply(cb_specs, summ_one)
rr_overall <- dplyr::bind_rows(rr_list)
print(rr_overall)
write.csv(rr_overall, file = file.path("outputs", "overall_RR_summary_center_MRV.csv"), row.names = FALSE)

# ==============================================================================
# 10) Lag-wise RR tables at p05 & p95 (center = MRV)
# ==============================================================================
get_probs <- function(rr_row) c(rr_row$p_lo_prob[1], rr_row$p_hi_prob[1])
probs_used <- lapply(rr_list, get_probs)

lag_one <- function(s, probs) {
  cen <- mrv_table$MRV[mrv_table$variable == s$var]
  lagwise_rr(s$lab, data[[s$var]], get(s$lab), coef_all, vcov_all, center_fn = cen, probs_used = probs)
}
lag_df_list <- Map(lag_one, cb_specs, probs_used)
lag_rr_all <- dplyr::bind_rows(lag_df_list)
write.csv(lag_rr_all, file = file.path("outputs", "lagwise_RR_p05_p95_center_MRV.csv"), row.names = FALSE)

for (nm in names(cb_specs)) {
  s <- cb_specs[[nm]]
  d <- lag_rr_all %>%
    dplyr::filter(variable == s$lab) %>%
    dplyr::mutate(exposure_quantile = factor(exposure_quantile, levels = c("p05","p95")))
  if (nrow(d) > 0) plot_lagcurve_publication(d, var_label = s$var)
}
