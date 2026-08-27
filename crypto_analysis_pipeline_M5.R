#############################################################################
# Cryptocurrency Analysis Pipeline - FIXED VERSION
# Addresses the "row names contain missing values" error in Section 21b
#############################################################################

# ---------------------------------------------------------------------------
# 0. Setup
# ---------------------------------------------------------------------------
required_packages <- c("quantmod", "TTR", "dplyr", "tidyr", "lubridate", "zoo",
                       "e1071", "tseries", "signal", "car", "MASS", "lmtest",
                       "rpart", "ranger", "gbm", "randomForest", "xgboost",
                       "glmnet", "caret", "rBayesianOptimization", "ncvreg")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "https://cloud.r-project.org")

# Load libraries
library(quantmod)
library(TTR)
library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)
library(e1071)
library(tseries)
library(signal)
library(car)
library(MASS)
library(lmtest)
library(rpart)
library(ranger)
library(gbm)
library(randomForest)
library(xgboost)
library(glmnet)
library(caret)
library(rBayesianOptimization)
library(ncvreg)

# Fix search path issues: reattach dplyr last so its generics win
if ("package:dplyr" %in% search()) suppressWarnings(detach("package:dplyr", unload = TRUE))
library(dplyr)

# Create output directory
out_dir <- "crypto_output"
if (!dir.exists(out_dir)) dir.create(out_dir)

# Load data collection + feature engineering
if (file.exists("crypto_feature_engineering.R")) {
  source("crypto_feature_engineering.R")
} else {
  stop("crypto_feature_engineering.R not found in working directory.")
}

# Define coin groups
GROUP1 <- c("BTC", "ETH", "BCH", "DASH")
GROUP2 <- c("XLM", "XRP", "XMR", "LTC")

#############################################################################
# 1. Descriptive statistics of the response variable (Close)
#############################################################################
descriptive_stats <- function(df, value_col = "Close") {
  x <- na.omit(df[[value_col]])
  x_sw <- if (length(x) > 5000) sample(x, 5000) else x
  sw <- tryCatch(shapiro.test(x_sw), error = function(e) list(statistic = NA, p.value = NA))
  jb <- tryCatch(tseries::jarque.bera.test(x), error = function(e) list(statistic = NA, p.value = NA))
  
  data.frame(
    N               = length(x),
    Mean            = mean(x),
    SD              = sd(x),
    Skewness        = e1071::skewness(x, type = 2),
    Kurtosis        = e1071::kurtosis(x, type = 2),
    Shapiro_W       = unname(sw$statistic),
    Shapiro_p       = unname(sw$p.value),
    JarqueBera_stat = unname(jb$statistic),
    JarqueBera_p    = unname(jb$p.value)
  )
}

run_descriptive_stats_all_coins <- function(data_list, value_col = "Close") {
  bind_rows(lapply(names(data_list), function(nm) {
    stats <- descriptive_stats(data_list[[nm]], value_col)
    stats$Symbol <- nm
    stats %>% relocate(Symbol)
  }))
}

#############################################################################
# 3. Noise reduction: Butterworth low-pass filter + ACF comparison
#############################################################################
apply_butterworth_filter <- function(df, value_col = "Close", order = 3, fc = 0.05) {
  x <- df[[value_col]]
  x <- zoo::na.locf(x, na.rm = FALSE)
  x <- zoo::na.locf(x, fromLast = TRUE, na.rm = FALSE)
  bw <- signal::butter(n = order, W = fc, type = "low")
  filtered <- as.numeric(signal::filtfilt(bw, x))
  df[[paste0(value_col, "_Filtered")]] <- filtered
  df
}

compare_acf_filtered_unfiltered <- function(df, value_col = "Close", filtered_col = NULL,
                                            lag_max = 30, plot = TRUE) {
  if (is.null(filtered_col)) filtered_col <- paste0(value_col, "_Filtered")
  acf_unf <- acf(df[[value_col]], lag.max = lag_max, plot = FALSE)
  acf_flt <- acf(df[[filtered_col]], lag.max = lag_max, plot = FALSE)
  if (plot) {
    par(mfrow = c(1, 2))
    plot(acf_unf, main = paste("ACF - Unfiltered", value_col))
    plot(acf_flt, main = paste("ACF - Filtered", value_col))
    par(mfrow = c(1, 1))
  }
  list(acf_unfiltered = acf_unf, acf_filtered = acf_flt)
}

summarize_acf_comparison <- function(acf_result, short_lags = c(1, 2), long_lags = c(5, 10, 20)) {
  get_acf_at <- function(acf_obj, lag_k) {
    idx <- which(acf_obj$lag[, 1, 1] == lag_k)
    if (length(idx) == 0) return(NA_real_)
    acf_obj$acf[idx, 1, 1]
  }
  bind_rows(lapply(c(short_lags, long_lags), function(k) {
    data.frame(
      Lag = k,
      Lag_Type = ifelse(k %in% short_lags, "short (noise)", "long (structure)"),
      ACF_Unfiltered = get_acf_at(acf_result$acf_unfiltered, k),
      ACF_Filtered   = get_acf_at(acf_result$acf_filtered, k)
    )
  })) %>%
    mutate(Delta = ACF_Filtered - ACF_Unfiltered)
}

#############################################################################
# 4. Outlier identification: 1.5 x IQR rule
#############################################################################
detect_outliers_iqr <- function(x, k = 1.5) {
  q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  lower <- q[1] - k * iqr
  upper <- q[2] + k * iqr
  which(x < lower | x > upper)
}

outlier_report <- function(df, cols) {
  bind_rows(lapply(cols, function(cn) {
    idx <- detect_outliers_iqr(df[[cn]])
    data.frame(Feature = cn, N_Outliers = length(idx),
               Pct_Outliers = 100 * length(idx) / nrow(df))
  }))
}

#############################################################################
# 5. Temporal properties: ADF test, ARCH-LM test, ACF
#############################################################################
arch_lm_test <- function(x, lags = 12) {
  x <- na.omit(x)
  n <- length(x)
  resid_sq <- (x - mean(x))^2
  dat <- data.frame(y = resid_sq[(lags + 1):n])
  for (i in 1:lags) dat[[paste0("lag", i)]] <- resid_sq[(lags + 1 - i):(n - i)]
  fit <- lm(y ~ ., data = dat)
  r2 <- summary(fit)$r.squared
  LM_stat <- (n - lags) * r2
  list(statistic = LM_stat, p.value = 1 - pchisq(LM_stat, df = lags))
}

temporal_properties <- function(df, value_col = "Close") {
  x <- na.omit(df[[value_col]])
  ret <- diff(log(x))
  adf_res <- tryCatch(tseries::adf.test(x), error = function(e) NULL)
  arch_res <- tryCatch(arch_lm_test(ret, lags = 12), error = function(e) list(statistic = NA, p.value = NA))
  acf_res <- acf(ret, lag.max = 20, plot = FALSE)
  list(ADF_statistic = if (!is.null(adf_res)) adf_res$statistic else NA,
       ADF_p = if (!is.null(adf_res)) adf_res$p.value else NA,
       ARCH_LM_stat = arch_res$statistic,
       ARCH_LM_p = arch_res$p.value,
       ACF = acf_res)
}

run_temporal_properties_all <- function(data_list, value_col = "Close") {
  bind_rows(lapply(names(data_list), function(nm) {
    r <- temporal_properties(data_list[[nm]], value_col)
    data.frame(Symbol = nm, ADF_stat = r$ADF_statistic, ADF_p = r$ADF_p,
               ARCH_LM_stat = r$ARCH_LM_stat, ARCH_LM_p = r$ARCH_LM_p)
  }))
}

#############################################################################
# 6. Correlation analysis: predictors vs. response
#############################################################################
correlation_analysis <- function(df, predictors, response = "Close") {
  cors <- sapply(predictors, function(p) {
    tryCatch(cor(df[[p]], df[[response]], use = "pairwise.complete.obs"),
             error = function(e) NA)
  })
  data.frame(Feature = predictors, Correlation = cors) %>%
    arrange(desc(abs(Correlation)))
}

#############################################################################
# 7. Multicollinearity assessment: VIF and PCA
#############################################################################
vif_analysis <- function(df, predictors, response = "Close") {
  md <- df[, c(response, predictors)]
  md <- md[complete.cases(md), ]
  
  form <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  fit <- lm(form, data = md)
  
  al <- tryCatch(alias(fit), error = function(e) NULL)
  aliased_vars <- if (!is.null(al) && !is.null(al$Complete)) rownames(al$Complete) else character(0)
  usable_predictors <- setdiff(predictors, aliased_vars)
  
  vif_out <- setNames(rep(NA_real_, length(predictors)), predictors)
  
  if (length(usable_predictors) >= 2) {
    form2 <- as.formula(paste(response, "~", paste(usable_predictors, collapse = " + ")))
    fit2 <- lm(form2, data = md)
    vif_vals <- tryCatch(car::vif(fit2), error = function(e) {
      message("VIF computation failed: ", conditionMessage(e))
      NULL
    })
    if (!is.null(vif_vals)) {
      if (is.matrix(vif_vals)) {
        vif_out[rownames(vif_vals)] <- vif_vals[, 1]
      } else {
        vif_out[names(vif_vals)] <- as.numeric(vif_vals)
      }
    }
  }
  
  data.frame(Feature = names(vif_out), VIF = as.numeric(vif_out),
             Aliased = names(vif_out) %in% aliased_vars) %>%
    dplyr::arrange(dplyr::desc(VIF))
}

pca_analysis <- function(df, predictors, variance_threshold = 0.95) {
  X <- df[, predictors, drop = FALSE]
  X <- X[complete.cases(X), , drop = FALSE]
  X_scaled <- scale(X)
  pca <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
  var_exp <- pca$sdev^2 / sum(pca$sdev^2)
  cum_var <- cumsum(var_exp)
  list(pca = pca, var_explained = var_exp, cum_var = cum_var,
       n_components_needed = which(cum_var >= variance_threshold)[1])
}

#############################################################################
# 8. Min-Max normalization + calendar-based test split
#############################################################################
min_max_normalize <- function(x, min_x = min(x, na.rm = TRUE), max_x = max(x, na.rm = TRUE)) {
  (x - min_x) / (max_x - min_x)
}

normalize_dataset <- function(df, cols) {
  params <- list()
  for (cn in cols) {
    mn <- min(df[[cn]], na.rm = TRUE)
    mx <- max(df[[cn]], na.rm = TRUE)
    df[[cn]] <- min_max_normalize(df[[cn]], mn, mx)
    params[[cn]] <- c(min = mn, max = mx)
  }
  list(data = df, params = params)
}

split_dataset <- function(df, test_start_date = as.Date("2025-10-01"),
                          train_frac = 0.7, val_frac = 0.3) {
  df <- df %>% arrange(Date)
  test <- df[df$Date >= test_start_date, ]
  rest <- df[df$Date < test_start_date, ]
  
  n_rest <- nrow(rest)
  n_train <- floor(train_frac * n_rest)
  train <- rest[seq_len(n_train), ]
  valid <- rest[(n_train + 1):n_rest, ]
  
  list(train = train, valid = valid, test = test)
}

apply_norm_params <- function(df, params) {
  for (cn in names(params)) {
    df[[cn]] <- min_max_normalize(df[[cn]], params[[cn]]["min"], params[[cn]]["max"])
  }
  df
}

#############################################################################
# 9. Multi-Strategy Integration model for feature selection
#############################################################################
# Deterministic stability selection with joint optimization
compute_fitness_ms <- function(bin_vec, X, y, feature_names, base_model = "huber", k = 1.345,
                                lambda2 = 0.1, lambda3 = 1, lambda4 = 0.1,
                                vif_threshold = 10, corr_threshold = 0.9, cv_folds = 3,
                                weights = NULL, groups = NULL) {
  S_idx <- which(bin_vec == 1)
  if (length(S_idx) == 0) return(Inf)   # lambda1 * I(Empty) -> infinite penalty

  p <- length(feature_names)
  sel_names <- feature_names[S_idx]
  Xs <- X[, sel_names, drop = FALSE]
  # Carry weights AND group (coin) identity alongside X/y through the same
  # complete.cases() filter so both stay row-aligned after any rows drop.
  dat <- data.frame(y = y, Xs,
                     .w = if (is.null(weights)) rep(1, length(y)) else weights,
                     .g = if (is.null(groups)) rep("ALL", length(y)) else as.character(groups))
  dat <- dat[complete.cases(dat), ]
  n <- nrow(dat)
  if (n < (ncol(dat) + 5) || n < cv_folds * 3) return(Inf)

  # `k` is Huber's psi.huber tuning constant -- when base_model == "huber",
  # this is JOINTLY evolved alongside the feature-selection bits (see the
  # extended chromosome in multi_strategy_feature_selection() below), so
  # the fitness landscape genuinely reflects (feature subset, hyperparameter)
  # pairs rather than features selected under a single fixed hyperparameter.
  #
  # NOTE on why per-ROW case weights (`.w`) alone are insufficient: they
  # were originally added to correct for unequal ROW COUNTS per coin, but
  # that is not the actual failure mode here (all four Group 1 coins cover
  # the same date range, so row counts are already nearly equal). The real
  # issue is that Huber's OWN robustness mechanism downweights large
  # residuals DURING FITTING based on residual size, not row count or
  # external weight -- a small, jointly-optimized k makes this far more
  # aggressive, and if one coin's true feature-Close relationship differs
  # from the others', its residuals under a shared coefficient set run
  # systematically larger and get silently discounted by Huber's fitting
  # process itself, regardless of `.w`. External weights cannot fix an
  # internal downweighting mechanism. Fixed below (see MSE_cv) by scoring
  # per-coin fit quality directly.
  fit_predict_ms <- function(train_dat, test_dat) {
    fit <- tryCatch(
      MASS::rlm(y ~ . - .w - .g, data = train_dat, psi = MASS::psi.huber, k = k, maxit = 50, weights = train_dat$.w),
      error = function(e) NULL)
    if (is.null(fit)) return(rep(NA_real_, nrow(test_dat)))
    tryCatch(predict(fit, newdata = test_dat), error = function(e) rep(NA_real_, nrow(test_dat)))
  }

  ## MSE_cv: k-fold cross-validated MSE, averaged PER-COIN (not pooled
  ## across rows). This is the actual fix: within each fold, compute each
  ## coin's own (weighted) MSE, then average ACROSS COINS with equal
  ## weight. A (features, k) combination that fits BTC/ETH well but
  ## BCH/DASH badly now shows up directly as a high score -- it can no
  ## longer hide behind a pooled average where a systematically-underfit
  ## coin's contribution gets diluted (or, with a small k, actively
  ## discounted by Huber's own downweighting during fitting).
  folds <- sample(rep(seq_len(cv_folds), length.out = n))
  mse_folds <- sapply(seq_len(cv_folds), function(fold_k) {
    train_fold <- dat[folds != fold_k, , drop = FALSE]
    test_fold  <- dat[folds == fold_k, , drop = FALSE]
    pred <- fit_predict_ms(train_fold, test_fold)
    sq_err <- (test_fold$y - pred)^2
    per_group_mse <- tapply(seq_len(nrow(test_fold)), test_fold$.g, function(idx) {
      stats::weighted.mean(sq_err[idx], w = test_fold$.w[idx], na.rm = TRUE)
    })
    mean(unlist(per_group_mse), na.rm = TRUE)   # equal weight per coin, not per row
  })
  mse_cv <- mean(mse_folds, na.rm = TRUE)
  if (!is.finite(mse_cv)) return(Inf)

  ## Degenerate-fit guard: several engineered "price-based" features
  ## (PriceChange = Close - Open, CO_Ratio = Close/Open, TypicalPrice =
  ## mean(High,Low,Close), CLV, Return, LogReturn) are computed from the
  ## SAME ROW's Close value. If the search selects one of these alongside
  ## its algebraic "partner" raw feature(s) (e.g. Open + PriceChange), a
  ## linear model can reconstruct Close EXACTLY -- not genuine skill, an
  ## algebraic identity (x + y = z when x, y are both already given). This
  ## is invisible to the VIF/correlation penalties above, which only check
  ## predictor-predictor relationships, never whether a combination
  ## reconstructs the RESPONSE. Real financial Close prices are noisy and
  ## non-stationary (recall the ADF/skewness/kurtosis diagnostics
  ## throughout this analysis), so a near-zero MSE_cv is diagnostic of
  ## leakage, not accuracy -- reject it the same way an empty feature set
  ## is rejected, forcing the search toward the next-best genuine solution.
  y_var <- var(dat$y, na.rm = TRUE)
  if (is.finite(y_var) && y_var > 0 && (mse_cv / y_var) < 1e-6) return(Inf)

  ## Penalty_VIF (fit an OLS on the full subset purely for VIF diagnostics,
  ## regardless of base_model, since car::vif() needs an lm object)
  vif_penalty <- 0
  if (length(S_idx) > 1) {
    vif_vals <- tryCatch({ f2 <- lm(y ~ . - .w - .g, data = dat); car::vif(f2) }, error = function(e) NULL)
    if (!is.null(vif_vals)) {
      max_vif <- max(vif_vals)
      if (max_vif > vif_threshold) vif_penalty <- max_vif - vif_threshold
    }
  }

  ## Penalty_Correlation
  corr_penalty <- 0
  if (length(S_idx) > 1) {
    cmat <- suppressWarnings(cor(Xs, use = "pairwise.complete.obs")); diag(cmat) <- 0
    viol <- abs(cmat) - corr_threshold; viol[viol < 0] <- 0
    corr_penalty <- sum(viol) / 2
  }

  ## Complexity/parsimony penalty (|S| / p): discourages accumulating large
  ## feature sets purely because they marginally lower CV-MSE, which is
  ## what let the feature count grow to ~30 without this term -- a larger
  ## pooled linear model is more prone to overfitting toward whichever
  ## coin(s) dominate the pooled training rows, which showed up as
  ## systematic per-coin bias (BCH/DASH underprediction).
  complexity_penalty <- length(S_idx) / p

  mse_cv + lambda2 * vif_penalty + lambda3 * corr_penalty + lambda4 * complexity_penalty
}

# Population Hamming diversity over the BINARY feature-selection genes only
# (columns 1:p) -- the continuous hyperparameter gene (column p+1, if
# present) is excluded, since Hamming distance isn't meaningful for a
# continuous value.
population_hamming_diversity <- function(pop_bits) {
  n <- nrow(pop_bits)
  if (n < 2) return(0)
  d <- as.matrix(dist(pop_bits, method = "manhattan"))   # Hamming distance for 0/1 vectors
  mean(d[upper.tri(d)]) / ncol(pop_bits)
}
#############################################################################
# Multi-Strategy Integration model for feature selection (RESTORED).
# This is the genuine stochastic algorithm -- Opposition-Based Learning
# (OBL) initialization, Genetic Algorithm operators (tournament selection,
# crossover, mutation), a Guided Search Mechanism (GSM), and Quantum-
# Inspired Diversification (QID) -- with Huber's tuning constant k jointly
# optimized alongside the feature-selection bits in a single extended
# chromosome. This was MISSING from this file entirely (called in three
# places -- Section 9's main universal selection, asset_specific_quantum_
# model(), and stability_selection() -- but never defined, which would
# have caused an immediate "could not find function" error the moment any
# of those ran) and has been restored here, replacing the deterministic
# correlation-screen + greedy-hillclimb reimplementation
# (multi_strategy_feature_selection_deterministic_joint) that Sections 13
# and 21's stability comparison were calling instead. That deterministic
# procedure has no population, no crossover/mutation, no OBL, no GSM, and
# no QID -- its only randomness comes from the outer bootstrap resample,
# so when the dominant correlated features survive resampling (which they
# reliably do, since ~63%-overlapping bootstrap resamples barely perturb
# correlation rankings), it deterministically reconverges to the same
# answer every replicate -- producing an artificial Jaccard stability
# ceiling of 1.0 that reflects an absence of genuine algorithmic
# stochasticity, not a real property of a robust selection procedure.
#############################################################################
multi_strategy_feature_selection <- function(X, y, feature_names,
                                              pop_size = 20, max_iter = 30,
                                              base_model = c("huber", "ols", "rf"),
                                              lambda2 = 0.1, lambda3 = 1, lambda4 = 0.1,
                                              vif_threshold = 10, corr_threshold = 0.9, cv_folds = 3,
                                              pc = 0.8, pm = NULL,
                                              guide_top_frac = 0.10, guide_rate = 0.15,
                                              diversity_delta = 0.15, qid_theta = 0.10 * pi, qid_frac = 0.3,
                                              k_bounds = c(0.5, 3.0), k_mutation_sd = 0.15,
                                              weights = NULL, groups = NULL, seed = 42) {
  set.seed(seed)
  base_model <- match.arg(base_model)
  p <- length(feature_names)
  if (is.null(pm)) pm <- 1 / p
  k_col <- p + 1   # index of the continuous Huber-k gene in the chromosome

  fitness_of <- function(chromosome) {
    compute_fitness_ms(chromosome[1:p], X, y, feature_names, base_model = base_model, k = chromosome[k_col],
                        lambda2 = lambda2, lambda3 = lambda3, lambda4 = lambda4,
                        vif_threshold = vif_threshold, corr_threshold = corr_threshold, cv_folds = cv_folds,
                        weights = weights, groups = groups)
  }
  clip_k <- function(k) min(max(k, k_bounds[1]), k_bounds[2])

  ## ---- Phase 1: Initialization via Opposition-Based Learning ----
  ## Binary feature genes: random Bernoulli(0.3), as before. Continuous k
  ## gene: uniform over k_bounds. Opposition mirrors both: bitwise-complement
  ## for the binary genes, and reflection about the bounds' midpoint
  ## (k_bar = k_min + k_max - k) for the continuous gene -- the standard
  ## OBL formula for continuous variables.
  pop <- matrix(NA_real_, nrow = pop_size, ncol = p + 1)
  pop[, 1:p] <- matrix(rbinom(pop_size * p, 1, 0.3), nrow = pop_size, ncol = p)
  for (i in seq_len(pop_size)) if (sum(pop[i, 1:p]) == 0) pop[i, sample(p, 1)] <- 1
  pop[, k_col] <- runif(pop_size, k_bounds[1], k_bounds[2])

  opp <- pop
  opp[, 1:p] <- 1 - pop[, 1:p]
  opp[, k_col] <- k_bounds[1] + k_bounds[2] - pop[, k_col]

  full_pop <- rbind(pop, opp)                       # 2N candidates
  full_fitness <- apply(full_pop, 1, fitness_of)     # evaluate all 2N
  keep_idx <- order(full_fitness)[seq_len(pop_size)] # retain top N
  pop <- full_pop[keep_idx, , drop = FALSE]
  fitness_vals <- full_fitness[keep_idx]

  best_idx <- which.min(fitness_vals)
  best_solution <- pop[best_idx, ]; best_fitness <- fitness_vals[best_idx]
  history <- numeric(max_iter)

  tournament_select <- function() {
    idx <- sample(seq_len(pop_size), min(3, pop_size))
    idx[which.min(fitness_vals[idx])]
  }

  for (gen in seq_len(max_iter)) {
    ## ---- Phase 4 (prep): priority vector from top 10% (bits + k) ----
    K <- max(1, round(guide_top_frac * pop_size))
    top_idx <- order(fitness_vals)[seq_len(K)]
    p_guide <- colMeans(pop[top_idx, 1:p, drop = FALSE])
    k_guide <- mean(pop[top_idx, k_col])

    new_pop <- matrix(NA_real_, nrow = pop_size, ncol = p + 1)
    new_pop[1, ] <- best_solution   # elitism

    i <- 2
    while (i <= pop_size) {
      parent1 <- pop[tournament_select(), ]
      parent2 <- pop[tournament_select(), ]

      ## -- Phase 3: uniform crossover (bits) + blend crossover (k) --
      if (runif(1) < pc) {
        mask <- runif(p) < 0.5
        child_bits <- ifelse(mask, parent1[1:p], parent2[1:p])
        child_k <- runif(1, min(parent1[k_col], parent2[k_col]), max(parent1[k_col], parent2[k_col]))
      } else {
        child_bits <- parent1[1:p]; child_k <- parent1[k_col]
      }

      ## -- Phase 3: standard bit-flip mutation (bits) + Gaussian mutation (k) --
      flip_idx <- runif(p) < pm
      child_bits[flip_idx] <- 1 - child_bits[flip_idx]
      if (runif(1) < pm) child_k <- clip_k(child_k + rnorm(1, 0, k_mutation_sd))

      ## -- Phase 4: Guided Search Mechanism -- bits flip toward p_guide
      ## (probabilistically); k pulled toward the elite mean k_guide
      guide_idx <- runif(p) < guide_rate
      if (any(guide_idx)) child_bits[guide_idx] <- rbinom(sum(guide_idx), 1, p_guide[guide_idx])
      if (runif(1) < guide_rate) child_k <- clip_k(k_guide + rnorm(1, 0, k_mutation_sd * 0.5))

      if (sum(child_bits) == 0) child_bits[sample(p, 1)] <- 1   # avoid empty chromosome
      new_pop[i, ] <- c(child_bits, child_k)
      i <- i + 1
    }

    pop <- new_pop
    fitness_vals <- apply(pop, 1, fitness_of)

    cur_best <- which.min(fitness_vals)
    if (fitness_vals[cur_best] < best_fitness) {
      best_fitness <- fitness_vals[cur_best]
      best_solution <- pop[cur_best, ]
    }

    ## ---- Phase 5: diversity guardrail + Quantum-Inspired Diversification ----
    diversity <- population_hamming_diversity(pop[, 1:p, drop = FALSE])
    if (diversity < diversity_delta) {
      n_targets <- max(1, round(qid_frac * pop_size))
      elite_row <- which.min(fitness_vals)
      candidate_rows <- setdiff(seq_len(pop_size), elite_row)
      target_idx <- sample(candidate_rows, min(n_targets, length(candidate_rows)))

      for (ti in target_idx) {
        for (j in seq_len(p)) {
          # Only rotate genes currently matching the (over-represented) best
          # solution -- push their qubit amplitude AWAY from it via a
          # quantum rotation gate U(delta_theta), injecting STRUCTURED
          # diversity rather than resetting to pure random noise.
          if (pop[ti, j] == best_solution[j]) {
            delta_theta <- ifelse(best_solution[j] == 1, -qid_theta, qid_theta)
            angle <- asin(sqrt(0.5)) + delta_theta
            prob <- min(max(sin(angle)^2, 0.01), 0.99)
            pop[ti, j] <- rbinom(1, 1, prob)
          }
        }
        # Continuous analogue for the k gene: if this individual's k has
        # converged close to the (over-represented) best k, push it away by
        # a structured step rather than resetting to pure random noise.
        if (abs(pop[ti, k_col] - best_solution[k_col]) < k_mutation_sd) {
          step <- if (pop[ti, k_col] >= best_solution[k_col]) k_mutation_sd else -k_mutation_sd
          pop[ti, k_col] <- clip_k(pop[ti, k_col] + step)
        }
      }
      fitness_vals <- apply(pop, 1, fitness_of)
      cur_best <- which.min(fitness_vals)
      if (fitness_vals[cur_best] < best_fitness) {
        best_fitness <- fitness_vals[cur_best]
        best_solution <- pop[cur_best, ]
      }
      message("  Gen ", gen, ": diversity ", round(diversity, 3), " < ", diversity_delta,
              " -> QID triggered on ", length(target_idx), " individuals")
    }

    history[gen] <- best_fitness
  }

  list(selected_features = feature_names[best_solution[1:p] == 1],
       best_k = best_solution[k_col],
       best_fitness = best_fitness, history = history, active_features = feature_names)
}

# Runs multi_strategy_feature_selection() N_RUNS times, independently, on
# the SAME (unresampled) data with different random seeds, and summarizes
# convergence across runs. This is a distinct diagnostic from bootstrap
# stability (Section 13/21), which holds the algorithm's randomness fixed
# per replicate and varies the DATA; here the data is held fixed and the
# algorithm's OWN stochasticity (population initialization, tournament
# draws, crossover/mutation, QID) is what varies. A stochastic metaheuristic
# reported from a single run has no basis for a convergence claim -- the
# standard remedy is exactly this: run it many times and report how much
# the outcome varies. The single run with lowest fitness is retained as the
# "official" reported model (consistent with prior single-run behavior),
# but is now reported ALONGSIDE the full distribution across runs so a
# reader can judge whether it is a representative outcome or an outlier.
run_multi_strategy_multiple <- function(X, y, feature_names, n_runs = 20, base_seed = 42, ...) {
  runs <- vector("list", n_runs)
  for (i in seq_len(n_runs)) {
    runs[[i]] <- multi_strategy_feature_selection(X = X, y = y, feature_names = feature_names,
                                                   seed = base_seed + i, ...)
    message("  Run ", i, "/", n_runs, ": best_fitness = ", round(runs[[i]]$best_fitness, 6),
            ", k = ", round(runs[[i]]$best_k, 4), ", n_features = ", length(runs[[i]]$selected_features))
  }

  fitness_vals <- sapply(runs, `[[`, "best_fitness")
  k_vals <- sapply(runs, `[[`, "best_k")
  n_features_vals <- sapply(runs, function(r) length(r$selected_features))
  best_run_idx <- which.min(fitness_vals)

  # Per-feature selection frequency across the N independent runs (NOT a
  # bootstrap resampling frequency -- same data every run, only the
  # algorithm's own randomness differs).
  sel_matrix <- matrix(0, nrow = n_runs, ncol = length(feature_names), dimnames = list(NULL, feature_names))
  for (i in seq_len(n_runs)) sel_matrix[i, runs[[i]]$selected_features] <- 1
  selection_freq <- colMeans(sel_matrix)

  # Coefficient of variation of best_fitness across runs is the primary
  # convergence diagnostic: a low CV means independent runs reliably find
  # solutions of similar quality (evidence the search has converged to
  # essentially the same region of the fitness landscape); a high CV means
  # the reported single-run result could easily have been a lucky or
  # unlucky draw rather than a representative outcome.
  fitness_cv <- sd(fitness_vals) / mean(fitness_vals)

  summary_table <- data.frame(
    Run = seq_len(n_runs), Best_Fitness = fitness_vals, Best_K = k_vals, N_Features = n_features_vals,
    Is_Best_Run = seq_len(n_runs) == best_run_idx
  )
  selection_freq_table <- data.frame(Feature = names(selection_freq), Selection_Frequency = as.numeric(selection_freq)) %>%
    arrange(desc(Selection_Frequency))

  message("Convergence diagnostics across ", n_runs, " independent runs:")
  message("  Best fitness: mean = ", round(mean(fitness_vals), 6), ", SD = ", round(sd(fitness_vals), 6),
          ", CV = ", round(fitness_cv, 4), " (lower CV = more consistent convergence)")
  message("  Best run (lowest fitness) is run ", best_run_idx, " with fitness = ", round(min(fitness_vals), 6))

  list(runs = runs, best_run = runs[[best_run_idx]], best_run_idx = best_run_idx,
       summary_table = summary_table, selection_freq_table = selection_freq_table,
       fitness_mean = mean(fitness_vals), fitness_sd = sd(fitness_vals), fitness_cv = fitness_cv)
}

# One-at-a-time (OAT) sensitivity analysis: varies each of lambda2 (VIF
# penalty), lambda3 (correlation penalty), lambda4 (complexity penalty),
# pop_size, and max_iter individually across a small grid around the
# baseline defaults used throughout this pipeline, holding all other
# parameters fixed at baseline, and runs the search ONCE per configuration
# (single runs, not the N_RUNS_MULTISTRATEGY multi-run diagnostic above --
# a full factorial or multi-run-per-configuration sweep here would be
# computationally prohibitive; OAT is the standard, tractable compromise
# for a search this expensive per evaluation). Reports how the number of
# selected features, best fitness, and best k each respond to varying one
# hyperparameter at a time, so a reader can judge whether the reported
# results are a fragile artifact of one specific configuration or robust
# across a reasonable neighborhood of it.
sensitivity_analysis_multi_strategy <- function(X, y, feature_names, weights = NULL, groups = NULL,
                                                 baseline = list(lambda2 = 0.1, lambda3 = 1, lambda4 = 0.1,
                                                                 pop_size = 20, max_iter = 30),
                                                 grids = list(lambda2 = c(0.05, 0.1, 0.2),
                                                              lambda3 = c(0.5, 1, 2),
                                                              lambda4 = c(0.05, 0.1, 0.2),
                                                              pop_size = c(10, 20, 30),
                                                              max_iter = c(15, 30, 45)),
                                                 seed = 42) {
  results <- list()
  for (param_name in names(grids)) {
    for (val in grids[[param_name]]) {
      cfg <- baseline
      cfg[[param_name]] <- val
      message("Sensitivity: varying ", param_name, " = ", val, " (others at baseline)...")
      res <- do.call(multi_strategy_feature_selection,
                      c(list(X = X, y = y, feature_names = feature_names, weights = weights, groups = groups,
                             seed = seed), cfg))
      results[[length(results) + 1]] <- data.frame(
        Varied_Parameter = param_name, Value = val,
        Lambda2 = cfg$lambda2, Lambda3 = cfg$lambda3, Lambda4 = cfg$lambda4,
        Pop_Size = cfg$pop_size, Max_Iter = cfg$max_iter,
        N_Features_Selected = length(res$selected_features),
        Best_Fitness = res$best_fitness, Best_K = res$best_k,
        Selected_Features = paste(res$selected_features, collapse = "; "),
        Is_Baseline_Value = identical(val, baseline[[param_name]])
      )
    }
  }
  bind_rows(results)
}

#==========================================

#############################################################################
# 10. Residual diagnostics
#############################################################################
residual_diagnostics <- function(fit, data, lb_lag = 10) {
  res <- residuals(fit)
  fitted_vals <- fitted(fit)
  X_design <- model.matrix(fit)
  H <- diag(X_design %*% solve(t(X_design) %*% X_design) %*% t(X_design))
  sigma_hat <- summary(fit)$sigma
  r_star <- res / (sigma_hat * sqrt(1 - H))
  
  p_params <- ncol(X_design)
  cooks_d <- (r_star^2 / p_params) * (H / (1 - H))
  
  res_sw <- if (length(res) > 5000) sample(res, 5000) else res
  
  shapiro_res <- tryCatch(shapiro.test(res_sw), error = function(e) NULL)
  jb_res <- tryCatch(tseries::jarque.bera.test(res), error = function(e) NULL)
  bp_res <- tryCatch(lmtest::bptest(lm(res ~ fitted_vals)), error = function(e) NULL)
  white_res <- tryCatch(lmtest::bptest(lm(res ~ fitted_vals), varformula = ~ fitted_vals + I(fitted_vals^2)),
                        error = function(e) NULL)
  lb_res <- tryCatch(Box.test(res, lag = lb_lag, type = "Ljung-Box"), error = function(e) NULL)
  dw_res <- tryCatch(lmtest::dwtest(fit), error = function(e) NULL)
  
  list(standardized_residuals = r_star, cooks_distance = cooks_d, leverage = H,
       shapiro = shapiro_res, jarque_bera = jb_res,
       breusch_pagan = bp_res, white_test = white_res,
       ljung_box = lb_res, durbin_watson = dw_res,
       outlier_idx = which(abs(r_star) > 3))
}

residual_diagnostics_table <- function(diag_res) {
  tests <- list(
    "Normality: Shapiro-Wilk" = diag_res$shapiro,
    "Normality: Jarque-Bera" = diag_res$jarque_bera,
    "Homoscedasticity: Breusch-Pagan" = diag_res$breusch_pagan,
    "Homoscedasticity: White" = diag_res$white_test,
    "Independence: Ljung-Box" = diag_res$ljung_box,
    "Independence: Durbin-Watson" = diag_res$durbin_watson
  )
  tbl <- bind_rows(lapply(names(tests), function(nm) {
    t <- tests[[nm]]
    if (is.null(t)) return(data.frame(Test = nm, Statistic = NA_real_, p_value = NA_real_))
    stat <- if (!is.null(t$statistic)) unname(t$statistic) else NA_real_
    pval <- if (!is.null(t$p.value)) t$p.value else NA_real_
    data.frame(Test = nm, Statistic = stat, p_value = pval)
  }))

  # Multiple testing correction: 6 diagnostic tests are evaluated on the
  # same fitted model, so judging each individually against the nominal
  # alpha=0.05 inflates the family-wise Type I error rate above 5% (up to
  # ~1-(0.95)^6 =~ 26% under independence). Two corrections are reported:
  # Bonferroni (conservative; controls family-wise error rate exactly under
  # independence, alpha/6 per test) and Benjamini-Hochberg (1995) (controls
  # the false discovery rate, less conservative and more standard when
  # tests examine related but distinct aspects of adequacy -- here,
  # normality, homoscedasticity, and independence -- rather than one tight
  # family of highly correlated hypotheses). Both raw and adjusted
  # conclusions are kept so a reader can see whether any decision actually
  # changes once correction is applied.
  tbl$p_value_bonferroni <- p.adjust(tbl$p_value, method = "bonferroni")
  tbl$p_value_BH <- p.adjust(tbl$p_value, method = "BH")
  concl <- function(p, label) ifelse(is.na(p), "n/a",
                                      ifelse(p < 0.05, paste0("reject H0 (", label, " p < 0.05)"), "fail to reject H0"))
  tbl$Conclusion_raw <- concl(tbl$p_value, "raw")
  tbl$Conclusion_bonferroni <- concl(tbl$p_value_bonferroni, "Bonferroni-adjusted")
  tbl$Conclusion_BH <- concl(tbl$p_value_BH, "BH-adjusted")
  tbl
}

plot_residual_diagnostics <- function(fit, data, response = "Close") {
  res <- residuals(fit)
  fitted_vals <- fitted(fit)
  diag_res <- residual_diagnostics(fit, data)
  r_star <- diag_res$standardized_residuals
  cooks_d <- diag_res$cooks_distance
  n <- length(res)
  
  par(mfrow = c(2, 2))
  
  qqnorm(res, main = "Q-Q Plot of Residuals")
  qqline(res, col = "red")
  
  plot(fitted_vals, res, xlab = "Fitted values", ylab = "Residuals",
       main = "Residuals vs Fitted")
  abline(h = 0, col = "red")
  lines(lowess(fitted_vals, res), col = "blue", lwd = 2)
  
  sqrt_r <- sqrt(abs(r_star))
  plot(fitted_vals, sqrt_r, xlab = "Fitted values",
       ylab = expression(sqrt("|Standardized Residuals|")), main = "Scale-Location")
  abline(h = mean(sqrt_r, na.rm = TRUE), col = "red", lty = 2)
  lines(lowess(fitted_vals, sqrt_r), col = "blue", lwd = 2)
  legend("topright", legend = c("Mean (ideal)", "Lowess trend"), col = c("red", "blue"),
         lty = c(2, 1), lwd = c(1, 2), bty = "n", cex = 0.7)
  
  cook_thr <- if (n < 500) 4 / n else stats::quantile(cooks_d, 0.99, na.rm = TRUE)
  influential <- cooks_d > cook_thr
  pt_col <- ifelse(influential, "red", adjustcolor("black", 0.4))
  size_ratio <- cooks_d / max(median(cooks_d, na.rm = TRUE), 1e-12)
  pt_cex <- pmin(0.4 + 0.15 * sqrt(size_ratio), 2.5)
  plot(data[[response]], res, xlab = response, ylab = "Residuals",
       main = "Residuals vs Close (size/color ~ Cook's D)",
       col = pt_col, cex = pt_cex, pch = 19)
  abline(h = 0, col = "gray40", lty = 2)
  legend("topright", legend = c(paste0("Top 1% Cook's D (>", signif(cook_thr, 3), ")"), "Remaining points"),
         col = c("red", "black"), pch = 19, bty = "n", cex = 0.7)
  
  par(mfrow = c(1, 1))
}

#############################################################################
# 12. Bootstrap uncertainty quantification
#############################################################################
residual_bootstrap <- function(fit, data, response = "Close", B = 200, seed = 123) {
  set.seed(seed)
  yhat <- fitted(fit)
  res <- residuals(fit)
  n <- length(yhat)
  coef_names <- names(coef(fit))
  boot_coefs <- matrix(NA, nrow = B, ncol = length(coef_names),
                       dimnames = list(NULL, coef_names))
  
  for (b in 1:B) {
    res_star <- sample(res, n, replace = TRUE)
    dat_star <- data
    dat_star[[response]] <- yhat + res_star
    fit_star <- tryCatch(MASS::rlm(formula(fit), data = dat_star,
                                   psi = MASS::psi.huber, maxit = 50),
                         error = function(e) NULL)
    if (!is.null(fit_star)) boot_coefs[b, ] <- coef(fit_star)
  }
  ci <- apply(boot_coefs, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  list(boot_coefs = boot_coefs, ci_95 = t(ci))
}

#############################################################################
# 13. Feature-selection stability analysis
#############################################################################
stability_selection <- function(X, y, feature_names, B = 30, pi_thr = 0.6, qfs_args = list()) {
  n <- nrow(X)
  p <- length(feature_names)
  selection_matrix <- matrix(0, nrow = B, ncol = p, dimnames = list(NULL, feature_names))
  selections <- list()
  
  for (b in 1:B) {
    # Meinshausen & Bühlmann (2010)'s actual prescribed procedure is
    # subsampling WITHOUT replacement at size floor(n/2), not a standard
    # full-size WITH-replacement bootstrap. This matters: any two standard
    # bootstrap resamples overlap by ~63% on average, which understates
    # genuine selection instability for every method being compared here --
    # a deterministic or near-deterministic procedure (e.g. AIC/BIC
    # stepwise) can look artificially stable simply because the data barely
    # changes between replicates. True subsampling perturbs the training
    # data more, giving a fairer, more discriminating comparison across all
    # six methods in Section 21.
    idx <- sample(1:n, size = floor(n / 2), replace = FALSE)
    qfs_args_b <- utils::modifyList(qfs_args, list(seed = 1000 + b))
    res <- do.call(multi_strategy_feature_selection,
                   c(list(X = X[idx, , drop = FALSE], y = y[idx], feature_names = feature_names), qfs_args_b))
    selections[[b]] <- res$selected_features
    selection_matrix[b, res$selected_features] <- 1
  }
  
  pi_hat <- colMeans(selection_matrix)
  jaccard <- function(s1, s2) {
    if (length(s1) == 0 && length(s2) == 0) return(1)
    length(intersect(s1, s2)) / length(union(s1, s2))
  }
  jvals <- c()
  if (B > 1) for (i in 1:(B - 1)) for (j in (i + 1):B) jvals <- c(jvals, jaccard(selections[[i]], selections[[j]]))

  # Standard error / CI for avg_jaccard_stability via bootstrap-of-REPLICATES
  # (not bootstrap-of-PAIRS): the B*(B-1)/2 pairwise Jaccard values are not
  # independent observations -- each of the B selections contributes to
  # B-1 overlapping pairs, so naively computing sd(jvals)/sqrt(length(jvals))
  # would understate the true uncertainty. Resampling at the level of the
  # actual unit of replication (the B selections themselves) and
  # recomputing the mean pairwise Jaccard among each resampled set gives a
  # statistically valid standard error, and is essentially free
  # computationally since it only reshuffles the already-computed
  # `selections` list rather than re-running the expensive search.
  n_boot_se <- 500
  boot_means <- replicate(n_boot_se, {
    boot_idx <- sample(seq_len(B), B, replace = TRUE)
    boot_sel <- selections[boot_idx]
    bvals <- c()
    if (B > 1) for (i in 1:(B - 1)) for (j in (i + 1):B) bvals <- c(bvals, jaccard(boot_sel[[i]], boot_sel[[j]]))
    if (length(bvals) == 0) NA_real_ else mean(bvals)
  })
  jaccard_se <- sd(boot_means, na.rm = TRUE)
  jaccard_ci <- quantile(boot_means, c(0.025, 0.975), na.rm = TRUE, type = 7)

  list(pi_hat = pi_hat, stable_features = names(pi_hat)[pi_hat >= pi_thr],
       avg_jaccard_stability = mean(jvals), jaccard_se = jaccard_se,
       jaccard_ci_lower = unname(jaccard_ci[1]), jaccard_ci_upper = unname(jaccard_ci[2]),
       selections = selections)
}


#=============================================
#---------------------------------------

#=========================================

#############################################################################
# 14. Cross-validation
#############################################################################
kfold_cv <- function(data, response, predictors, K = 10, seed = 1) {
  set.seed(seed)
  n <- nrow(data)
  folds <- sample(rep(1:K, length.out = n))
  form <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  bind_rows(lapply(1:K, function(k) {
    fit <- MASS::rlm(form, data = data[folds != k, ], psi = MASS::psi.huber, maxit = 50)
    pred <- predict(fit, newdata = data[folds == k, ])
    actual <- data[[response]][folds == k]
    data.frame(Fold = k, MAE = mean(abs(actual - pred)), RMSE = sqrt(mean((actual - pred)^2)),
               R2 = 1 - sum((actual - pred)^2) / sum((actual - mean(actual))^2))
  }))
}

expanding_window_cv <- function(data, response, predictors, initial_frac = 0.5, horizon = 30, step = 30) {
  data <- data %>% arrange(Date)
  n <- nrow(data)
  form <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  t0 <- floor(initial_frac * n)
  results <- list()
  i <- 1
  t <- t0
  while ((t + horizon) <= n) {
    fit <- MASS::rlm(form, data = data[1:t, ], psi = MASS::psi.huber, maxit = 50)
    test_d <- data[(t + 1):(t + horizon), ]
    pred <- predict(fit, newdata = test_d)
    actual <- test_d[[response]]
    results[[i]] <- data.frame(Window = i, TrainEnd = t, MAE = mean(abs(actual - pred)),
                               RMSE = sqrt(mean((actual - pred)^2)),
                               R2 = 1 - sum((actual - pred)^2) / sum((actual - mean(actual))^2))
    t <- t + step
    i <- i + 1
  }
  bind_rows(results)
}

#############################################################################
# 15. Out-of-sample validation + generalizability
#############################################################################
out_of_sample_validation <- function(fit, test_data, response) {
  pred <- predict(fit, newdata = test_data)
  actual <- test_data[[response]]
  data.frame(MAE = mean(abs(actual - pred)), RMSE = sqrt(mean((actual - pred)^2)),
             R2 = 1 - sum((actual - pred)^2) / sum((actual - mean(actual))^2))
}

generalizability_test <- function(fit, group2_valid_list, response, selected_features) {
  bind_rows(lapply(names(group2_valid_list), function(nm) {
    dat <- group2_valid_list[[nm]][, c(response, selected_features)]
    dat <- dat[complete.cases(dat), ]
    res <- out_of_sample_validation(fit, dat, response)
    res$Symbol <- nm
    res %>% relocate(Symbol)
  }))
}

#############################################################################
# 16. Model comparison framework
#############################################################################
fit_huber <- function(train_d, response, predictors, k = 1.345, maxit = 100,
                      max_drop_attempts = 5, weights = NULL) {
  active_predictors <- predictors
  attempt <- 0
  repeat {
    form <- as.formula(paste(response, "~", paste(active_predictors, collapse = " + ")))
    rlm_args <- list(formula = form, data = train_d, psi = MASS::psi.huber, k = k, maxit = maxit)
    if (!is.null(weights)) rlm_args$weights <- weights
    fit <- tryCatch(do.call(MASS::rlm, rlm_args), error = function(e) e)
    if (!inherits(fit, "error")) return(fit)
    
    is_singular <- grepl("singular", conditionMessage(fit), ignore.case = TRUE)
    attempt <- attempt + 1
    if (!is_singular || attempt > max_drop_attempts || length(active_predictors) <= 1) {
      stop(fit)
    }
    
    vif_vals <- tryCatch({
      md <- train_d[, c(response, active_predictors)]
      md <- md[complete.cases(md), ]
      f2 <- lm(as.formula(paste(response, "~", paste(active_predictors, collapse = "+"))), data = md)
      car::vif(f2)
    }, error = function(e) NULL)
    
    if (is.null(vif_vals) || length(vif_vals) == 0) {
      al <- tryCatch(alias(lm(as.formula(paste(response, "~", paste(active_predictors, collapse = "+"))),
                              data = train_d)), error = function(e) NULL)
      aliased <- if (!is.null(al) && !is.null(al$Complete)) rownames(al$Complete) else NULL
      drop_feat <- if (!is.null(aliased) && length(aliased) > 0) aliased[1] else active_predictors[length(active_predictors)]
    } else {
      drop_feat <- names(which.max(vif_vals))
    }
    
    message("  fit_huber: singular fit, dropping most collinear feature '", drop_feat,
            "' and retrying (", length(active_predictors) - 1, " features remain)")
    active_predictors <- setdiff(active_predictors, drop_feat)
  }
}

fit_xgboost <- function(train_d, response, predictors, nrounds = 300, max_depth = 4, eta = 0.05) {
  X <- as.matrix(train_d[, predictors])
  y <- train_d[[response]]
  dtrain <- xgboost::xgb.DMatrix(data = X, label = y)
  params <- list(max_depth = max_depth, eta = eta, objective = "reg:squarederror")
  xgboost::xgb.train(params = params, data = dtrain, nrounds = nrounds, verbose = 0)
}

predict_xgboost <- function(model, newdata, predictors) {
  X <- as.matrix(newdata[, predictors])
  predict(model, xgboost::xgb.DMatrix(data = X))
}

adaboost_r2 <- function(X, y, n_estimators = 50, max_depth = 4) {
  X <- as.data.frame(X)
  n <- nrow(X)
  w <- rep(1 / n, n)
  models <- list()
  betas <- c()
  dat <- data.frame(y = y, X)
  for (m in 1:n_estimators) {
    idx <- sample(1:n, n, replace = TRUE, prob = w)
    tree <- rpart::rpart(y ~ ., data = dat[idx, ], control = rpart::rpart.control(maxdepth = max_depth))
    pred <- predict(tree, newdata = dat)
    err <- abs(y - pred)
    Lmax <- max(err)
    if (Lmax == 0) Lmax <- 1e-8
    L <- err / Lmax
    avg_loss <- sum(w * L)
    if (avg_loss >= 0.5) break
    beta <- avg_loss / (1 - avg_loss)
    w <- w * beta^(1 - L)
    w <- w / sum(w)
    models[[length(models) + 1]] <- tree
    betas <- c(betas, beta)
  }
  list(models = models, betas = betas)
}

predict_adaboost_r2 <- function(ab_model, newdata) {
  newdata <- as.data.frame(newdata)
  preds <- sapply(ab_model$models, function(tr) predict(tr, newdata = newdata))
  if (is.null(dim(preds))) preds <- matrix(preds, ncol = length(ab_model$models))
  weights <- log(1 / ab_model$betas)
  weights[!is.finite(weights)] <- 0
  apply(preds, 1, function(row) {
    ord <- order(row)
    cw <- cumsum(weights[ord])
    med_idx <- which(cw >= 0.5 * sum(weights))[1]
    row[ord][med_idx]
  })
}

fit_extra_trees <- function(train_d, response, predictors, num.trees = 300, min.node.size = 5, mtry = NULL) {
  form <- as.formula(paste(response, "~", paste(predictors, collapse = "+")))
  if (is.null(mtry)) mtry <- max(1, floor(sqrt(length(predictors))))
  ranger::ranger(form, data = train_d, splitrule = "extratrees",
                 num.trees = num.trees, min.node.size = min.node.size, mtry = mtry)
}

fit_gbm <- function(train_d, response, predictors, n.trees = 300, interaction.depth = 3, shrinkage = 0.05) {
  gbm::gbm(as.formula(paste(response, "~", paste(predictors, collapse = "+"))),
           data = train_d, distribution = "gaussian", n.trees = n.trees,
           interaction.depth = interaction.depth, shrinkage = shrinkage, verbose = FALSE)
}

fit_rf <- function(train_d, response, predictors, ntree = 300, mtry = NULL, nodesize = 5) {
  form <- as.formula(paste(response, "~", paste(predictors, collapse = "+")))
  if (is.null(mtry)) {
    randomForest::randomForest(form, data = train_d, ntree = ntree, nodesize = nodesize)
  } else {
    randomForest::randomForest(form, data = train_d, ntree = ntree, mtry = mtry, nodesize = nodesize)
  }
}

fit_svm <- function(train_d, response, predictors, cost = 1, gamma = NULL) {
  form <- as.formula(paste(response, "~", paste(predictors, collapse = "+")))
  if (is.null(gamma)) {
    e1071::svm(form, data = train_d, type = "eps-regression", kernel = "radial", cost = cost)
  } else {
    e1071::svm(form, data = train_d, type = "eps-regression", kernel = "radial", cost = cost, gamma = gamma)
  }
}

predict_voting_ensemble <- function(gbm_m, rf_m, svm_m, newdata) {
  preds <- cbind(gbm = predict(gbm_m, newdata = newdata, n.trees = gbm_m$n.trees),
                 rf = predict(rf_m, newdata = newdata),
                 svm = predict(svm_m, newdata = newdata))
  rowMeans(preds)
}

fit_stacking_ensemble <- function(train_d, response, predictors,
                                  svm_params = list(), rf_params = list(), gbm_params = list()) {
  svm_m <- do.call(fit_svm, c(list(train_d = train_d, response = response, predictors = predictors), svm_params))
  rf_m <- do.call(fit_rf, c(list(train_d = train_d, response = response, predictors = predictors), rf_params))
  gbm_m <- do.call(fit_gbm, c(list(train_d = train_d, response = response, predictors = predictors), gbm_params))
  base_preds <- data.frame(svm = predict(svm_m, newdata = train_d),
                           rf = predict(rf_m, newdata = train_d),
                           gbm = predict(gbm_m, newdata = train_d, n.trees = gbm_m$n.trees),
                           y = train_d[[response]])
  meta_model <- lm(y ~ svm + rf + gbm, data = base_preds)
  list(svm = svm_m, rf = rf_m, gbm = gbm_m, meta = meta_model)
}

predict_stacking_ensemble <- function(stack, newdata) {
  base_preds <- data.frame(svm = predict(stack$svm, newdata = newdata),
                           rf = predict(stack$rf, newdata = newdata),
                           gbm = predict(stack$gbm, newdata = newdata, n.trees = stack$gbm$n.trees))
  predict(stack$meta, newdata = base_preds)
}

#############################################################################
# 16b. Hyperparameter tuning
#############################################################################
tune_all_models <- function(train_d, valid_d, response, predictors,
                            init_points = 4, iters_n = 6, fixed_huber_k = NULL) {
  message("Tuning hyperparameters via Bayesian optimization on Group 1 TRAIN -> VALIDATION...")
  actual_valid <- valid_d[[response]]
  p <- length(predictors)
  best <- list()
  
  eval_mae <- function(pred) {
    if (is.null(pred) || all(is.na(pred))) return(1e6)
    mean(abs(actual_valid - pred), na.rm = TRUE)
  }
  
  bayes_tune <- function(objective_fn, bounds, label) {
    result <- tryCatch(
      rBayesianOptimization::BayesianOptimization(
        FUN = objective_fn, bounds = bounds,
        init_points = max(init_points, length(bounds) + 2), n_iter = iters_n,
        acq = "ucb", kappa = 2.576, verbose = FALSE),
      error = function(e) {
        message("  ", label, ": Bayesian optimization failed (", conditionMessage(e),
                "); falling back to bounds midpoint")
        NULL
      })
    if (is.null(result)) return(setNames(lapply(bounds, mean), names(bounds)))
    setNames(lapply(names(bounds), function(nm) result$Best_Par[[nm]]), names(bounds))
  }
  
  if (!is.null(fixed_huber_k)) {
    best$Huber <- list(k = fixed_huber_k)
    message("  Huber: k = ", round(fixed_huber_k, 3), " (carried forward from joint feature-selection optimization)")
  } else {
    huber_obj <- function(k) {
      fit <- tryCatch(fit_huber(train_d, response, predictors, k = k), error = function(e) NULL)
      list(Score = -(if (is.null(fit)) 1e6 else eval_mae(predict(fit, newdata = valid_d))), Pred = 0)
    }
    huber_best <- bayes_tune(huber_obj, list(k = c(0.5, 3.0)), "Huber")
    best$Huber <- list(k = huber_best$k)
    message("  Huber: k = ", round(best$Huber$k, 3))
  }
  
  ada_obj <- function(n_estimators, max_depth) {
    ne <- round(n_estimators)
    md <- round(max_depth)
    fit <- tryCatch(adaboost_r2(train_d[, predictors, drop = FALSE], train_d[[response]],
                                n_estimators = ne, max_depth = md),
                    error = function(e) NULL)
    list(Score = -(if (is.null(fit)) 1e6 else eval_mae(predict_adaboost_r2(fit, valid_d[, predictors, drop = FALSE]))),
         Pred = 0)
  }
  ada_best <- bayes_tune(ada_obj, list(n_estimators = c(20, 100), max_depth = c(2, 8)), "AdaBoost")
  best$AdaBoost <- list(n_estimators = round(ada_best$n_estimators), max_depth = round(ada_best$max_depth))
  message("  AdaBoost: n_estimators = ", best$AdaBoost$n_estimators, ", max_depth = ", best$AdaBoost$max_depth)
  
  et_obj <- function(num_trees, min_node_size) {
    nt <- round(num_trees)
    mns <- round(min_node_size)
    fit <- tryCatch(fit_extra_trees(train_d, response, predictors, num.trees = nt, min.node.size = mns),
                    error = function(e) NULL)
    list(Score = -(if (is.null(fit)) 1e6 else eval_mae(predict(fit, data = valid_d)$predictions)), Pred = 0)
  }
  et_best <- bayes_tune(et_obj, list(num_trees = c(100, 600), min_node_size = c(2, 15)), "ExtraTrees")
  best$ExtraTrees <- list(num.trees = round(et_best$num_trees), min.node.size = round(et_best$min_node_size))
  message("  ExtraTrees: num.trees = ", best$ExtraTrees$num.trees, ", min.node.size = ", best$ExtraTrees$min.node.size)
  
  xgb_obj <- function(max_depth, eta) {
    md <- round(max_depth)
    fit <- tryCatch(fit_xgboost(train_d, response, predictors, max_depth = md, eta = eta), error = function(e) NULL)
    list(Score = -(if (is.null(fit)) 1e6 else eval_mae(predict_xgboost(fit, valid_d, predictors))), Pred = 0)
  }
  xgb_best <- bayes_tune(xgb_obj, list(max_depth = c(2, 8), eta = c(0.01, 0.2)), "XGBoost")
  best$XGBoost <- list(max_depth = round(xgb_best$max_depth), eta = xgb_best$eta)
  message("  XGBoost: max_depth = ", best$XGBoost$max_depth, ", eta = ", round(best$XGBoost$eta, 4))
  
  gbm_obj <- function(interaction_depth, shrinkage) {
    id <- round(interaction_depth)
    fit <- tryCatch(fit_gbm(train_d, response, predictors, interaction.depth = id, shrinkage = shrinkage),
                    error = function(e) NULL)
    list(Score = -(if (is.null(fit)) 1e6 else eval_mae(predict(fit, newdata = valid_d, n.trees = fit$n.trees))),
         Pred = 0)
  }
  gbm_best <- bayes_tune(gbm_obj, list(interaction_depth = c(1, 6), shrinkage = c(0.01, 0.2)), "GBM")
  best$GBM <- list(interaction.depth = round(gbm_best$interaction_depth), shrinkage = gbm_best$shrinkage)
  message("  GBM: interaction.depth = ", best$GBM$interaction.depth, ", shrinkage = ", round(best$GBM$shrinkage, 4))
  
  rf_obj <- function(mtry, nodesize) {
    mt <- max(1, round(mtry))
    ns <- round(nodesize)
    fit <- tryCatch(fit_rf(train_d, response, predictors, mtry = mt, nodesize = ns), error = function(e) NULL)
    list(Score = -(if (is.null(fit)) 1e6 else eval_mae(predict(fit, newdata = valid_d))), Pred = 0)
  }
  rf_best <- bayes_tune(rf_obj, list(mtry = c(1, max(2, p)), nodesize = c(2, 15)), "RandomForest")
  best$RandomForest <- list(mtry = max(1, round(rf_best$mtry)), nodesize = round(rf_best$nodesize))
  message("  RandomForest: mtry = ", best$RandomForest$mtry, ", nodesize = ", best$RandomForest$nodesize)
  
  svm_obj <- function(cost, gamma) {
    fit <- tryCatch(fit_svm(train_d, response, predictors, cost = cost, gamma = gamma), error = function(e) NULL)
    list(Score = -(if (is.null(fit)) 1e6 else eval_mae(predict(fit, newdata = valid_d))), Pred = 0)
  }
  svm_best <- bayes_tune(svm_obj, list(cost = c(0.01, 20), gamma = c(0.001, 1)), "SVM")
  best$SVM <- list(cost = svm_best$cost, gamma = svm_best$gamma)
  message("  SVM: cost = ", round(best$SVM$cost, 3), ", gamma = ", round(best$SVM$gamma, 5))
  
  best
}

hyperparams_to_table <- function(hyperparams) {
  bind_rows(lapply(names(hyperparams), function(m) {
    p <- hyperparams[[m]]
    data.frame(Model = m, Parameter = names(p), Value = as.character(unlist(p)))
  }))
}

compare_models <- function(train_d, test_d, response, predictors, hyperparams = NULL) {
  results <- list()
  actual <- test_d[[response]]
  hp <- function(name, defaults) {
    if (!is.null(hyperparams) && !is.null(hyperparams[[name]])) {
      modifyList(defaults, hyperparams[[name]])
    } else {
      defaults
    }
  }
  
  t0 <- Sys.time()
  huber_params <- hp("Huber", list(k = 1.345, maxit = 100))
  huber_fit <- do.call(fit_huber, c(list(train_d = train_d, response = response, predictors = predictors), huber_params))
  t1 <- Sys.time()
  results[["Huber"]] <- list(pred = predict(huber_fit, newdata = test_d),
                             time = as.numeric(difftime(t1, t0, units = "secs")),
                             mem = as.numeric(object.size(huber_fit)) / 1024)
  
  t0 <- Sys.time()
  ada_params <- hp("AdaBoost", list(n_estimators = 50, max_depth = 4))
  ab_fit <- do.call(adaboost_r2, c(list(X = train_d[, predictors, drop = FALSE], y = train_d[[response]]), ada_params))
  t1 <- Sys.time()
  results[["AdaBoost_R2"]] <- list(pred = predict_adaboost_r2(ab_fit, test_d[, predictors, drop = FALSE]),
                                   time = as.numeric(difftime(t1, t0, units = "secs")),
                                   mem = as.numeric(object.size(ab_fit)) / 1024)
  
  t0 <- Sys.time()
  et_params <- hp("ExtraTrees", list(num.trees = 300, min.node.size = 5))
  et_fit <- do.call(fit_extra_trees, c(list(train_d = train_d, response = response, predictors = predictors), et_params))
  t1 <- Sys.time()
  results[["ExtraTrees"]] <- list(pred = predict(et_fit, data = test_d)$predictions,
                                  time = as.numeric(difftime(t1, t0, units = "secs")),
                                  mem = as.numeric(object.size(et_fit)) / 1024)
  
  t0 <- Sys.time()
  xgb_params <- hp("XGBoost", list(max_depth = 4, eta = 0.05))
  xgb_fit <- do.call(fit_xgboost, c(list(train_d = train_d, response = response, predictors = predictors), xgb_params))
  t1 <- Sys.time()
  results[["XGBoost"]] <- list(pred = predict_xgboost(xgb_fit, test_d, predictors),
                               time = as.numeric(difftime(t1, t0, units = "secs")),
                               mem = as.numeric(object.size(xgb_fit)) / 1024)
  
  t0 <- Sys.time()
  gbm_params <- hp("GBM", list(interaction.depth = 3, shrinkage = 0.05))
  gbm_fit <- do.call(fit_gbm, c(list(train_d = train_d, response = response, predictors = predictors), gbm_params))
  t1 <- Sys.time()
  results[["GBM"]] <- list(pred = predict(gbm_fit, newdata = test_d, n.trees = gbm_fit$n.trees),
                           time = as.numeric(difftime(t1, t0, units = "secs")),
                           mem = as.numeric(object.size(gbm_fit)) / 1024)
  
  t0 <- Sys.time()
  rf_params <- hp("RandomForest", list(mtry = NULL, nodesize = 5))
  rf_fit <- do.call(fit_rf, c(list(train_d = train_d, response = response, predictors = predictors), rf_params))
  t1 <- Sys.time()
  results[["RandomForest"]] <- list(pred = predict(rf_fit, newdata = test_d),
                                    time = as.numeric(difftime(t1, t0, units = "secs")),
                                    mem = as.numeric(object.size(rf_fit)) / 1024)
  
  t0 <- Sys.time()
  svm_params <- hp("SVM", list(cost = 1, gamma = NULL))
  svm_fit <- do.call(fit_svm, c(list(train_d = train_d, response = response, predictors = predictors), svm_params))
  t1 <- Sys.time()
  results[["SVM"]] <- list(pred = predict(svm_fit, newdata = test_d),
                           time = as.numeric(difftime(t1, t0, units = "secs")),
                           mem = as.numeric(object.size(svm_fit)) / 1024)
  
  t0 <- Sys.time()
  voting_pred <- predict_voting_ensemble(gbm_fit, rf_fit, svm_fit, test_d)
  t1 <- Sys.time()
  results[["VotingEnsemble"]] <- list(pred = voting_pred,
                                      time = as.numeric(difftime(t1, t0, units = "secs")),
                                      mem = (as.numeric(object.size(gbm_fit)) + as.numeric(object.size(rf_fit)) +
                                               as.numeric(object.size(svm_fit))) / 1024)
  
  t0 <- Sys.time()
  stack_fit <- fit_stacking_ensemble(train_d, response, predictors,
                                     svm_params = svm_params, rf_params = rf_params, gbm_params = gbm_params)
  t1 <- Sys.time()
  results[["StackingEnsemble"]] <- list(pred = predict_stacking_ensemble(stack_fit, test_d),
                                        time = as.numeric(difftime(t1, t0, units = "secs")),
                                        mem = as.numeric(object.size(stack_fit)) / 1024)
  
  # Bootstrap standard errors for MAE/RMSE/R2: resamples (actual,
  # predicted) PAIRS from the test set with replacement (case resampling),
  # recomputes each metric per replicate, and uses the SD across replicates
  # as the SE. This does NOT refit any model -- it reuses the already-
  # computed test-set predictions -- so it properly reflects sampling
  # variability in the held-out test set itself, rather than presenting
  # point estimates as if computed on infinite data, at negligible
  # additional computational cost.
  n_boot_metrics <- 1000
  n_test <- length(actual)
  compute_metrics <- function(a, p) {
    c(MAE = mean(abs(a - p)), RMSE = sqrt(mean((a - p)^2)),
      R2 = 1 - sum((a - p)^2) / sum((a - mean(a))^2))
  }
  summary_df <- bind_rows(lapply(names(results), function(nm) {
    p <- results[[nm]]$pred
    point_metrics <- compute_metrics(actual, p)
    boot_metrics <- replicate(n_boot_metrics, {
      idx <- sample(seq_len(n_test), n_test, replace = TRUE)
      tryCatch(compute_metrics(actual[idx], p[idx]), error = function(e) c(MAE = NA, RMSE = NA, R2 = NA))
    })
    data.frame(Model = nm,
               MAE = unname(point_metrics["MAE"]), MAE_SE = sd(boot_metrics["MAE", ], na.rm = TRUE),
               RMSE = unname(point_metrics["RMSE"]), RMSE_SE = sd(boot_metrics["RMSE", ], na.rm = TRUE),
               R2 = unname(point_metrics["R2"]), R2_SE = sd(boot_metrics["R2", ], na.rm = TRUE),
               TrainTime_sec = results[[nm]]$time, Memory_KB = results[[nm]]$mem)
  })) %>% arrange(MAE)
  
  list(summary = summary_df, predictions = lapply(results, function(r) r$pred), actual = actual)
}

#############################################################################
# 17. Statistical properties of selected features
#############################################################################
feature_category_map <- c(
  Open = "Market data", High = "Market data", Low = "Market data", Volume = "Market data",
  Year = "Temporal", Month = "Temporal", DayOfWeek = "Temporal", Quarter = "Temporal",
  Return = "Price-based", LogReturn = "Price-based", PriceChange = "Price-based", HL_Range = "Price-based",
  CO_Ratio = "Price-based", SMA_5 = "Price-based", SMA_10 = "Price-based", SMA_20 = "Price-based",
  SMA_50 = "Price-based", EMA_12 = "Price-based", EMA_26 = "Price-based",
  Price_to_SMA20 = "Price-based", TypicalPrice = "Price-based",
  RollSD_5 = "Volatility", RollSD_10 = "Volatility", RollSD_20 = "Volatility", ATR_14 = "Volatility",
  BB_Upper = "Volatility", BB_Lower = "Volatility", BB_Width = "Volatility", ParkinsonVol = "Volatility",
  RSI_14 = "Momentum", MACD = "Momentum", MACD_Signal = "Momentum", MACD_Hist = "Momentum",
  Stoch_K = "Momentum", Stoch_D = "Momentum", ROC_10 = "Momentum", Momentum_10 = "Momentum",
  WilliamsR = "Momentum", CCI = "Momentum",
  Volume_SMA_5 = "Volume-based", Volume_SMA_20 = "Volume-based", Volume_Ratio = "Volume-based",
  OBV = "Volume-based", VPT = "Volume-based", CMF = "Volume-based", Volume_ROC = "Volume-based",
  DollarVolume = "Volume-based",
  HL_Spread = "Market microstructure", Amihud = "Market microstructure", RollSpread = "Market microstructure",
  PriceImpact = "Market microstructure", KyleLambda = "Market microstructure",
  IntradayVol = "Market microstructure", CLV = "Market microstructure",
  Gold = "Macro/Cross-asset", Gold_Return = "Macro/Cross-asset",
  CrudeOil = "Macro/Cross-asset", CrudeOil_Return = "Macro/Cross-asset",
  SP500 = "Macro/Cross-asset", SP500_Return = "Macro/Cross-asset",
  EURUSD = "Macro/Cross-asset", EURUSD_Return = "Macro/Cross-asset",
  GBPUSD = "Macro/Cross-asset", GBPUSD_Return = "Macro/Cross-asset",
  Nikkei = "Macro/Cross-asset", Nikkei_Return = "Macro/Cross-asset"
)

get_feature_category <- function(feature_names) {
  cats <- feature_category_map[feature_names]
  cats[is.na(cats)] <- "Unclassified"
  unname(cats)
}

wald_test_bootstrap <- function(fit, boot_coefs) {
  point_est <- coef(fit)
  boot_se <- apply(boot_coefs, 2, sd, na.rm = TRUE)
  wald_stat <- point_est / boot_se
  p_value <- 2 * (1 - pnorm(abs(wald_stat)))
  data.frame(Feature = names(point_est), Estimate = as.numeric(point_est),
             Bootstrap_SE = as.numeric(boot_se), Wald_Statistic = as.numeric(wald_stat),
             p_value = p_value, Significant_0.05 = p_value < 0.05)
}

feature_statistical_properties <- function(selected_features, group1_train, response = "Close",
                                           vif_table, stability_pi_hat, wald_table,
                                           corr_threshold = 0.3) {
  corr_vals <- sapply(selected_features, function(f) {
    tryCatch(cor(group1_train[[f]], group1_train[[response]], use = "pairwise.complete.obs"),
             error = function(e) NA)
  })
  vif_lookup <- setNames(vif_table$VIF, vif_table$Feature)
  
  data.frame(
    Feature = selected_features,
    Category = get_feature_category(selected_features),
    Correlation = as.numeric(corr_vals),
    Predictive_Signal_rho_gt_0.3 = abs(as.numeric(corr_vals)) > corr_threshold,
    VIF = vif_lookup[selected_features],
    Multicollinearity_Flag_VIF_gt_10 = vif_lookup[selected_features] > 10,
    Bootstrap_Selection_Freq = as.numeric(stability_pi_hat[selected_features]),
    Wald_Statistic = wald_table$Wald_Statistic[match(selected_features, wald_table$Feature)],
    Wald_p_value = wald_table$p_value[match(selected_features, wald_table$Feature)],
    Significant_0.05 = wald_table$Significant_0.05[match(selected_features, wald_table$Feature)]
  )
}

#############################################################################
# 18. Paired t-tests
#############################################################################
paired_model_ttests <- function(predictions, actual, baseline = "Huber") {
  base_err <- abs(actual - predictions[[baseline]])
  others <- setdiff(names(predictions), baseline)
  bind_rows(lapply(others, function(m) {
    err_m <- abs(actual - predictions[[m]])
    tt <- t.test(base_err, err_m, paired = TRUE)
    data.frame(Comparison = paste(baseline, "vs", m),
               Mean_AbsError_Huber = mean(base_err), Mean_AbsError_Other = mean(err_m),
               Mean_Diff = mean(base_err - err_m), t_statistic = unname(tt$statistic),
               df = unname(tt$parameter), p_value = tt$p.value,
               CI_Lower = tt$conf.int[1], CI_Upper = tt$conf.int[2],
               Significant_0.05 = tt$p.value < 0.05)
  }))
}

#############################################################################
# 19. Visual comparison
#############################################################################
plot_actual_vs_predicted_coins <- function(fit, split_list, selected_features, coins,
                                           split_name = "test", response = "Close") {
  par(mfrow = c(2, 2))
  for (coin in coins) {
    d <- split_list[[coin]][[split_name]]
    d <- d[complete.cases(d[, c(selected_features, response)]), ]
    pred <- predict(fit, newdata = d)
    actual <- d[[response]]
    plot(actual, pred, xlab = "Actual Close (normalized)", ylab = "Predicted Close (normalized)",
         main = paste("Huber: Actual vs Predicted -", coin),
         pch = 19, col = adjustcolor("steelblue", 0.5))
    abline(0, 1, col = "red", lwd = 2, lty = 2)
  }
  par(mfrow = c(1, 1))
}

plot_actual_vs_predicted_group1 <- function(fit, split_list, selected_features, response = "Close") {
  plot_actual_vs_predicted_coins(fit, split_list, selected_features, coins = GROUP1,
                                 split_name = "test", response = response)
}

per_coin_bias <- function(fit, split_list, selected_features, coins, split_name = "test", response = "Close") {
  bind_rows(lapply(coins, function(coin) {
    d <- split_list[[coin]][[split_name]]
    d <- d[complete.cases(d[, c(selected_features, response)]), ]
    pred <- predict(fit, newdata = d)
    actual <- d[[response]]
    err <- actual - pred
    data.frame(Symbol = coin, N = nrow(d), Mean_Error = mean(err), Median_Error = median(err),
               MAE = mean(abs(err)), Pct_Underpredicted = mean(err > 0) * 100)
  }))
}

prediction_error_distribution <- function(fit, test_d, selected_features, response = "Close") {
  test_d <- test_d[complete.cases(test_d[, c(selected_features, response)]), ]
  pred <- predict(fit, newdata = test_d)
  err <- test_d[[response]] - pred
  
  probs <- list("68%" = c(0.16, 0.84), "95%" = c(0.025, 0.975), "99%" = c(0.005, 0.995))
  bounds <- bind_rows(lapply(names(probs), function(lbl) {
    q <- quantile(err, probs = probs[[lbl]])
    data.frame(Level = lbl, Lower = q[1], Upper = q[2])
  }))
  
  hist(err, breaks = 50, col = "lightgray", border = "white",
       main = "Prediction Error Distribution (Huber, Group 1 Test)", xlab = "Error (Actual - Predicted)")
  cols <- c("darkgreen", "orange", "red")
  for (i in seq_len(nrow(bounds))) abline(v = c(bounds$Lower[i], bounds$Upper[i]), col = cols[i], lty = 2, lwd = 2)
  legend("topright", legend = bounds$Level, col = cols, lty = 2, lwd = 2, title = "Coverage")
  
  bounds
}

#############################################################################
# 20. Bootstrap uncertainty quantification for test-set predictions
#############################################################################
bootstrap_prediction_intervals <- function(train_d, test_d, response, predictors, B = 1000, seed = 2024,
                                           symbol_col = "Symbol", block_length = 10, calib_d = NULL) {
  set.seed(seed)
  form <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  n_train <- nrow(train_d)
  n_test <- nrow(test_d)
  n_calib <- if (!is.null(calib_d)) nrow(calib_d) else 0
  boot_point_preds <- matrix(NA, nrow = B, ncol = n_test)
  boot_preds <- matrix(NA, nrow = B, ncol = n_test)
  boot_point_preds_calib <- if (n_calib > 0) matrix(NA, nrow = B, ncol = n_calib) else NULL
  
  has_symbol <- symbol_col %in% names(train_d) && symbol_col %in% names(test_d)
  test_symbols <- if (has_symbol) test_d[[symbol_col]] else rep("ALL", n_test)
  
  block_bootstrap_noise <- function(resid_pool_ordered, target_length, block_length) {
    n_pool <- length(resid_pool_ordered)
    if (n_pool == 0) return(rep(NA_real_, target_length))
    bl <- min(block_length, n_pool)
    noise <- numeric(0)
    while (length(noise) < target_length) {
      start <- sample.int(n_pool, 1)
      idx <- ((start - 1 + seq_len(bl) - 1) %% n_pool) + 1
      noise <- c(noise, resid_pool_ordered[idx])
    }
    noise[seq_len(target_length)]
  }
  
  for (b in 1:B) {
    idx <- sample(1:n_train, n_train, replace = TRUE)
    fit_b <- tryCatch(MASS::rlm(form, data = train_d[idx, ], psi = MASS::psi.huber, maxit = 50),
                      error = function(e) NULL)
    if (is.null(fit_b)) next
    point_b <- predict(fit_b, newdata = test_d)
    boot_point_preds[b, ] <- point_b
    
    if (!is.null(calib_d)) {
      boot_point_preds_calib[b, ] <- tryCatch(predict(fit_b, newdata = calib_d),
                                              error = function(e) rep(NA_real_, n_calib))
    }
    
    resid_full_b <- tryCatch(train_d[[response]] - predict(fit_b, newdata = train_d),
                             error = function(e) NULL)
    if (is.null(resid_full_b)) next
    
    noise_b <- rep(NA_real_, n_test)
    if (has_symbol) {
      for (sym in unique(test_symbols)) {
        test_idx_sym <- which(test_symbols == sym)
        train_idx_sym <- which(train_d[[symbol_col]] == sym)
        pool <- resid_full_b[train_idx_sym]
        if (length(pool) == 0) pool <- resid_full_b
        noise_b[test_idx_sym] <- block_bootstrap_noise(pool, length(test_idx_sym), block_length)
      }
    } else {
      noise_b <- block_bootstrap_noise(resid_full_b, n_test, block_length)
    }
    boot_preds[b, ] <- point_b + noise_b
  }
  
  list(boot_preds = boot_preds,
       point_pred = colMeans(boot_point_preds, na.rm = TRUE),
       point_pred_calib = if (!is.null(calib_d)) colMeans(boot_point_preds_calib, na.rm = TRUE) else NULL,
       ci_lower = apply(boot_preds, 2, quantile, probs = 0.025, na.rm = TRUE),
       ci_upper = apply(boot_preds, 2, quantile, probs = 0.975, na.rm = TRUE))
}

conformalize_prediction_intervals <- function(boot_result, calib_d, test_d, response = "Close",
                                              symbol_col = "Symbol", alpha = 0.05) {
  point_pred_calib <- boot_result$point_pred_calib
  actual_calib <- calib_d[[response]]
  half_width_raw_test <- (boot_result$ci_upper - boot_result$ci_lower) / 2
  
  compute_scale_factor <- function(actual_c, pred_c, ref_half_width) {
    scores <- abs(actual_c - pred_c) / ref_half_width
    n <- length(scores)
    if (n == 0 || !is.finite(ref_half_width) || ref_half_width == 0) return(1)
    q_level <- min(1, ceiling((n + 1) * (1 - alpha)) / n)
    as.numeric(quantile(scores, probs = q_level, na.rm = TRUE, type = 7))
  }
  
  has_symbol <- symbol_col %in% names(calib_d) && symbol_col %in% names(test_d)
  if (has_symbol) {
    coins <- unique(calib_d[[symbol_col]])
    ref_half_width_by_coin <- setNames(sapply(coins, function(sym) {
      idx_test <- which(test_d[[symbol_col]] == sym)
      if (length(idx_test) == 0) return(mean(half_width_raw_test, na.rm = TRUE))
      mean(half_width_raw_test[idx_test], na.rm = TRUE)
    }), coins)
    
    scale_by_coin <- setNames(sapply(coins, function(sym) {
      idx_calib <- which(calib_d[[symbol_col]] == sym)
      compute_scale_factor(actual_calib[idx_calib], point_pred_calib[idx_calib], ref_half_width_by_coin[sym])
    }), coins)
  } else {
    ref_half_width_by_coin <- c(ALL = mean(half_width_raw_test, na.rm = TRUE))
    scale_by_coin <- c(ALL = compute_scale_factor(actual_calib, point_pred_calib, ref_half_width_by_coin["ALL"]))
  }
  
  scale_by_coin <- pmax(scale_by_coin, 1)
  
  list(scale_by_coin = scale_by_coin, ref_half_width_by_coin = ref_half_width_by_coin)
}

apply_conformal_scaling <- function(boot_result, target_d, scale_info, symbol_col = "Symbol") {
  half_width_raw <- (boot_result$ci_upper - boot_result$ci_lower) / 2
  scale_vec <- rep(scale_info$scale_by_coin[1], length(half_width_raw))
  if (symbol_col %in% names(target_d) && length(scale_info$scale_by_coin) > 1) {
    for (sym in names(scale_info$scale_by_coin)) {
      idx <- which(target_d[[symbol_col]] == sym)
      if (length(idx) > 0) scale_vec[idx] <- scale_info$scale_by_coin[sym]
    }
  }
  new_half_width <- half_width_raw * scale_vec
  list(point_pred = boot_result$point_pred,
       ci_lower = boot_result$point_pred - new_half_width,
       ci_upper = boot_result$point_pred + new_half_width,
       scale_applied = scale_vec)
}

plot_bootstrap_prediction_ci <- function(test_d, boot_result, response = "Close",
                                         date_col = "Date", title = "") {
  ord <- order(test_d[[date_col]])
  x <- test_d[[date_col]][ord]
  actual <- test_d[[response]][ord]
  pred <- boot_result$point_pred[ord]
  lower <- boot_result$ci_lower[ord]
  upper <- boot_result$ci_upper[ord]
  
  plot(x, actual, type = "n", xlab = "Date", ylab = "Close (normalized)",
       main = paste("Bootstrap 95% Prediction CI -", title), ylim = range(c(actual, lower, upper), na.rm = TRUE))
  polygon(c(x, rev(x)), c(upper, rev(lower)), col = adjustcolor("gray", 0.5), border = NA)
  lines(x, pred, col = "red", lwd = 2)
  points(x, actual, col = "blue", pch = 19, cex = 0.6)
  legend("topleft", legend = c("Actual", "Point Prediction", "95% CI"),
         col = c("blue", "red", "gray"), pch = c(19, NA, 15), lty = c(NA, 1, NA),
         lwd = c(NA, 2, NA), bty = "n")
}

plot_bootstrap_prediction_ci_group1 <- function(test_d, boot_result, group1 = GROUP1,
                                                response = "Close", date_col = "Date",
                                                symbol_col = "Symbol") {
  par(mfrow = c(2, 2))
  for (coin in group1) {
    idx <- which(test_d[[symbol_col]] == coin)
    sub_test <- test_d[idx, ]
    sub_result <- list(point_pred = boot_result$point_pred[idx],
                       ci_lower = boot_result$ci_lower[idx],
                       ci_upper = boot_result$ci_upper[idx])
    plot_bootstrap_prediction_ci(sub_test, sub_result, response = response,
                                 date_col = date_col, title = coin)
  }
  par(mfrow = c(1, 1))
}

bootstrap_coefficient_table <- function(fit, boot_coefs) {
  point_est <- coef(fit)
  se <- apply(boot_coefs, 2, sd, na.rm = TRUE)
  ci <- apply(boot_coefs, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  data.frame(Feature = names(point_est), `Point Estimate` = as.numeric(point_est),
             `Bootstrap SE` = as.numeric(se), `95% CI Lower` = ci[1, ], `95% CI Upper` = ci[2, ],
             check.names = FALSE)
}

bootstrap_ci_summary <- function(test_d, boot_result, response = "Close", symbol_col = "Symbol") {
  ci_width <- boot_result$ci_upper - boot_result$ci_lower
  actual <- test_d[[response]]
  covered <- (actual >= boot_result$ci_lower) & (actual <= boot_result$ci_upper)
  
  overall <- data.frame(
    Metric = c("Mean CI width", "Median CI width", "Min CI width", "Max CI width", "Coverage rate"),
    Value = c(round(mean(ci_width, na.rm = TRUE), 5), round(median(ci_width, na.rm = TRUE), 5),
              round(min(ci_width, na.rm = TRUE), 5), round(max(ci_width, na.rm = TRUE), 5),
              paste0(round(100 * mean(covered, na.rm = TRUE), 1), "%"))
  )
  
  by_coin <- NULL
  if (symbol_col %in% names(test_d)) {
    by_coin <- bind_rows(lapply(unique(test_d[[symbol_col]]), function(coin) {
      idx <- which(test_d[[symbol_col]] == coin)
      data.frame(Symbol = coin,
                 Mean_CI_Width = mean(ci_width[idx], na.rm = TRUE),
                 Median_CI_Width = median(ci_width[idx], na.rm = TRUE),
                 Min_CI_Width = min(ci_width[idx], na.rm = TRUE),
                 Max_CI_Width = max(ci_width[idx], na.rm = TRUE),
                 Coverage_Rate = mean(covered[idx], na.rm = TRUE))
    }))
  }
  
  list(overall = overall, by_coin = by_coin)
}

#############################################################################
# 21. Feature-selection stability comparison
#############################################################################
lasso_cv_selection <- function(X, y, feature_names) {
  Xm <- as.matrix(X[, feature_names, drop = FALSE])
  cvfit <- glmnet::cv.glmnet(Xm, y, alpha = 1)
  coefs <- as.matrix(coef(cvfit, s = "lambda.min"))
  setdiff(rownames(coefs)[coefs[, 1] != 0], "(Intercept)")
}

# SCAD-penalized regression (Fan & Li, 2001). Unlike LASSO's constant
# penalty rate, SCAD's penalty derivative decays to zero for large
# coefficients, so it does not shrink genuinely large, well-estimated
# coefficients the way LASSO does -- addressing LASSO's known bias toward
# over-penalizing strong signals. Fit via ncvreg (the standard R
# implementation), with lambda chosen by cross-validation.
scad_selection <- function(X, y, feature_names) {
  Xm <- as.matrix(X[, feature_names, drop = FALSE])
  cvfit <- tryCatch(ncvreg::cv.ncvreg(Xm, y, penalty = "SCAD"), error = function(e) NULL)
  if (is.null(cvfit)) return(character(0))
  coefs <- tryCatch(as.numeric(coef(cvfit, lambda = cvfit$lambda.min)), error = function(e) NULL)
  if (is.null(coefs)) return(character(0))
  nm <- names(coef(cvfit, lambda = cvfit$lambda.min))
  setdiff(nm[coefs != 0], "(Intercept)")
}

# Adaptive LASSO (Zou, 2006). Standard two-step procedure: (1) fit an
# initial consistent estimator (ridge, alpha=0, since OLS is unreliable
# once collinear candidate features are involved) to obtain a coefficient
# magnitude for each feature; (2) re-fit LASSO with per-feature penalty
# weights = 1/|beta_initial|^gamma, so features with larger initial
# coefficients are penalized LESS and near-zero-coefficient features are
# penalized MORE (approaching exclusion). This adaptive weighting is what
# gives the method its oracle property (asymptotically consistent variable
# selection) that plain LASSO does not have.
adaptive_lasso_selection <- function(X, y, feature_names, gamma = 1) {
  Xm <- as.matrix(X[, feature_names, drop = FALSE])
  init_fit <- tryCatch(glmnet::cv.glmnet(Xm, y, alpha = 0), error = function(e) NULL)
  if (is.null(init_fit)) return(character(0))
  init_coefs <- as.numeric(coef(init_fit, s = "lambda.min"))[-1]
  adapt_weights <- 1 / (abs(init_coefs)^gamma + 1e-8)
  cvfit <- tryCatch(glmnet::cv.glmnet(Xm, y, alpha = 1, penalty.factor = adapt_weights), error = function(e) NULL)
  if (is.null(cvfit)) return(character(0))
  coefs <- as.matrix(coef(cvfit, s = "lambda.min"))
  setdiff(rownames(coefs)[coefs[, 1] != 0], "(Intercept)")
}

rfe_selection <- function(X, y, feature_names, size = NULL) {
  if (is.null(size)) size <- max(3, round(length(feature_names) / 3))
  ctrl <- caret::rfeControl(functions = caret::lmFuncs, method = "cv", number = 5)
  rfe_fit <- tryCatch(caret::rfe(x = X[, feature_names, drop = FALSE], y = y, sizes = size, rfeControl = ctrl),
                      error = function(e) NULL)
  if (is.null(rfe_fit)) return(character(0))
  
  var_imp <- tryCatch(caret::varImp(rfe_fit), error = function(e) NULL)
  if (is.null(var_imp) || nrow(var_imp) == 0) return(caret::predictors(rfe_fit))
  ranked <- rownames(var_imp)[order(-var_imp$Overall)]
  ranked[1:min(size, length(ranked))]
}

rf_importance_selection <- function(X, y, feature_names, k = NULL) {
  if (is.null(k)) k <- max(3, round(length(feature_names) / 3))
  dat <- data.frame(y = y, X[, feature_names, drop = FALSE])
  rf_fit <- randomForest::randomForest(y ~ ., data = dat, ntree = 200, importance = TRUE)
  imp <- randomForest::importance(rf_fit)[, "%IncMSE"]
  names(sort(imp, decreasing = TRUE))[1:min(k, length(imp))]
}

forward_stepwise_selection <- function(X, y, feature_names, k = NULL) {
  if (is.null(k)) k <- max(3, round(length(feature_names) / 3))
  dat <- data.frame(y = y, X[, feature_names, drop = FALSE])
  dat <- dat[complete.cases(dat), ]
  if (nrow(dat) < length(feature_names) + 5) return(character(0))
  
  null_model <- lm(y ~ 1, data = dat)
  full_formula <- as.formula(paste("y ~", paste(feature_names, collapse = " + ")))
  
  step_model <- tryCatch(
    stats::step(null_model, scope = list(lower = ~1, upper = full_formula),
                direction = "forward", steps = k, trace = 0),
    error = function(e) NULL)
  if (is.null(step_model)) return(character(0))
  
  setdiff(names(coef(step_model)), "(Intercept)")
}

bic_stepwise_selection <- function(X, y, feature_names, k = NULL) {
  if (is.null(k)) k <- max(3, round(length(feature_names) / 3))
  dat <- data.frame(y = y, X[, feature_names, drop = FALSE])
  dat <- dat[complete.cases(dat), ]
  if (nrow(dat) < length(feature_names) + 5) return(character(0))
  
  n <- nrow(dat)
  null_model <- lm(y ~ 1, data = dat)
  full_formula <- as.formula(paste("y ~", paste(feature_names, collapse = " + ")))
  
  step_model <- tryCatch(
    stats::step(null_model, scope = list(lower = ~1, upper = full_formula),
                direction = "forward", steps = k, k = log(n), trace = 0),
    error = function(e) NULL)
  if (is.null(step_model)) return(character(0))
  
  setdiff(names(coef(step_model)), "(Intercept)")
}

generic_stability_selection <- function(X, y, feature_names, selector_fn, B = 15, pi_thr = 0.6, ...) {
  n <- nrow(X)
  selections <- list()
  selmat <- matrix(0, nrow = B, ncol = length(feature_names), dimnames = list(NULL, feature_names))
  for (b in 1:B) {
    # Same M&B (2010) subsampling as stability_selection() above -- without
    # replacement, size floor(n/2) -- so every method in Section 21's
    # comparison is evaluated under the IDENTICAL resampling scheme. Using
    # a fuller/gentler resampling for the benchmarks and a stricter one for
    # Multi-Strategy Integration would make the stability comparison
    # apples-to-oranges regardless of which scheme is "more correct."
    idx <- sample(1:n, size = floor(n / 2), replace = FALSE)
    sel <- tryCatch(selector_fn(X[idx, , drop = FALSE], y[idx], feature_names, ...), error = function(e) character(0))
    selections[[b]] <- sel
    if (length(sel)) selmat[b, sel] <- 1
  }
  pi_hat <- colMeans(selmat)
  jaccard <- function(s1, s2) {
    if (length(s1) == 0 && length(s2) == 0) return(1)
    length(intersect(s1, s2)) / length(union(s1, s2))
  }
  jvals <- c()
  if (B > 1) for (i in 1:(B - 1)) for (j in (i + 1):B) jvals <- c(jvals, jaccard(selections[[i]], selections[[j]]))

  # Same bootstrap-of-replicates SE/CI as stability_selection() -- see the
  # detailed rationale there. Reused identically here for consistency
  # across every method in the Section 21 comparison.
  n_boot_se <- 500
  boot_means <- replicate(n_boot_se, {
    boot_idx <- sample(seq_len(B), B, replace = TRUE)
    boot_sel <- selections[boot_idx]
    bvals <- c()
    if (B > 1) for (i in 1:(B - 1)) for (j in (i + 1):B) bvals <- c(bvals, jaccard(boot_sel[[i]], boot_sel[[j]]))
    if (length(bvals) == 0) NA_real_ else mean(bvals)
  })
  jaccard_se <- sd(boot_means, na.rm = TRUE)
  jaccard_ci <- quantile(boot_means, c(0.025, 0.975), na.rm = TRUE, type = 7)

  list(pi_hat = pi_hat, stable_features = names(pi_hat)[pi_hat >= pi_thr],
       avg_jaccard_stability = mean(jvals), jaccard_se = jaccard_se,
       jaccard_ci_lower = unname(jaccard_ci[1]), jaccard_ci_upper = unname(jaccard_ci[2]),
       selections = selections)
}

compare_feature_selection_stability <- function(X, y, feature_names, B = 15, pi_thr = 0.6, k = NULL,
                                                quantum_stability_result = NULL) {
  if (is.null(k)) k <- max(3, round(length(feature_names) / 3))
  
  quantum_stab <- if (!is.null(quantum_stability_result)) quantum_stability_result else
    stability_selection(X, y, feature_names, B = B, pi_thr = pi_thr)
  
  lasso_stab <- generic_stability_selection(X, y, feature_names, lasso_cv_selection, B = B, pi_thr = pi_thr)
  rfe_stab <- generic_stability_selection(X, y, feature_names, rfe_selection, B = B, pi_thr = pi_thr, size = k)
  rfimp_stab <- generic_stability_selection(X, y, feature_names, rf_importance_selection, B = B, pi_thr = pi_thr, k = k)
  stepwise_stab <- generic_stability_selection(X, y, feature_names, forward_stepwise_selection, B = B, pi_thr = pi_thr, k = k)
  bic_stab <- generic_stability_selection(X, y, feature_names, bic_stepwise_selection, B = B, pi_thr = pi_thr, k = k)
  scad_stab <- generic_stability_selection(X, y, feature_names, scad_selection, B = B, pi_thr = pi_thr)
  adalasso_stab <- generic_stability_selection(X, y, feature_names, adaptive_lasso_selection, B = B, pi_thr = pi_thr)
  
  data.frame(
    Method = c("Multi-Strategy Integration", "LASSO (CV)", "RFE", "RandomForest Importance",
               "Forward Stepwise (AIC)", "Forward Stepwise (BIC)", "SCAD", "Adaptive LASSO"),
    Avg_Jaccard_Stability = c(quantum_stab$avg_jaccard_stability, lasso_stab$avg_jaccard_stability,
                              rfe_stab$avg_jaccard_stability, rfimp_stab$avg_jaccard_stability,
                              stepwise_stab$avg_jaccard_stability, bic_stab$avg_jaccard_stability,
                              scad_stab$avg_jaccard_stability, adalasso_stab$avg_jaccard_stability),
    Jaccard_SE = c(quantum_stab$jaccard_se, lasso_stab$jaccard_se, rfe_stab$jaccard_se,
                   rfimp_stab$jaccard_se, stepwise_stab$jaccard_se, bic_stab$jaccard_se,
                   scad_stab$jaccard_se, adalasso_stab$jaccard_se),
    Jaccard_CI_Lower = c(quantum_stab$jaccard_ci_lower, lasso_stab$jaccard_ci_lower, rfe_stab$jaccard_ci_lower,
                         rfimp_stab$jaccard_ci_lower, stepwise_stab$jaccard_ci_lower, bic_stab$jaccard_ci_lower,
                         scad_stab$jaccard_ci_lower, adalasso_stab$jaccard_ci_lower),
    Jaccard_CI_Upper = c(quantum_stab$jaccard_ci_upper, lasso_stab$jaccard_ci_upper, rfe_stab$jaccard_ci_upper,
                         rfimp_stab$jaccard_ci_upper, stepwise_stab$jaccard_ci_upper, bic_stab$jaccard_ci_upper,
                         scad_stab$jaccard_ci_upper, adalasso_stab$jaccard_ci_upper),
    N_Stable_Features = c(length(quantum_stab$stable_features), length(lasso_stab$stable_features),
                          length(rfe_stab$stable_features), length(rfimp_stab$stable_features),
                          length(stepwise_stab$stable_features), length(bic_stab$stable_features),
                          length(scad_stab$stable_features), length(adalasso_stab$stable_features))
  ) %>% arrange(desc(Avg_Jaccard_Stability))
}

#############################################################################
# 22. Group 2: universal vs. asset-specific features
#############################################################################
asset_specific_quantum_model <- function(train_d, valid_d, predictor_cols, response = "Close",
                                         pop_size = 20, max_iter = 30, seed = 42) {
  y <- train_d[[response]]
  note <- "ok"
  
  qfs_res <- tryCatch(
    multi_strategy_feature_selection(X = train_d[, predictor_cols], y = y, feature_names = predictor_cols,
                                     pop_size = pop_size, max_iter = max_iter, seed = seed),
    error = function(e) {
      note <<- paste("Multi-strategy feature selection failed:", conditionMessage(e))
      NULL
    })
  sel_feats <- if (!is.null(qfs_res)) qfs_res$selected_features else character(0)
  
  if (length(sel_feats) > 0) {
    keep <- sapply(sel_feats, function(f) length(unique(na.omit(train_d[[f]]))) > 1)
    if (any(!keep)) note <- paste0(note, "; dropped zero-variance: ", paste(sel_feats[!keep], collapse = ","))
    sel_feats <- sel_feats[keep]
  }
  if (length(sel_feats) == 0) {
    sel_feats <- predictor_cols[1]
    note <- paste0(note, "; quantum selection empty, fell back to single feature: ", sel_feats)
  }
  
  joint_k <- if (!is.null(qfs_res)) qfs_res$best_k else 1.345
  
  fit <- tryCatch(fit_huber(train_d, response, sel_feats, k = joint_k, maxit = 100),
                  error = function(e) {
                    note <<- paste0(note, "; Huber fit failed: ", conditionMessage(e))
                    NULL
                  })
  if (is.null(fit)) return(list(selected_features = sel_feats, R2 = NA_real_, MAE = NA_real_, note = note))
  
  valid_d2 <- valid_d[complete.cases(valid_d[, c(response, sel_feats)]), ]
  if (nrow(valid_d2) < 2) {
    note <- paste0(note, "; too few complete validation rows: ", nrow(valid_d2))
    return(list(selected_features = sel_feats, R2 = NA_real_, MAE = NA_real_, note = note))
  }
  
  pred <- tryCatch(predict(fit, newdata = valid_d2),
                   error = function(e) {
                     note <<- paste0(note, "; prediction failed: ", conditionMessage(e))
                     rep(NA_real_, nrow(valid_d2))
                   })
  actual <- valid_d2[[response]]
  if (all(is.na(pred))) {
    note <- paste0(note, "; all predictions NA")
    return(list(selected_features = sel_feats, R2 = NA_real_, MAE = NA_real_, note = note))
  }
  
  list(selected_features = sel_feats,
       R2 = 1 - sum((actual - pred)^2, na.rm = TRUE) / sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE),
       MAE = mean(abs(actual - pred), na.rm = TRUE), note = note)
}

compare_universal_vs_asset_specific <- function(final_huber_fit, split_list, predictor_cols,
                                                selected_features, group2 = GROUP2, response = "Close",
                                                qfs_pop_size = 20, qfs_max_iter = 30) {
  bind_rows(lapply(group2, function(coin) {
    train_d <- split_list[[coin]]$train
    valid_d <- split_list[[coin]]$valid
    
    valid_u <- valid_d[complete.cases(valid_d[, c(response, selected_features)]), ]
    pred_u <- predict(final_huber_fit, newdata = valid_u)
    actual_u <- valid_u[[response]]
    r2_u <- 1 - sum((actual_u - pred_u)^2) / sum((actual_u - mean(actual_u))^2)
    mae_u <- mean(abs(actual_u - pred_u))
    
    as_res <- asset_specific_quantum_model(train_d, valid_d, predictor_cols, response,
                                           pop_size = qfs_pop_size, max_iter = qfs_max_iter,
                                           seed = 42 + which(group2 == coin))
    
    data.frame(Symbol = coin, R2_Universal = r2_u, MAE_Universal = mae_u,
               R2_AssetSpecific = as_res$R2, MAE_AssetSpecific = as_res$MAE,
               N_AssetSpecific_Features = length(as_res$selected_features),
               AssetSpecific_Features = paste(as_res$selected_features, collapse = "; "),
               N_Train_Rows = nrow(train_d), N_Valid_Rows = nrow(valid_u),
               AssetSpecific_Note = as_res$note)
  }))
}

hypothesis_test_universal_vs_specific <- function(comparison_df) {
  complete_df <- comparison_df[complete.cases(comparison_df[, c("R2_Universal", "R2_AssetSpecific")]), ]
  
  if (nrow(complete_df) < 2) {
    warning("Only ", nrow(complete_df), " of ", nrow(comparison_df),
            " Group 2 coins have a complete (R2_Universal, R2_AssetSpecific) pair.")
    return(data.frame(
      Mean_R2_Universal = mean(comparison_df$R2_Universal, na.rm = TRUE),
      Mean_R2_AssetSpecific = mean(comparison_df$R2_AssetSpecific, na.rm = TRUE),
      Mean_MAE_Universal = mean(comparison_df$MAE_Universal, na.rm = TRUE),
      Mean_MAE_AssetSpecific = mean(comparison_df$MAE_AssetSpecific, na.rm = TRUE),
      t_statistic = NA_real_, df = NA_real_, p_value = NA_real_,
      Conclusion = paste0("Insufficient complete R2 pairs (", nrow(complete_df), " of ",
                          nrow(comparison_df), " coins) to run the paired t-test.")
    ))
  }
  
  tt <- t.test(complete_df$R2_Universal, complete_df$R2_AssetSpecific,
               paired = TRUE, alternative = "less")
  data.frame(
    Mean_R2_Universal = mean(comparison_df$R2_Universal, na.rm = TRUE),
    Mean_R2_AssetSpecific = mean(comparison_df$R2_AssetSpecific, na.rm = TRUE),
    Mean_MAE_Universal = mean(comparison_df$MAE_Universal, na.rm = TRUE),
    Mean_MAE_AssetSpecific = mean(comparison_df$MAE_AssetSpecific, na.rm = TRUE),
    t_statistic = unname(tt$statistic), df = unname(tt$parameter), p_value = tt$p.value,
    Conclusion = ifelse(tt$p.value < 0.05,
                        "Reject H0: universal features underperform asset-specific quantum-selected features (p < 0.05)",
                        "Fail to reject H0: no significant evidence universal features underperform asset-specific quantum-selected features")
  )
}

#############################################################################
# MAIN PIPELINE EXECUTION
#############################################################################

cat("\n=== Starting Cryptocurrency Analysis Pipeline ===\n")

## ---- 1. Descriptive statistics ----
cat("\n1. Computing descriptive statistics...\n")
desc_stats_table <- run_descriptive_stats_all_coins(feature_data_list, "Close")
write.csv(desc_stats_table, file.path(out_dir, "01_descriptive_stats.csv"), row.names = FALSE)

## ---- 3. Butterworth filtering + ACF comparison ----
cat("\n3. Applying Butterworth filter...\n")
feature_data_list <- lapply(feature_data_list, apply_butterworth_filter,
                            value_col = "Close", order = 3, fc = 0.05)
acf_comparison_results <- lapply(feature_data_list, compare_acf_filtered_unfiltered,
                                 value_col = "Close", plot = FALSE)

for (nm in names(feature_data_list)) {
  png(file.path(out_dir, paste0("03_acf_comparison_", nm, ".png")), width = 900, height = 450)
  compare_acf_filtered_unfiltered(feature_data_list[[nm]], value_col = "Close", plot = TRUE)
  dev.off()
}

acf_summary_table <- bind_rows(lapply(names(acf_comparison_results), function(nm) {
  s <- summarize_acf_comparison(acf_comparison_results[[nm]])
  s$Symbol <- nm
  s %>% relocate(Symbol)
}))
write.csv(acf_summary_table, file.path(out_dir, "03_acf_comparison_summary.csv"), row.names = FALSE)

## ---- Column bookkeeping ----
all_feature_cols <- setdiff(names(feature_data_list[[1]]), c("Date", "Symbol"))
outlier_check_cols <- setdiff(all_feature_cols, c("Year", "Month", "DayOfWeek", "Quarter"))

# Exclude target-leakage features
leakage_features <- c("PriceChange", "CO_Ratio", "TypicalPrice", "CLV", "Return", "LogReturn")
predictor_cols <- setdiff(all_feature_cols, c("Close", "Close_Filtered", leakage_features))
message("Excluded target-leakage features from predictor pool: ", paste(leakage_features, collapse = ", "))

## ---- 4. Outlier detection ----
cat("\n4. Detecting outliers...\n")
outlier_summary <- lapply(feature_data_list, outlier_report, cols = outlier_check_cols)
for (nm in names(outlier_summary)) {
  write.csv(outlier_summary[[nm]], file.path(out_dir, paste0("04_outliers_", nm, ".csv")), row.names = FALSE)
}

## Group 1 subset
group1_feature_data_list <- feature_data_list[GROUP1]

## ---- 5. Temporal properties ----
cat("\n5. Computing temporal properties...\n")
temporal_props_table <- run_temporal_properties_all(group1_feature_data_list, "Close")
write.csv(temporal_props_table, file.path(out_dir, "05_temporal_properties.csv"), row.names = FALSE)

## ---- 6. Correlation analysis ----
cat("\n6. Computing correlations...\n")
correlation_results <- lapply(group1_feature_data_list, correlation_analysis,
                              predictors = predictor_cols, response = "Close")
for (nm in names(correlation_results)) {
  write.csv(correlation_results[[nm]], file.path(out_dir, paste0("06_correlation_", nm, ".csv")), row.names = FALSE)
}

## ---- 7. VIF + PCA ----
cat("\n7. Computing VIF and PCA...\n")
vif_results <- lapply(group1_feature_data_list, function(df) vif_analysis(df, predictor_cols, "Close"))
for (nm in names(vif_results)) {
  write.csv(vif_results[[nm]], file.path(out_dir, paste0("07_vif_", nm, ".csv")), row.names = FALSE)
}

pca_results <- lapply(group1_feature_data_list, pca_analysis, predictors = predictor_cols)
for (nm in names(pca_results)) {
  pr <- pca_results[[nm]]
  pca_table <- data.frame(Component = seq_along(pr$var_explained),
                           VarExplained = pr$var_explained, CumVar = pr$cum_var,
                           N_Components_Needed_95pct = pr$n_components_needed)
  write.csv(pca_table, file.path(out_dir, paste0("07_pca_", nm, ".csv")), row.names = FALSE)
}
# (pr$pca, the full prcomp object, is intentionally not written to CSV --
# not naturally tabular; the loadings/scores can be reconstructed from
# pipeline_results.rds$pca_results[[coin]]$pca if needed.)

## ---- 8. Split + normalize ----
cat("\n8. Splitting data and normalizing...\n")
split_list_raw <- lapply(feature_data_list, split_dataset)

split_list <- lapply(split_list_raw, function(spl) {
  norm <- normalize_dataset(spl$train, c(predictor_cols, "Close"))
  list(train = norm$data,
       valid = apply_norm_params(spl$valid, norm$params),
       test = apply_norm_params(spl$test, norm$params),
       close_norm_params = norm$params[["Close"]])  # c(min=..., max=...), TRAIN-period range only
})

# Per-coin Close min/max (from TRAINING data only, matching how normalization
# was actually fit -- needed to convert normalized-scale RMSE back to raw
# dollar terms: RMSE_raw = RMSE_normalized * (max - min). Only the RANGE
# matters for RMSE (a dispersion statistic); the min itself does not enter.
close_norm_params <- lapply(split_list, function(spl) spl$close_norm_params)
write.csv(data.frame(Symbol = names(close_norm_params),
                      Train_Close_Min = sapply(close_norm_params, `[`, "min"),
                      Train_Close_Max = sapply(close_norm_params, `[`, "max"),
                      Train_Close_Range = sapply(close_norm_params, function(p) p["max"] - p["min"])),
          file.path(out_dir, "08_close_normalization_params.csv"), row.names = FALSE)

group1_train <- bind_rows(lapply(GROUP1, function(s) split_list[[s]]$train))
group1_valid <- bind_rows(lapply(GROUP1, function(s) split_list[[s]]$valid))
group1_test <- bind_rows(lapply(GROUP1, function(s) split_list[[s]]$test))
group1_train <- group1_train[complete.cases(group1_train[, c(predictor_cols, "Close")]), ]
group1_valid <- group1_valid[complete.cases(group1_valid[, c(predictor_cols, "Close")]), ]
group1_test <- group1_test[complete.cases(group1_test[, c(predictor_cols, "Close")]), ]
group2_valid_list <- setNames(lapply(GROUP2, function(s) split_list[[s]]$valid), GROUP2)

## ---- 9. Multi-Strategy Integration feature selection ----
cat("\n9. Running Multi-Strategy Integration feature selection...\n")
coin_row_counts <- table(group1_train$Symbol)
message("Group 1 training row counts by coin: ", paste(names(coin_row_counts), coin_row_counts, sep = "=", collapse = ", "))
huber_case_weights <- as.numeric(1 / coin_row_counts[group1_train$Symbol])
huber_case_weights <- huber_case_weights / mean(huber_case_weights)

N_RUNS_MULTISTRATEGY <- 20  # number of independent runs for convergence diagnostics; reduce if runtime is prohibitive
message("Running Multi-Strategy Integration ", N_RUNS_MULTISTRATEGY, " times independently for convergence diagnostics...")
qfs_multi <- run_multi_strategy_multiple(X = group1_train[, predictor_cols], y = group1_train$Close,
                                          feature_names = predictor_cols, n_runs = N_RUNS_MULTISTRATEGY,
                                          pop_size = 20, max_iter = 30,
                                          weights = huber_case_weights, groups = group1_train$Symbol)
qfs_result <- qfs_multi$best_run   # the single best-fitness run across N_RUNS_MULTISTRATEGY, kept as the "official" model
selected_features <- qfs_result$selected_features
joint_best_k <- qfs_result$best_k
message("Selected meta-features (from best of ", N_RUNS_MULTISTRATEGY, " runs): ", paste(selected_features, collapse = ", "))
message("Jointly-optimized Huber k: ", round(joint_best_k, 4))

write.csv(qfs_multi$summary_table, file.path(out_dir, "09_multirun_summary.csv"), row.names = FALSE)
write.csv(qfs_multi$selection_freq_table, file.path(out_dir, "09_multirun_feature_selection_frequency.csv"),
          row.names = FALSE)

# Convergence plot: overlay ALL N_RUNS_MULTISTRATEGY fitness trajectories.
# A properly converging search should show every run's curve flattening to
# a similar final fitness level well before the final generation --
# trajectories that flatten at visibly DIFFERENT levels, or that are still
# decreasing at the final generation, would indicate the search has not
# genuinely converged (either an insufficient max_iter, or a fitness
# landscape with multiple, substantially different local optima).
png(file.path(out_dir, "09_ga_convergence.png"), width = 900, height = 600)
all_histories <- sapply(qfs_multi$runs, function(r) r$history)  # matrix: generation x run
matplot(all_histories, type = "l", lty = 1, lwd = 1, col = adjustcolor("steelblue", alpha.f = 0.4),
        xlab = "Generation", ylab = "Best Fitness (lower = better)",
        main = paste0("Multi-Strategy Integration: Convergence Across ", N_RUNS_MULTISTRATEGY, " Independent Runs"))
lines(qfs_multi$best_run$history, lwd = 2.5, col = "firebrick")
legend("topright", legend = c(paste0("Individual runs (n=", N_RUNS_MULTISTRATEGY, ")"), "Best run (reported model)"),
       col = c("steelblue", "firebrick"), lty = 1, lwd = c(1, 2.5), bty = "n")
dev.off()
write.csv(data.frame(Generation = seq_along(qfs_result$history), Best_Fitness_Reported_Run = qfs_result$history),
          file.path(out_dir, "09_ga_convergence.csv"), row.names = FALSE)

# Sensitivity analysis: how much do the penalty weights and population
# parameters actually matter? Run once (single run per configuration,
# not the N_RUNS_MULTISTRATEGY multi-run diagnostic, to keep this
# tractable -- 15 configurations total across the 5 varied parameters).
message("Running sensitivity analysis over penalty weights and population parameters...")
sensitivity_results <- sensitivity_analysis_multi_strategy(
  X = group1_train[, predictor_cols], y = group1_train$Close, feature_names = predictor_cols,
  weights = huber_case_weights, groups = group1_train$Symbol)
write.csv(sensitivity_results, file.path(out_dir, "09_sensitivity_analysis.csv"), row.names = FALSE)

final_huber_fit <- fit_huber(group1_train, "Close", selected_features, k = joint_best_k, maxit = 100,
                             weights = huber_case_weights)

## ---- 9b. Validation check ----
cat("\n9b. Validating on Group 1 validation set...\n")
validation_check <- out_of_sample_validation(final_huber_fit, group1_valid, "Close")
write.csv(validation_check, file.path(out_dir, "09b_validation_check_group1.csv"), row.names = FALSE)
message("Group 1 VALIDATION check -- MAE: ", round(validation_check$MAE, 4),
        ", R2: ", round(validation_check$R2, 4))

## ---- 10 & 11. Residual diagnostics ----
cat("\n10-11. Computing residual diagnostics...\n")
res_diag <- residual_diagnostics(final_huber_fit, group1_train)
res_diag_table <- residual_diagnostics_table(res_diag)
write.csv(res_diag_table, file.path(out_dir, "10_residual_diagnostic_tests.csv"), row.names = FALSE)

png(file.path(out_dir, "11_residual_diagnostics.png"), width = 900, height = 900)
plot_residual_diagnostics(final_huber_fit, group1_train, response = "Close")
dev.off()

## ---- 12. Bootstrap ----
cat("\n12. Running residual bootstrap (B=1000)...\n")
boot_res <- residual_bootstrap(final_huber_fit, group1_train, response = "Close", B = 1000)
write.csv(as.data.frame(boot_res$ci_95), file.path(out_dir, "12_bootstrap_ci.csv"))

## ---- 13. Stability selection ----
cat("\n13. Running stability selection...\n")
stability_res <- stability_selection(X = group1_train[, predictor_cols], y = group1_train$Close,
                                      feature_names = predictor_cols, B = 30, pi_thr = 0.6,
                                      qfs_args = list(pop_size = 15, max_iter = 15))

write.csv(data.frame(Feature = names(stability_res$pi_hat), SelectionProb = stability_res$pi_hat),
          file.path(out_dir, "13_stability_selection.csv"), row.names = FALSE)

## ---- 14. Cross-validation ----
cat("\n14. Running cross-validation...\n")
cv10 <- kfold_cv(group1_train, "Close", selected_features, K = 10)
expanding_cv <- expanding_window_cv(group1_train, "Close", selected_features)
write.csv(cv10, file.path(out_dir, "14_kfold_cv.csv"), row.names = FALSE)
write.csv(expanding_cv, file.path(out_dir, "14_expanding_window_cv.csv"), row.names = FALSE)

## ---- 15. Out-of-sample validation ----
cat("\n15. Running out-of-sample validation...\n")
oos_result <- out_of_sample_validation(final_huber_fit, group1_test, "Close")
generalizability_result <- generalizability_test(final_huber_fit, group2_valid_list, "Close", selected_features)
write.csv(oos_result, file.path(out_dir, "15_oos_validation_group1.csv"), row.names = FALSE)
write.csv(generalizability_result, file.path(out_dir, "15_generalizability_group2.csv"), row.names = FALSE)

## ---- 16. Model comparison ----
cat("\n16. Comparing models...\n")
tuned_hyperparams <- tune_all_models(group1_train, group1_valid, "Close", selected_features,
                                     fixed_huber_k = joint_best_k)
write.csv(hyperparams_to_table(tuned_hyperparams), file.path(out_dir, "16b_tuned_hyperparameters.csv"), row.names = FALSE)

model_comparison_valid_full <- compare_models(group1_train, group1_valid, "Close", selected_features,
                                              hyperparams = tuned_hyperparams)
model_comparison_valid <- model_comparison_valid_full$summary
write.csv(model_comparison_valid, file.path(out_dir, "16_model_comparison_VALIDATION.csv"), row.names = FALSE)
best_model_by_validation <- model_comparison_valid$Model[which.min(model_comparison_valid$MAE)]
message("Best model by Group 1 VALIDATION MAE: ", best_model_by_validation)

model_comparison_full <- compare_models(group1_train, group1_test, "Close", selected_features,
                                        hyperparams = tuned_hyperparams)
model_comparison <- model_comparison_full$summary
write.csv(model_comparison, file.path(out_dir, "16_model_comparison_TEST.csv"), row.names = FALSE)

## ---- 16b2. Per-coin RMSE, converted from [0,1]-normalized scale back to
## raw dollar prices. The pooled model_comparison table above reports ONE
## RMSE across all four Group1 coins' test rows combined -- not comparable
## to published single-coin studies (e.g. "BTC Prediction on test data"),
## and not in dollar units even if it were isolated per coin, since RMSE
## computed on min-max-normalized Close is not scale-invariant across
## different price levels. Both issues are fixed here: RMSE is computed
## SEPARATELY per coin (using group1_test$Symbol to subset each model's
## pooled predictions/actuals), then converted to dollars via
## RMSE_raw = RMSE_normalized * (Train_Close_Max - Train_Close_Min) for
## THAT coin specifically -- only the training-period RANGE matters for
## RMSE (a dispersion statistic), not the min itself, and it must be the
## TRAINING range since that is what normalization was actually fit on.
## NOTE: if a coin's test-period Close exceeds its training-period max (or
## falls below its training min), the model is being evaluated on
## extrapolated price levels it never saw during training -- flagged
## explicitly below rather than left implicit in the dollar figure.
rmse_by_coin_dollars <- bind_rows(lapply(names(model_comparison_full$predictions), function(model_name) {
  pred <- model_comparison_full$predictions[[model_name]]
  actual <- model_comparison_full$actual
  bind_rows(lapply(GROUP1, function(coin) {
    idx <- group1_test$Symbol == coin
    if (sum(idx) == 0) return(NULL)
    rmse_norm <- sqrt(mean((actual[idx] - pred[idx])^2, na.rm = TRUE))
    rng <- close_norm_params[[coin]]
    train_range <- unname(rng["max"] - rng["min"])
    test_close_actual_range <- range(group1_test$Close[idx] * train_range + unname(rng["min"]))
    extrapolated <- test_close_actual_range[1] < unname(rng["min"]) || test_close_actual_range[2] > unname(rng["max"])
    data.frame(Model = model_name, Symbol = coin,
               RMSE_normalized = rmse_norm,
               Train_Close_Range_USD = train_range,
               RMSE_raw_USD = rmse_norm * train_range,
               Test_Period_Extrapolates_Beyond_Training_Range = extrapolated)
  }))
}))
write.csv(rmse_by_coin_dollars, file.path(out_dir, "16_rmse_by_coin_raw_dollars.csv"), row.names = FALSE)
if (any(rmse_by_coin_dollars$Test_Period_Extrapolates_Beyond_Training_Range)) {
  message("NOTE: for coin(s) flagged TRUE above, the test period's actual price range falls outside ",
          "the training period's range -- dollar RMSE for those rows reflects performance partly in ",
          "an extrapolated price regime, not purely interpolation within the training range.")
}

## ---- 16c. Alternative feature sets ----
cat("\n16c. Comparing alternative feature sets...\n")
k_alt <- length(selected_features)
rf_features_alt <- rf_importance_selection(group1_train[, predictor_cols], group1_train$Close,
                                           predictor_cols, k = k_alt)
stepwise_features_alt <- forward_stepwise_selection(group1_train[, predictor_cols], group1_train$Close,
                                                    predictor_cols, k = k_alt)
bic_features_alt <- bic_stepwise_selection(group1_train[, predictor_cols], group1_train$Close,
                                           predictor_cols, k = k_alt)

# Genetic algorithm for feature selection
genetic_feature_selection <- function(X, y, feature_names,
                                      pop_size = 30, max_iter = 40,
                                      pc = 0.8, pm = 0.02, elitism = 2, tournament_k = 3,
                                      w1 = 0.7, w2 = 0.2, w3 = 0.1, seed = 99) {
  set.seed(seed)
  p <- length(feature_names)
  
  compute_fitness_ga <- function(S_idx, X, y, feature_names,
                                 w1 = 0.7, w2 = 0.2, w3 = 0.1,
                                 corr_threshold = 0.9, vif_threshold = 10) {
    if (length(S_idx) == 0) return(Inf)
    sel_names <- feature_names[S_idx]
    Xs <- X[, sel_names, drop = FALSE]
    dat <- data.frame(y = y, Xs)
    dat <- dat[complete.cases(dat), ]
    if (nrow(dat) < (ncol(dat) + 5)) return(Inf)
    
    fit <- tryCatch(lm(y ~ ., data = dat), error = function(e) NULL)
    if (is.null(fit)) return(Inf)
    mae <- mean(abs(residuals(fit)))
    
    complexity <- length(S_idx) / p
    
    vif_penalty <- 0
    if (length(S_idx) > 1) {
      vif_vals <- tryCatch(car::vif(fit), error = function(e) NULL)
      if (!is.null(vif_vals)) {
        max_vif <- max(vif_vals)
        if (max_vif > vif_threshold) vif_penalty <- max_vif - vif_threshold
      }
    }
    
    corr_penalty <- 0
    if (length(S_idx) > 1) {
      cmat <- suppressWarnings(cor(Xs, use = "pairwise.complete.obs"))
      diag(cmat) <- 0
      viol <- abs(cmat) - corr_threshold
      viol[viol < 0] <- 0
      corr_penalty <- sum(viol) / 2
    }
    
    w1 * mae + w2 * complexity + w3 * vif_penalty + corr_penalty
  }
  
  fitness_of <- function(bin_vec) {
    S_idx <- which(bin_vec == 1)
    compute_fitness_ga(S_idx, X, y, feature_names, w1, w2, w3)
  }
  
  pop <- matrix(rbinom(pop_size * p, 1, 0.3), nrow = pop_size, ncol = p)
  for (i in 1:pop_size) if (sum(pop[i, ]) == 0) pop[i, sample(p, 1)] <- 1
  
  fitness_vals <- apply(pop, 1, fitness_of)
  best_idx <- which.min(fitness_vals)
  best_solution <- pop[best_idx, ]
  best_fitness <- fitness_vals[best_idx]
  history <- numeric(max_iter)
  
  tournament_select <- function() {
    idx <- sample(seq_len(pop_size), tournament_k)
    idx[which.min(fitness_vals[idx])]
  }
  
  for (gen in 1:max_iter) {
    new_pop <- matrix(NA, nrow = pop_size, ncol = p)
    
    elite_idx <- order(fitness_vals)[seq_len(elitism)]
    new_pop[seq_len(elitism), ] <- pop[elite_idx, ]
    
    i <- elitism + 1
    while (i <= pop_size) {
      parent1 <- pop[tournament_select(), ]
      parent2 <- pop[tournament_select(), ]
      
      if (runif(1) < pc) {
        cut <- sample(seq_len(p - 1), 1)
        child <- c(parent1[seq_len(cut)], parent2[(cut + 1):p])
      } else {
        child <- parent1
      }
      
      mutate_idx <- runif(p) < pm
      child[mutate_idx] <- 1 - child[mutate_idx]
      if (sum(child) == 0) child[sample(p, 1)] <- 1
      
      new_pop[i, ] <- child
      i <- i + 1
    }
    
    pop <- new_pop
    fitness_vals <- apply(pop, 1, fitness_of)
    
    cur_best <- which.min(fitness_vals)
    if (fitness_vals[cur_best] < best_fitness) {
      best_fitness <- fitness_vals[cur_best]
      best_solution <- pop[cur_best, ]
    }
    history[gen] <- best_fitness
  }
  
  list(selected_features = feature_names[best_solution == 1],
       best_fitness = best_fitness, history = history)
}

ga_result_alt <- genetic_feature_selection(group1_train[, predictor_cols], group1_train$Close,
                                           predictor_cols, pop_size = 30, max_iter = 40)
ga_features_alt <- ga_result_alt$selected_features
message("Genetic algorithm selected ", length(ga_features_alt), " features (best fitness = ",
        round(ga_result_alt$best_fitness, 5), ")")

write.csv(data.frame(Feature = rf_features_alt), file.path(out_dir, "16c_rf_importance_features.csv"), row.names = FALSE)
write.csv(data.frame(Feature = stepwise_features_alt), file.path(out_dir, "16c_forward_stepwise_features.csv"), row.names = FALSE)
write.csv(data.frame(Feature = bic_features_alt), file.path(out_dir, "16c_bic_stepwise_features.csv"), row.names = FALSE)
write.csv(data.frame(Feature = ga_features_alt), file.path(out_dir, "16c_genetic_algorithm_features.csv"), row.names = FALSE)

# VIF check on all feature sets
vif_quantum_set <- vif_analysis(group1_train, selected_features, "Close")
vif_rf_set <- vif_analysis(group1_train, rf_features_alt, "Close")
vif_stepwise_set <- vif_analysis(group1_train, stepwise_features_alt, "Close")
vif_ga_set <- vif_analysis(group1_train, ga_features_alt, "Close")

vif_by_feature_set <- bind_rows(
  vif_quantum_set %>% mutate(FeatureSet = "Multi-Strategy (Huber-tailored)"),
  vif_rf_set %>% mutate(FeatureSet = "RF Importance (top-k)"),
  vif_stepwise_set %>% mutate(FeatureSet = "Forward Stepwise"),
  vif_ga_set %>% mutate(FeatureSet = "Genetic Algorithm")
) %>% relocate(FeatureSet)
write.csv(vif_by_feature_set, file.path(out_dir, "16c_vif_by_feature_set.csv"), row.names = FALSE)

vif_summary_by_set <- vif_by_feature_set %>%
  group_by(FeatureSet) %>%
  summarise(Max_VIF = if (all(is.na(VIF))) NA_real_ else max(VIF, na.rm = TRUE),
            Mean_VIF = if (all(is.na(VIF))) NA_real_ else mean(VIF, na.rm = TRUE),
            N_Aliased = sum(Aliased, na.rm = TRUE), N_VIF_gt_10 = sum(VIF > 10, na.rm = TRUE),
            N_Features = n(),
            .groups = "drop")
write.csv(vif_summary_by_set, file.path(out_dir, "16c_vif_summary_by_feature_set.csv"), row.names = FALSE)

if (any(is.na(vif_summary_by_set$Max_VIF))) {
  degenerate_sets <- vif_summary_by_set$FeatureSet[is.na(vif_summary_by_set$Max_VIF)]
  message("NOTE: VIF is undefined (NA) for feature set(s) with fewer than 2 usable predictors: ",
          paste(degenerate_sets, collapse = ", "),
          " -- check the corresponding *_features.csv to confirm how many features were actually selected.")
}

comparison_huber_features <- model_comparison %>% mutate(FeatureSet = "Multi-Strategy (Huber-tailored)")
comparison_rf_features <- compare_models(group1_train, group1_test, "Close", rf_features_alt,
                                         hyperparams = tuned_hyperparams)$summary %>%
  mutate(FeatureSet = "RF Importance (top-k)")
comparison_stepwise_features <- compare_models(group1_train, group1_test, "Close", stepwise_features_alt,
                                               hyperparams = tuned_hyperparams)$summary %>%
  mutate(FeatureSet = "Forward Stepwise")
comparison_ga_features <- compare_models(group1_train, group1_test, "Close", ga_features_alt,
                                         hyperparams = tuned_hyperparams)$summary %>%
  mutate(FeatureSet = "Genetic Algorithm")

feature_set_comparison <- bind_rows(comparison_huber_features, comparison_rf_features,
                                    comparison_stepwise_features, comparison_ga_features) %>%
  relocate(FeatureSet, .after = Model) %>%
  arrange(Model, FeatureSet)
write.csv(feature_set_comparison, file.path(out_dir, "16c_model_comparison_by_feature_set.csv"), row.names = FALSE)

feature_set_ranks <- feature_set_comparison %>%
  group_by(FeatureSet) %>%
  mutate(MAE_Rank = rank(MAE)) %>%
  ungroup() %>%
  arrange(FeatureSet, MAE_Rank)
write.csv(feature_set_ranks, file.path(out_dir, "16c_model_ranks_by_feature_set.csv"), row.names = FALSE)

huber_rank_by_set <- feature_set_ranks %>% filter(Model == "Huber") %>% select(FeatureSet, MAE, MAE_Rank)
message("Huber's MAE rank across feature sets (1 = best of 9):")
for (i in seq_len(nrow(huber_rank_by_set))) {
  message("  ", huber_rank_by_set$FeatureSet[i], ": rank ", huber_rank_by_set$MAE_Rank[i],
          " (MAE = ", round(huber_rank_by_set$MAE[i], 5), ")")
}
write.csv(huber_rank_by_set, file.path(out_dir, "16c_huber_rank_by_feature_set.csv"), row.names = FALSE)

## ---- 17. Statistical properties ----
cat("\n17. Computing feature statistical properties...\n")
vif_selected <- vif_analysis(group1_train, selected_features, "Close")
wald_table <- wald_test_bootstrap(final_huber_fit, boot_res$boot_coefs)
feature_props <- feature_statistical_properties(selected_features, group1_train, "Close",
                                                vif_selected, stability_res$pi_hat, wald_table)
write.csv(feature_props, file.path(out_dir, "17_feature_statistical_properties.csv"), row.names = FALSE)

## ---- 18. Paired t-tests ----
cat("\n18. Running paired t-tests...\n")
paired_tests <- paired_model_ttests(model_comparison_full$predictions, model_comparison_full$actual,
                                    baseline = "Huber")
write.csv(paired_tests, file.path(out_dir, "18_paired_ttests.csv"), row.names = FALSE)

## ---- 19. Actual vs predicted ----
cat("\n19. Generating actual vs predicted plots...\n")
png(file.path(out_dir, "19_actual_vs_predicted_group1.png"), width = 900, height = 900)
plot_actual_vs_predicted_group1(final_huber_fit, split_list, selected_features, response = "Close")
dev.off()

group1_bias <- per_coin_bias(final_huber_fit, split_list, selected_features, coins = GROUP1, split_name = "test")
write.csv(group1_bias, file.path(out_dir, "19_per_coin_bias.csv"), row.names = FALSE)

png(file.path(out_dir, "19_prediction_error_distribution.png"), width = 800, height = 600)
error_dist_bounds <- prediction_error_distribution(final_huber_fit, group1_test, selected_features, response = "Close")
dev.off()
write.csv(error_dist_bounds, file.path(out_dir, "19_prediction_error_bounds.csv"), row.names = FALSE)

## ---- 20. Bootstrap prediction intervals ----
cat("\n20. Computing bootstrap prediction intervals (B=1000)...\n")
boot_pred_result <- bootstrap_prediction_intervals(group1_train, group1_test, "Close", selected_features,
                                                   B = 1000, block_length = 10, calib_d = group1_valid)

scale_info <- conformalize_prediction_intervals(boot_pred_result, group1_valid, group1_test, response = "Close")
message("Conformal calibration scale factors by coin:")
for (nm in names(scale_info$scale_by_coin)) {
  message("  ", nm, ": ", round(scale_info$scale_by_coin[nm], 3))
}
write.csv(data.frame(Symbol = names(scale_info$scale_by_coin),
                      ConformalScaleFactor = as.numeric(scale_info$scale_by_coin)),
          file.path(out_dir, "20_conformal_scale_factors.csv"), row.names = FALSE)
boot_pred_calibrated <- apply_conformal_scaling(boot_pred_result, group1_test, scale_info)

png(file.path(out_dir, "20_bootstrap_prediction_ci.png"), width = 1200, height = 900)
plot_bootstrap_prediction_ci_group1(group1_test, boot_pred_calibrated, group1 = GROUP1, response = "Close")
dev.off()

bootstrap_coef_table <- bootstrap_coefficient_table(final_huber_fit, boot_res$boot_coefs)
write.csv(bootstrap_coef_table, file.path(out_dir, "20_bootstrap_coefficient_table.csv"), row.names = FALSE)

ci_summary_raw <- bootstrap_ci_summary(group1_test, boot_pred_result, response = "Close")
ci_summary <- bootstrap_ci_summary(group1_test, boot_pred_calibrated, response = "Close")
write.csv(ci_summary_raw$overall, file.path(out_dir, "20_bootstrap_ci_summary_RAW.csv"), row.names = FALSE)
write.csv(ci_summary_raw$by_coin, file.path(out_dir, "20_bootstrap_ci_summary_by_coin_RAW.csv"), row.names = FALSE)
write.csv(ci_summary$overall, file.path(out_dir, "20_bootstrap_ci_summary.csv"), row.names = FALSE)
write.csv(ci_summary$by_coin, file.path(out_dir, "20_bootstrap_ci_summary_by_coin.csv"), row.names = FALSE)

## ---- 21. Feature-selection stability comparison ----
cat("\n21. Comparing feature selection stability...\n")
stability_comparison <- compare_feature_selection_stability(
  X = group1_train[, predictor_cols], y = group1_train$Close, feature_names = predictor_cols,
  B = 15, pi_thr = 0.6, k = length(selected_features), quantum_stability_result = stability_res)


write.csv(stability_comparison, file.path(out_dir, "21_stability_comparison.csv"), row.names = FALSE)

## ---- 21b. Consolidated feature-selection method comparison ----
cat("\n21b. Creating feature selection method comparison table...\n")
lasso_features_final <- tryCatch(
  lasso_cv_selection(group1_train[, predictor_cols], group1_train$Close, predictor_cols),
  error = function(e) character(0))
if (length(lasso_features_final) == 0) lasso_features_final <- predictor_cols[1]

rfe_features_final <- tryCatch(
  rfe_selection(group1_train[, predictor_cols], group1_train$Close, predictor_cols, size = k_alt),
  error = function(e) character(0))
if (length(rfe_features_final) == 0) rfe_features_final <- predictor_cols[1]

eval_huber_test_mae <- function(features) {
  fit <- tryCatch(fit_huber(group1_train, "Close", features, k = tuned_hyperparams$Huber$k),
                  error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  pred <- predict(fit, newdata = group1_test)
  mean(abs(group1_test$Close - pred), na.rm = TRUE)
}
lasso_mae_final <- eval_huber_test_mae(lasso_features_final)
rfe_mae_final <- eval_huber_test_mae(rfe_features_final)

scad_features_final <- tryCatch(
  scad_selection(group1_train[, predictor_cols], group1_train$Close, predictor_cols),
  error = function(e) character(0))
if (length(scad_features_final) == 0) scad_features_final <- predictor_cols[1]
scad_mae_final <- eval_huber_test_mae(scad_features_final)

adalasso_features_final <- tryCatch(
  adaptive_lasso_selection(group1_train[, predictor_cols], group1_train$Close, predictor_cols),
  error = function(e) character(0))
if (length(adalasso_features_final) == 0) adalasso_features_final <- predictor_cols[1]
adalasso_mae_final <- eval_huber_test_mae(adalasso_features_final)

get_fs_mae <- function(feature_set_label) {
  val <- feature_set_comparison$MAE[feature_set_comparison$Model == "Huber" &
                                      feature_set_comparison$FeatureSet == feature_set_label]
  if (length(val) == 0) return(NA_real_)
  val
}

stability_lookup <- setNames(stability_comparison$Avg_Jaccard_Stability, stability_comparison$Method)

# NOTE: lookup keys below must match stability_comparison$Method /
# feature_set_comparison$FeatureSet EXACTLY (case- and suffix-sensitive --
# e.g. "Forward Stepwise (AIC)" here vs. bare "Forward Stepwise" in
# feature_set_comparison, since those two tables use different naming
# conventions for the same method). A mismatch here previously produced a
# fabricated 0.0 (masked via silent NA-replacement) rather than a visible
# failure -- see the warning checks after these vectors are built.
method_names <- c("Optimization-based (Ours)", "LASSO (CV-selected lambda)", 
                  "Recursive Feature Elimination", "Random Forest Importance", 
                  "Forward Stepwise (AIC)", "SCAD", "Adaptive LASSO")

n_features <- c(length(selected_features), length(lasso_features_final), 
                length(rfe_features_final), length(rf_features_alt), 
                length(stepwise_features_alt), length(scad_features_final),
                length(adalasso_features_final))

stability_vals <- c(stability_lookup["Multi-Strategy Integration"], 
                    stability_lookup["LASSO (CV)"],
                    stability_lookup["RFE"], 
                    stability_lookup["RandomForest Importance"],
                    stability_lookup["Forward Stepwise (AIC)"],
                    stability_lookup["SCAD"],
                    stability_lookup["Adaptive LASSO"])

mae_vals <- c(get_fs_mae("Multi-Strategy (Huber-tailored)"), 
              lasso_mae_final, rfe_mae_final,
              get_fs_mae("RF Importance (top-k)"), 
              get_fs_mae("Forward Stepwise"),
              scad_mae_final, adalasso_mae_final)

# Surface, don't mask, any lookup failure -- a silently-substituted 0 is
# indistinguishable from a genuine zero-stability result (0 is a valid
# point on the Jaccard scale), and an Inf substituted for a failed MAE
# lookup can silently corrupt any downstream min()/comparison over this
# column. Print exactly which method(s) failed to match and why, so a
# naming mismatch between stability_comparison$Method and this table's
# lookup keys (the actual bug that produced Forward Stepwise's phantom
# 0.0 in this table) is caught immediately rather than only discoverable
# by manually cross-checking two separate output tables against each other.
if (any(is.na(stability_vals))) {
  message("WARNING: stability lookup failed for: ", paste(method_names[is.na(stability_vals)], collapse = ", "),
          " -- check that stability_comparison$Method's labels match the keys used above.")
}
if (any(is.na(mae_vals))) {
  message("WARNING: MAE lookup failed for: ", paste(method_names[is.na(mae_vals)], collapse = ", "),
          " -- check that feature_set_comparison$FeatureSet's labels match the keys used above.")
}

stability_se_lookup <- setNames(stability_comparison$Jaccard_SE, stability_comparison$Method)
stability_lo_lookup <- setNames(stability_comparison$Jaccard_CI_Lower, stability_comparison$Method)
stability_hi_lookup <- setNames(stability_comparison$Jaccard_CI_Upper, stability_comparison$Method)
stability_lookup_names <- c("Multi-Strategy Integration", "LASSO (CV)", "RFE", "RandomForest Importance",
                            "Forward Stepwise (AIC)", "SCAD", "Adaptive LASSO")
stability_se_vals <- stability_se_lookup[stability_lookup_names]
stability_lo_vals <- stability_lo_lookup[stability_lookup_names]
stability_hi_vals <- stability_hi_lookup[stability_lookup_names]

feature_selection_summary_table <- data.frame(
  Method = method_names,
  N_Features_Selected = n_features,
  Mean_Jaccard_Stability = as.numeric(stability_vals),
  Jaccard_SE = as.numeric(stability_se_vals),
  Jaccard_CI_Lower = as.numeric(stability_lo_vals),
  Jaccard_CI_Upper = as.numeric(stability_hi_vals),
  Mean_Test_MAE = as.numeric(mae_vals),
  row.names = NULL  # Explicitly set row.names to NULL
)

write.csv(feature_selection_summary_table, file.path(out_dir, "21b_feature_selection_method_comparison.csv"),
          row.names = FALSE)

## ---- 22. Group 2: universal vs. asset-specific ----
cat("\n22. Testing universal vs. asset-specific features on Group 2...\n")
universal_vs_specific <- compare_universal_vs_asset_specific(final_huber_fit, split_list, predictor_cols,
                                                             selected_features, group2 = GROUP2, response = "Close")
hypothesis_result <- hypothesis_test_universal_vs_specific(universal_vs_specific)
write.csv(universal_vs_specific, file.path(out_dir, "22_universal_vs_assetspecific.csv"), row.names = FALSE)
write.csv(hypothesis_result, file.path(out_dir, "22_hypothesis_test.csv"), row.names = FALSE)

png(file.path(out_dir, "22_actual_vs_predicted_group2.png"), width = 900, height = 900)
plot_actual_vs_predicted_coins(final_huber_fit, split_list, selected_features, coins = GROUP2,
                               split_name = "valid", response = "Close")
dev.off()

## ---- Save full results ----
cat("\nSaving full results bundle...\n")
saveRDS(list(desc_stats = desc_stats_table, temporal_props = temporal_props_table,
             correlation_results = correlation_results, vif_results = vif_results,
             pca_results = pca_results, qfs_result = qfs_result,
             selected_features = selected_features, final_huber_fit = final_huber_fit,
             validation_check = validation_check,
             residual_diagnostics = res_diag, residual_diagnostics_table = res_diag_table, bootstrap = boot_res,
             stability = stability_res, cv10 = cv10, expanding_cv = expanding_cv,
             oos_result = oos_result, generalizability_result = generalizability_result,
             model_comparison_validation = model_comparison_valid,
             best_model_by_validation = best_model_by_validation,
             tuned_hyperparams = tuned_hyperparams,
             feature_set_comparison = feature_set_comparison, huber_rank_by_feature_set = huber_rank_by_set,
             vif_by_feature_set = vif_by_feature_set, vif_summary_by_feature_set = vif_summary_by_set,
             model_comparison = model_comparison, model_predictions = model_comparison_full$predictions,
             close_norm_params = close_norm_params, rmse_by_coin_dollars = rmse_by_coin_dollars,
             feature_statistical_properties = feature_props, paired_ttests = paired_tests,
             prediction_error_bounds = error_dist_bounds, group1_per_coin_bias = group1_bias,
             bootstrap_prediction_ci = boot_pred_result, bootstrap_prediction_ci_calibrated = boot_pred_calibrated,
             conformal_scale_factors = scale_info, bootstrap_ci_summary_raw = ci_summary_raw,
             bootstrap_ci_summary = ci_summary,
             bootstrap_coefficient_table = bootstrap_coef_table,
             stability_comparison = stability_comparison,
             feature_selection_summary_table = feature_selection_summary_table,
             universal_vs_asset_specific = universal_vs_specific,
             hypothesis_test_result = hypothesis_result),
        file.path(out_dir, "pipeline_results.rds"))

cat("\n================ PIPELINE COMPLETE ================\n")
cat("Selected meta-features (", length(selected_features), "):\n", sep = "")
print(selected_features)
cat("\nGroup 1 VALIDATION check (selected features/model, never used for fitting):\n")
print(validation_check)
cat("\nModel comparison on Group 1 VALIDATION (used for model selection):\n")
print(model_comparison_valid)
cat("\nBest model by VALIDATION MAE:", best_model_by_validation, "\n")
cat("\nModel comparison on Group 1 TEST (final unbiased reporting, sorted by MAE):\n")
print(model_comparison)
cat("\nSupplementary: Huber's MAE rank across feature sets (1 = best of 9):\n")
print(huber_rank_by_set)
cat("\nSupplementary: VIF summary by feature set (collinearity check):\n")
print(vif_summary_by_set)
cat("\nPaired t-tests (Huber vs. each benchmark):\n")
print(paired_tests)
cat("\nFeature-selection stability comparison:\n")
print(stability_comparison)
cat("\nUniversal vs. asset-specific hypothesis test (Group 2):\n")
print(hypothesis_result)
cat("=====================================================\n")


