# Cryptocurrency Price Prediction: Joint Feature Selection & Robust Regression

A statistical framework for cryptocurrency closing-price prediction that jointly optimizes feature selection and Huber M-estimation hyperparameters within a single stochastic search
(Multi-Strategy Integration: opposition-based learning + genetic operators + guided search + quantum-inspired diversification).

## Repository contents

| File | Purpose |
|---|---|
| `crypto_feature_engineering.R` | Downloads daily OHLCV data for 8 cryptocurrencies from Yahoo Finance, engineers 66 candidate features across 8 categories (price-based, volatility, momentum, volume, macro/cross-asset, etc.), and writes `all_coins_features.csv`. |
| `crypto_analysis_pipeline_M.R` | Full analysis pipeline: descriptive statistics, temporal diagnostics, correlation/VIF/PCA, the joint feature-selection + hyperparameter search, model comparison against 8 benchmarks, bootstrap uncertainty quantification, feature-selection stability analysis, and out-of-sample validation on unseen coins. |

## Requirements

- **R** ≥ 4.0 (tested on R 4.x)
- **Internet access** for the first script (Yahoo Finance data download)
- No manual package installation needed — both scripts auto-install any missing CRAN packages
  on first run via a `required_packages` check at the top of each file.

Key packages used: `quantmod`, `TTR`, `dplyr`, `MASS` (Huber M-estimation), `glmnet` (LASSO/adaptive
LASSO), `ncvreg` (SCAD), `caret` (RFE), `randomForest`, `ranger`, `gbm`, `xgboost`, `e1071` (SVM),
`rBayesianOptimization`, `car`, `lmtest`, `tseries`.

## How to run

Clone the repository, then run the two scripts **in order** from the repository root:

```bash
git clone https://github.com/anuajayi/crypto-price-prediction.git
cd crypto-price-prediction

Rscript crypto_feature_engineering.R
Rscript crypto_analysis_pipeline_M.R
```

`crypto_feature_engineering.R` must run first — the analysis pipeline reads its output
(`all_coins_features.csv`) as input.

## What gets produced

Running the analysis pipeline creates a `crypto_output/` directory containing:

- Numbered CSVs for every analysis stage (descriptive stats, temporal properties,
  correlation/VIF/PCA per coin, feature selection results, model comparison tables, bootstrap
  coefficient/prediction-interval tables, stability comparisons, generalizability results)
- PNG plots (residual diagnostics, actual-vs-predicted scatter, bootstrap prediction
  intervals, GA convergence curves)
- pipeline_results.rds — every intermediate and final object from the run, for loading
  back into an R session (`readRDS("crypto_output/pipeline_results.rds")`) without re-running
  the full pipeline.

None of these generated files are tracked in version control (see `.gitignore`) — they're
regenerated fresh each run.

## Runtime note

The pipeline includes a 20-run convergence diagnostic and a 15-configuration sensitivity
analysis for the joint feature-selection search, both of which run the search's full
population-based optimization multiple times. Expect this to be the most time-consuming part
of a full run; reduce `N_RUNS_MULTISTRATEGY` or the sensitivity analysis's parameter grids near
the top of Section 9 in `crypto_analysis_pipeline_M.R` if runtime is a constraint.

## Methodology summary

- Data: 8 cryptocurrencies (BTC, ETH, BCH, DASH, XLM, XRP, XMR, LTC), split into Group 1
  (model development) and Group 2 (out-of-sample generalizability testing on unseen assets).
- Feature selection: a stochastic Multi-Strategy Integration search jointly selects a
  sparse feature subset and Huber's robustness hyperparameter (k) in a single optimization,
  rather than fixing one and tuning the other sequentially.
- Validation: residual diagnostics (6 tests, with multiple-testing correction), bootstrap
  coefficient and prediction intervals, k-fold and expanding-window cross-validation,
  feature-selection stability analysis (with bootstrap standard errors) against 7 alternative
  methods (LASSO, adaptive LASSO, SCAD, RFE, Random Forest importance, forward stepwise AIC/BIC).
- Comparison: 8 benchmark models (AdaBoost, Extra Trees, GBM, XGBoost, SVM, Random Forest,
  and voting/stacking ensembles), all with Bayesian-optimized hyperparameters.
