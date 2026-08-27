#############################################################################
# Cryptocurrency Data Collection & Feature Engineering
#
# Coins:   BTC, ETH, BCH, DASH, XLM, XRP, XMR, LTC
# Source:  Yahoo Finance (via quantmod::getSymbols)
# Period:  2020-08-20  to  2026-06-30
#
# For each coin, extracts OHLCV (+ Date) and engineers 62 additional features
# across 8 categories:
#   Category 1: Market data features         -  4 (satisfied by the raw
#                                                Open/High/Low/Volume columns already present
#   Category 2: Temporal features            -  4 variables
#   Category 3: Price-based features         - 13 variables
#   Category 4: Volatility features          -  8 variables
#   Category 5: Momentum indicators          - 10 variables
#   Category 6: Volume-based features        -  8 variables
#   Category 7: Market microstructure feat.  -  7 variables
#   Category 8: Macro / cross-asset feat.    - 12 variables (6 series x level + daily return)
#   -----------------------------------------------------
#   TOTAL NEW ENGINEERED FEATURES            - 66 variables
#
# Category 1 ("market data") satisfied directly by the raw Open/High/Low/Close/Volume columns 
#
# Category 8 (macro/cross-asset) pulls Gold (GC=F), Crude Oil WTI (CL=F),
# S&P 500 (^GSPC), EUR/USD (EURUSD=X), GBP/USD (GBPUSD=X), and Nikkei 225
# (^N225) from Yahoo Finance. Crypto trades 7 days/week; these markets
# don't, so each series is forward-filled across non-trading days onto a
# full daily calendar, then both its level and daily return are lagged by
# 1 day before merging into each coin's feature set -- ensuring only
# already-observed macro information is used (avoids look-ahead bias from
# markets that close at different times/timezones than crypto's daily bar).
#############################################################################

# ---------------------------------------------------------------------------
# 0. Setup
# ---------------------------------------------------------------------------
required_packages <- c("quantmod", "TTR", "dplyr", "tidyr", "lubridate", "zoo")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "https://cloud.r-project.org")

library(quantmod)
library(TTR)
library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)

# ---------------------------------------------------------------------------
# 1. Configuration
# ---------------------------------------------------------------------------
symbols    <- c("BTC-USD", "ETH-USD", "BCH-USD", "DASH-USD",
                 "XLM-USD", "XRP-USD", "XMR-USD", "LTC-USD")
coin_names <- c("BTC", "ETH", "BCH", "DASH", "XLM", "XRP", "XMR", "LTC")

start_date <- as.Date("2020-08-20")
end_date   <- as.Date("2026-06-30")

out_dir <- "crypto_output"
if (!dir.exists(out_dir)) dir.create(out_dir)

# ---------------------------------------------------------------------------
# 2. Data collection
# ---------------------------------------------------------------------------
get_coin_data <- function(sym) {
  message("Downloading: ", sym)
  dat <- tryCatch(
    getSymbols(sym, src = "yahoo", from = start_date, to = end_date,
               auto.assign = FALSE),
    error = function(e) {
      message("  Failed to download ", sym, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(dat) || nrow(dat) == 0) return(NULL)

  df <- data.frame(Date = index(dat), coredata(dat))
  colnames(df) <- c("Date", "Open", "High", "Low", "Close", "Volume", "Adjusted")
  df %>% arrange(Date) %>% distinct(Date, .keep_all = TRUE)
}

raw_data_list <- setNames(lapply(symbols, get_coin_data), coin_names)
raw_data_list <- raw_data_list[!sapply(raw_data_list, is.null)]

if (length(raw_data_list) == 0) stop("No data was downloaded for any symbol.")

# ---------------------------------------------------------------------------
# 2b. Macro / cross-asset indicators
#     Gold, Crude Oil (WTI), S&P 500, EUR/USD, GBP/USD, Nikkei 225.
#     Crypto trades 7 days/week; these markets don't -- so each series is
#     forward-filled across non-trading days onto a full daily calendar.
#     Both the level and the daily return are then lagged by 1 day so that
#     a given crypto trading day only ever uses already-observed macro
#     information (these markets close hours before/after crypto's daily
#     bar in different timezones, so same-day values would leak information
#     crypto markets hadn't seen yet).
# ---------------------------------------------------------------------------
macro_symbols <- c("GC=F", "CL=F", "^GSPC", "EURUSD=X", "GBPUSD=X", "^N225")
macro_names   <- c("Gold", "CrudeOil", "SP500", "EURUSD", "GBPUSD", "Nikkei")

get_macro_series <- function(sym) {
  message("Downloading macro series: ", sym)
  dat <- tryCatch(
    getSymbols(sym, src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE),
    error = function(e) {
      message("  Failed to download ", sym, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  data.frame(Date = as.Date(index(dat)), Close = as.numeric(quantmod::Cl(dat))) %>%
    dplyr::arrange(Date) %>% dplyr::distinct(Date, .keep_all = TRUE)
}

full_calendar <- data.frame(Date = seq(start_date, end_date, by = "day"))

prepare_macro_series <- function(sym, out_name, full_calendar) {
  raw <- get_macro_series(sym)
  if (is.null(raw)) {
    message("  Skipping macro series ", out_name, " (download failed)")
    out <- full_calendar
    out[[out_name]] <- NA_real_
    out[[paste0(out_name, "_Return")]] <- NA_real_
    return(out)
  }
  merged <- dplyr::left_join(full_calendar, raw, by = "Date")
  # Forward-fill across weekends/holidays (macro market closed that day);
  # back-fill any leading gap before the macro series' first trading day.
  merged$Close <- zoo::na.locf(merged$Close, na.rm = FALSE)
  merged$Close <- zoo::na.locf(merged$Close, fromLast = TRUE, na.rm = FALSE)

  ret <- merged$Close / dplyr::lag(merged$Close) - 1

  # Lag both level and return by 1 day -- today's crypto features use
  # yesterday's (already-observed) macro close/return, not same-day.
  merged[[out_name]] <- dplyr::lag(merged$Close, 1)
  merged[[paste0(out_name, "_Return")]] <- dplyr::lag(ret, 1)

  merged[, c("Date", out_name, paste0(out_name, "_Return"))]
}

macro_features <- Reduce(
  function(x, y) dplyr::full_join(x, y, by = "Date"),
  Map(prepare_macro_series, macro_symbols, macro_names, MoreArgs = list(full_calendar = full_calendar))
)

# ---------------------------------------------------------------------------
# 3. Feature engineering (50 new features / 7 categories, + 12 macro
#    cross-asset features merged in below)
# ---------------------------------------------------------------------------
engineer_features <- function(df) {
  df <- df %>% arrange(Date)

  o <- df$Open; h <- df$High; l <- df$Low; c <- df$Close; v <- df$Volume

  ## ---- Category 1: Market data features (0 new) -----------------------
  ## Satisfied directly by the raw Open/High/Low/Close/Volume columns already present in df.

  ## ---- Category 2: Temporal features (4) ------------------------------
  df$Year      <- year(df$Date)
  df$Month     <- month(df$Date)
  df$DayOfWeek <- wday(df$Date)     # 1 = Sunday ... 7 = Saturday
  df$Quarter   <- quarter(df$Date)

  ## ---- Category 3: Price-based features (13) ---------------------------
  df$Return         <- c / lag(c) - 1
  df$LogReturn      <- log(c / lag(c))
  df$PriceChange    <- c - o
  df$HL_Range        <- h - l
  df$CO_Ratio        <- c / o
  df$SMA_5           <- SMA(c, n = 5)
  df$SMA_10          <- SMA(c, n = 10)
  df$SMA_20          <- SMA(c, n = 20)
  df$SMA_50          <- SMA(c, n = 50)
  df$EMA_12          <- EMA(c, n = 12)
  df$EMA_26          <- EMA(c, n = 26)
  df$Price_to_SMA20  <- c / df$SMA_20
  df$TypicalPrice    <- (h + l + c) / 3

  ## ---- Category 4: Volatility features (8) ------------------------------
  df$RollSD_5  <- rollapply(df$Return, 5,  sd, fill = NA, align = "right")
  df$RollSD_10 <- rollapply(df$Return, 10, sd, fill = NA, align = "right")
  df$RollSD_20 <- rollapply(df$Return, 20, sd, fill = NA, align = "right")

  atr <- ATR(cbind(High = h, Low = l, Close = c), n = 14)
  df$ATR_14 <- atr[, "atr"]

  bb <- BBands(cbind(High = h, Low = l, Close = c), n = 20, sd = 2)
  df$BB_Upper <- bb[, "up"]
  df$BB_Lower <- bb[, "dn"]
  df$BB_Width <- (bb[, "up"] - bb[, "dn"]) / bb[, "mavg"]

  ohlc_mat <- cbind(Open = o, High = h, Low = l, Close = c)
  park <- tryCatch(
    volatility(ohlc_mat, n = 10, calc = "parkinson"),
    error = function(e) rep(NA_real_, length(c))
  )
  df$ParkinsonVol <- as.numeric(park)

  ## ---- Category 5: Momentum indicators (10) ------------------------------
  df$RSI_14 <- RSI(c, n = 14)

  macd <- MACD(c, nFast = 12, nSlow = 26, nSig = 9)
  df$MACD        <- macd[, "macd"]
  df$MACD_Signal <- macd[, "signal"]
  df$MACD_Hist   <- macd[, "macd"] - macd[, "signal"]

  stoch_vals <- stoch(cbind(High = h, Low = l, Close = c),
                       nFastK = 14, nFastD = 3, nSlowD = 3)
  df$Stoch_K <- stoch_vals[, "fastK"]
  df$Stoch_D <- stoch_vals[, "slowD"]

  df$ROC_10      <- ROC(c, n = 10, type = "discrete")
  df$Momentum_10 <- momentum(c, n = 10)
  df$WilliamsR   <- WPR(cbind(High = h, Low = l, Close = c), n = 14)
  df$CCI         <- CCI(cbind(High = h, Low = l, Close = c), n = 20)

  ## ---- Category 6: Volume-based features (8) -----------------------------
  df$Volume_SMA_5  <- SMA(v, n = 5)
  df$Volume_SMA_20 <- SMA(v, n = 20)
  df$Volume_Ratio  <- v / df$Volume_SMA_20
  df$OBV           <- OBV(c, v)

  vpt <- numeric(length(c))
  vpt[1] <- 0
  for (i in 2:length(c)) {
    ret <- (c[i] - c[i - 1]) / c[i - 1]
    vpt[i] <- vpt[i - 1] + v[i] * ret
  }
  df$VPT <- vpt

  mfm <- ((c - l) - (h - c)) / (h - l)          # money flow multiplier
  mfm[is.nan(mfm) | is.infinite(mfm)] <- 0
  mfv <- mfm * v
  df$CMF <- rollapply(mfv, 20, sum, fill = NA, align = "right") /
            rollapply(v,   20, sum, fill = NA, align = "right")

  df$Volume_ROC   <- ROC(v, n = 10, type = "discrete")
  df$DollarVolume <- c * v

  ## ---- Category 7: Market microstructure features (7) --------------------
  df$HL_Spread <- (h - l) / c
  df$Amihud    <- abs(df$Return) / df$DollarVolume

  price_diff <- c(NA, diff(c))
  df$RollSpread <- rollapply(price_diff, width = 20, align = "right", fill = NA,
    FUN = function(x) {
      if (length(x) < 3 || any(is.na(x))) return(NA_real_)
      cv <- cov(x[-1], x[-length(x)])
      if (cv < 0) 2 * sqrt(-cv) else NA_real_
    })

  df$PriceImpact <- abs(df$Return) / v
  df$KyleLambda  <- abs(df$PriceChange) / v
  df$IntradayVol <- (h - l) / o
  df$CLV         <- ((c - l) - (h - c)) / (h - l)

  df <- df %>% dplyr::select(-Adjusted)
  return(df)
}

# ---------------------------------------------------------------------------
# 3b. Missing-data treatment for lagged / warm-up technical indicators
#     (e.g. a 26-day EMA has no value for the first 25 rows). We forward-fill
#     using the first available value to cover that warm-up window, then drop
#     any rows that still contain NA (covers indicators outside the regex
#     below, or genuine source-data gaps) so saved CSVs are NA-free.
# ---------------------------------------------------------------------------
treat_missing_lagged_features <- function(df) {
  # Ratio-style features (CO_Ratio, HL_Spread, Amihud, PriceImpact,
  # KyleLambda, ...) divide by price or volume and produce +/-Inf on
  # zero-volume or zero-open-price days -- more common for smaller-cap
  # coins. is.na(Inf) is FALSE in R, so Inf silently passes through
  # complete.cases() and later corrupts Min-Max normalization (an Inf in a
  # column makes max = Inf, collapsing the whole column). Convert Inf/-Inf
  # to NA first so it's treated like any other missing value below.
  num_cols <- names(df)[sapply(df, is.numeric)]
  for (cn in num_cols) {
    x <- df[[cn]]
    x[is.infinite(x)] <- NA
    df[[cn]] <- x
  }

  lag_feature_cols <- names(df)[grepl(
    "SMA|EMA|RSI|MACD|ATR|BB_|ROC|Momentum|Stoch|CCI|WilliamsR|CMF|RollSD|
     ParkinsonVol|RollSpread|OBV|VPT|Return|LogReturn|Price_to_SMA20",
    names(df))]
  for (cn in lag_feature_cols) {
    x <- df[[cn]]
    x <- zoo::na.locf(x, na.rm = FALSE)                    # forward-fill
    x <- zoo::na.locf(x, fromLast = TRUE, na.rm = FALSE)    # cover leading NAs
    df[[cn]] <- x
  }
  # Ensure rows with any remaining NAs (including former Inf values, and
  # non-lag features like CO_Ratio/Amihud/etc. that had a genuine Inf on a
  # zero-volume/zero-price day) are deleted.
  df[complete.cases(df), ]
}

# ---------------------------------------------------------------------------
# 4. Apply feature engineering to every coin
# ---------------------------------------------------------------------------
feature_data_list <- lapply(names(raw_data_list), function(nm) {
  message("Engineering features for: ", nm)
  df <- engineer_features(raw_data_list[[nm]])
  df <- dplyr::left_join(df, macro_features, by = "Date")   # merge in 12 macro/cross-asset features
  df <- treat_missing_lagged_features(df)
  df$Symbol <- nm
  df %>% relocate(Symbol, .after = Date)
})
names(feature_data_list) <- names(raw_data_list)

# Save one CSV per coin
for (nm in names(feature_data_list)) {
  fp <- file.path(out_dir, paste0(nm, "_features.csv"))
  write.csv(feature_data_list[[nm]], fp, row.names = FALSE)
  message("Saved: ", fp)
}

# Combine into a single long-format master dataset
master_df <- bind_rows(feature_data_list)
write.csv(master_df, file.path(out_dir, "all_coins_features.csv"), row.names = FALSE)
message("Saved combined dataset: ", file.path(out_dir, "all_coins_features.csv"))

# ---------------------------------------------------------------------------
# 5. Quick summary
# ---------------------------------------------------------------------------
n_features <- ncol(master_df) - 2   # exclude Date, Symbol
cat("\n================ SUMMARY ================\n")
cat("Coins collected      :", paste(names(feature_data_list), collapse = ", "), "\n")
cat("Date range requested  :", as.character(start_date), "to", as.character(end_date), "\n")
cat("Total rows (all coins):", nrow(master_df), "\n")
cat("Total columns         :", ncol(master_df), "\n")
cat("Engineered features   :", n_features, "(target: 50 new + 12 macro + 5 raw OHLCV + Date/Symbol)\n")
cat("===========================================\n")
