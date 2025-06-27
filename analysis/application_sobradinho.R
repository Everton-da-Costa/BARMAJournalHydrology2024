# =============================================================================
# Application Sobradinho data set
# =============================================================================
# PAPER: Dynamic beta modeling of seasonal hydro-environmental time series
# Test inferences and link function selection in dynamic beta modeling of
# seasonal hydro-environmental time series with temporary abnormal regimes
#
# DESCRIPTION: The Sobradinho data set represents the useful volume of
# Sobradinho reservoirs. This code is used to analyze the data as described in
# the associated paper.
#
# AUTHOR: Everton da Costa
# UPDATE: 2024-05-15
#
#
# =============================================================================
# Clear the work space
rm(list = ls())

# Clear the console
cat("\014")

# -----------------------------------------------------------------------------
# Load functions
# -----------------------------------------------------------------------------
source(here::here("R", "barma.R"))

# -----------------------------------------------------------------------------
# Load required packages
# -----------------------------------------------------------------------------
library(forecast)       # For time series forecasting
library(zoo)            # For working with yearmon dates
library(moments)

# -----------------------------------------------------------------------------
# Load the Sobradinho dataset
# -----------------------------------------------------------------------------

# Load the Sobradinho dataset containing monthly time series data
# The dataset has two columns:
#   - time: Month names with year (e.g., "Jan 1999")
#   - y: Numeric values representing the time series measurements
sobradinho_data <- read.csv(here::here("data", "sobradinho.csv"))
sobradinho_df <- data.frame(sobradinho_data)

# Generate monthly dates using yearmon format (YYYY-MM)
# This avoids locale issues with month names in different languages
# Creates 301 monthly periods starting from January 1999
sobradinho_df$time <- seq(from = as.yearmon("Jan 1999"),
                         by = 1/12, length.out = 301)

# Create the time series
sobradinho_ts <- ts(sobradinho_df$y,
                   start = c(1999, 01),
                   frequency = 12)

# -----------------------------------------------------------------------------
# Descriptive statistics
# -----------------------------------------------------------------------------
# Some plots of useful volume from Sobradinho reservoir
plot(stl(sobradinho_ts, s.window = 12))

# Seasonality plots
monthplot(sobradinho_ts, ylab = " ")

# Time series plot
plot(sobradinho_ts,
     ylim = c(0,1),
     ylab = "", xlab = "Year",
     xaxt = "n"
)
axis(1, seq(1999, 2024, 2))

# Autocorrelation function
acf(sobradinho_ts)

# Partial autocorrelation function
pacf(sobradinho_ts)

# Descriptive statistics: Useful volume of the sobradinho reservoir
descriptive_df <- data.frame(
  min = min(sobradinho_ts),
  max = max(sobradinho_ts),
  median = median(sobradinho_ts),
  mean = mean(sobradinho_ts),
  sd = sd(sobradinho_ts),
  skewness = moments::skewness(sobradinho_ts),
  excess_kurtosis = moments::kurtosis(sobradinho_ts
  )
)
descriptive_df <- round(descriptive_df, 4)
descriptive_df

# -----------------------------------------------------------------------------
# Split the data into training and test sets
# -----------------------------------------------------------------------------
stop_train_date <- as.yearmon("2023-01")

y_train_loc <- sobradinho_df$time <= stop_train_date
y_test_loc <- sobradinho_df$time > stop_train_date
y_train_df <- sobradinho_df[y_train_loc, ]
y_test_df <- sobradinho_df[y_test_loc, ]

# Convert the training and test sets to time series
y_train_ts <- ts(y_train_df$y, end = stop_train_date, frequency = 12)
y_test_ts <- ts(y_test_df$y, start = stop_train_date + 1/12, frequency = 12)


# -----------------------------------------------------------------------------
# Create regressors for the training and test sets
# -----------------------------------------------------------------------------
len_y_train <- length(y_train_ts)
len_y_test <- length(y_test_ts)
vec_train <- 1:length(y_train_ts)
vec_test <- (max(vec_train) + 1):(max(vec_train) + len_y_test)

# Create a dummy variable for the first dry period
start_dry_period_A = as.yearmon("1999-01")
end_dry_period_A = as.yearmon("2003-12")

# Create a dummy variable for the second dry period
start_dry_period_B = as.yearmon("2012-10")
end_dry_period_B  = as.yearmon("2020-01")

dry_period_time_A <- (sobradinho_df$time >= start_dry_period_A &
                        sobradinho_df$time <= end_dry_period_A)

dry_period_time_B <- (sobradinho_df$time >= start_dry_period_B &
                        sobradinho_df$time <= end_dry_period_B)

dry_period_time <- dry_period_time_A | dry_period_time_B

dummy_dry <- as.numeric(dry_period_time)
dummy_dry_train <- dummy_dry[vec_train]
dummy_dry_test <- dummy_dry[vec_test]
y_train_df$dummy_dry <- dummy_dry_train

# Create harmonic terms and their interactions with the dry period dummy
hs_train <- sin(2 * pi * vec_train / 12)
hc_train <- cos(2 * pi * vec_train / 12)
hs_test <- sin(2 * pi * vec_test / 12)
hc_test <- cos(2 * pi * vec_test / 12)


# Create monthly dummy variables and their interactions with the dry period
d04_train <- as.numeric(ifelse(cycle(y_train_ts) == 4, 1, 0))
d04_test <- as.numeric(ifelse(cycle(y_test_ts) == 4, 1, 0))


# Construct the regressor matrices for the training and test sets
X_train <- cbind(hs = hs_train,
                 hc = hc_train,
                 dummy_dry = dummy_dry_train,
                 d04 = d04_train)

X_test <- cbind(hs = hs_test,
                hc = hc_test,
                dummy_dry = dummy_dry_test,
                d04 = d04_test)

# -----------------------------------------------------------------------------
# Fit the \beta ARMA model with AR(1, 12, 13) MA(1, 25) and the regressors
# -----------------------------------------------------------------------------
ar_vec <- c(1, 12, 13)
ma_vec <- c(1, 25)
fit_BARMA <- barma(y_train_ts,
                   ar = ar_vec,
                   ma = ma_vec,
                   link = "cloglog",
                   X = X_train,
                   X_hat = X_test)

# Print the fitted model
fit_BARMA$model

# ---------------------------------------------------------------------------
# Fit SARIMA with the training time series
# ---------------------------------------------------------------------------
fit_SARIMA <- forecast::Arima(y = y_train_ts,
                              order = c(3, 0, 1),
                              seasonal = c(2, 0, 0))

# package ‘forecast’ version 8.22.0
# > forecast::auto.arima(y = y_train_ts,
# +                                    d = 0, D = 0,
# +                                    allowdrift = FALSE)
# Series: y_train_ts
# ARIMA(3,0,1)(2,0,0)[12] with non-zero mean
#
# Coefficients:
#          ar1      ar2     ar3      ma1    sar1    sar2    mean
#       2.3044  -1.8427  0.5264  -0.7975  0.2702  0.2551  0.5055
# s.e.  0.1351   0.2046  0.0834   0.1319  0.0629  0.0669  0.1048
#
# sigma^2 = 0.003055:  log likelihood = 426.88
# AIC=-837.76   AICc=-837.24   BIC=-808.42

# ---------------------------------------------------------------------------
# Fit SARIMAX with the training time series
# ---------------------------------------------------------------------------
fit_SARIMAX <- forecast::Arima(y = y_train_ts,
                               order = c(3, 0, 0),
                               seasonal = c(2, 0, 0),
                               xreg =  cbind(
                                 dummy_dry = dummy_dry_train,
                                 d04 = d04_train))

# package ‘forecast’ version 8.22.0
# > forecast::auto.arima(y = y_train_ts,
# +                                     d = 0, D = 0,
# +                                     allowdrift = FALSE,
# +                                     xreg =  cbind(
# +                                       dummy_dry = dummy_dry_train,
# +                                       d04 = d04_train))
# Series: y_train_ts
# Regression with ARIMA(3,0,0)(2,0,0)[12] errors
#
# Coefficients:
#          ar1      ar2     ar3    sar1    sar2  intercept  dummy_dry     d04
#       1.5003  -0.6737  0.0846  0.3173  0.2527     0.5160    -0.0414  0.0102
# s.e.  0.0647   0.1038  0.0597  0.0576  0.0667     0.0784     0.0272  0.0126
#
# sigma^2 = 0.003066:  log likelihood = 426.83
# AIC=-835.66   AICc=-835.01   BIC=-802.66

# ---------------------------------------------------------------------------
# Fit ETS with the training time series
# ---------------------------------------------------------------------------
fit_ETS <- forecast::ets(y_train_ts, model = c("AAA"), damped = TRUE)

# package ‘forecast’ version 8.22.0
# > forecast::ets(y_train_ts)
# ETS(A,Ad,A)
#
# Call:
#   forecast::ets(y = y_train_ts)
#
# Smoothing parameters:
#   alpha = 0.9802
#   beta  = 0.0062
#   gamma = 0.0196
#   phi   = 0.8191
#
# Initial states:
#   l = 0.5297
#   b = -0.0486
#   s = -0.18 -0.2013 -0.1682 -0.0828 -0.0098 0.0225
#         0.1021 0.1358 0.2009 0.1592 0.0615 -0.0402
#
# sigma:  0.0613
#
#      AIC      AICc       BIC
# 42.19790  44.73123 108.19358

# =============================================================================
# Forecast
# =============================================================================

# \beta ARMA forecast values
forecast_BARMA <- fit_BARMA$forecast

# SARIMA forecast values
forecast_SARIMA <- forecast(fit_SARIMA, h = len_y_test)$mean

# SARIMAx forecast values
forecast_SARIMAX <- forecast(fit_SARIMAX,
                             xreg =  cbind(dummy_dry = dummy_dry_test,
                                           d04 = d04_test),
                             h = len_y_test)$mean

# ETS forecast values
forecast_ETS <- forecast(fit_ETS, h = len_y_test)$mean

# =============================================================================
# Evaluation metrics
# =============================================================================

# ---------------------------------------------------------------------------
# Calculate evaluation metrics for the test set
# ---------------------------------------------------------------------------
evaluation_metrics <- function(actual_values, forecasts) {
  # Calculate Mean Absolute Error (MAE)
  mae <- mean(abs(forecasts - actual_values))

  # Calculate Root Mean Squared Error (RMSE)
  rmse <- sqrt(mean((forecasts - actual_values)^2))

  return(data.frame(MAE = mae,
                    RMSE = rmse))
}


# ---------------------------------------------------------------------------
# Initialize result matrices
MAE_BARMA <- RMSE_BARMA <- matrix(nrow = 1, ncol = len_y_test)
MAE_SARIMA <- RMSE_SARIMA <- matrix(nrow = 1, ncol = len_y_test)
MAE_SARIMAX <- RMSE_SARIMAX <- matrix(nrow = 1, ncol = len_y_test)
MAE_ETS <- RMSE_ETS <- matrix(nrow = 1, ncol = len_y_test)


# ---------------------------------------------------------------------------
# Evaluate models on the test set
# ---------------------------------------------------------------------------
for (i in 1:len_y_test) {
  # Calculate the actual values from the test set
  actual_values <- y_test_ts[1:i]

  # Get the current forecasts for each model
  forecast_current_BARMA <- forecast_BARMA[1:i]
  forecast_current_SARIMA <- forecast_SARIMA[1:i]
  forecast_current_SARIMAX <- forecast_SARIMAX[1:i]
  forecast_current_ETS <- forecast_ETS[1:i]

  # Calculate evaluation metrics for each model
  metrics_BARMA <- evaluation_metrics(actual_values, forecast_current_BARMA)
  metrics_SARIMA <- evaluation_metrics(actual_values, forecast_current_SARIMA)
  metrics_SARIMAX <- evaluation_metrics(actual_values, forecast_current_SARIMAX)
  metrics_ETS <- evaluation_metrics(actual_values, forecast_current_ETS)

  # Store the evaluation metrics in the result matrices
  MAE_BARMA[i] <- metrics_BARMA$MAE
  MAE_SARIMA[i] <- metrics_SARIMA$MAE
  MAE_SARIMAX[i] <- metrics_SARIMAX$MAE
  MAE_ETS[i] <- metrics_ETS$MAE

  RMSE_SARIMA[i] <- metrics_SARIMA$RMSE
  RMSE_SARIMAX[i] <- metrics_SARIMAX$RMSE
  RMSE_BARMA[i] <- metrics_BARMA$RMSE
  RMSE_ETS[i] <- metrics_ETS$RMSE
}


# ---------------------------------------------------------------------------
# Consolidate the results
MAE_results <- rbind(
  MAE_BARMA,
  MAE_SARIMA,
  MAE_SARIMAX,
  MAE_ETS
)

RMSE_results <- rbind(
  RMSE_BARMA,
  RMSE_SARIMA,
  RMSE_SARIMAX,
  RMSE_ETS
)

# Round the results to two decimal places
MAE_results <- round(MAE_results * 100, 2)
RMSE_results <- round(RMSE_results * 100, 2)

rownames(MAE_results) <- c("MAE_BARMA",
                           "MAE_SARIMA",
                           "MAE_SARIMAX",
                           "MAE_ETS"
)

rownames(RMSE_results) <- c("RMSE_BARMA",
                            "RMSE_SARIMA",
                            "RMSE_SARIMAX",
                            "RMSE_ETS"
)

colnames(MAE_results) <- colnames(RMSE_results) <- 1:12

# Print the results
print(MAE_results)
print(RMSE_results)


# =============================================================================
# Plot comparision
# =============================================================================
# data
# X11()
plot(x = 1:12,
     y = y_test_ts,
     xlab = "Forecast horizon", ylab = "",
     xaxt = "n",
     ylim = c(0.35, 1),
     type = "p",
     pch = 19)
axis(1, at = 1:12, labels = 1:12)

# forecasts
lines(x = 1:12, y = forecast_BARMA)
lines(x = 1:12, y = forecast_SARIMA, lty = "longdash")
lines(x = 1:12, y = forecast_SARIMAX, lty = "twodash")
lines(x = 1:12, y = forecast_ETS, lty = "dotted")

legend(
  legend = c("data",
             expression(paste(beta, "ARMA")),
             "SARIMA",
             "SARIMAX",
             "ETS"),

  col = c(1, 1, 1, 1, 1,1 ),
  lty = c(NA, "solid","longdash", "twodash", "dotted"),
  pch = c(19, NA, NA, NA, NA, NA),
  "bottomleft",
  bty = "n"
)


# -----------------------------------------------------------------------------
# Fit the \beta ARMA model with AR(1, 23) MA(12, 2, 23), regressors using
# loglog link function
# -----------------------------------------------------------------------------
fit_BARMA_logit <- barma(y_train_ts,
                         ar = ar_vec,
                         ma = ma_vec,
                         link = "loglog",
                         X = X_train,
                         X_hat = X_test)


# -----------------------------------------------------------------------------
# Fit the \beta ARMA model with AR(1, 23) MA(12, 2, 23), regressors using
# cloglog link function
# -----------------------------------------------------------------------------
fit_BARMA_cloglog <- barma(y_train_ts,
                           ar = ar_vec,
                           ma = ma_vec,
                           link = "cloglog",
                           X = X_train,
                           X_hat = X_test)


# -----------------------------------------------------------------------------
# Print the fitted model
fit_BARMA_logit$model
fit_BARMA$model
fit_BARMA_cloglog$model
