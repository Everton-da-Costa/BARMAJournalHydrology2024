# =============================================================================
# Application Tres Marias data set
# =============================================================================
# PAPER: Dynamic beta modeling of seasonal hydro-environmental time series
# Test inferences and link function selection in dynamic beta modeling of
# seasonal hydro-environmental time series with temporary abnormal regimes
#
# DESCRIPTION: The Tres Marias data set represents the useful volume of Tres
# Marias reservoirs. This code is used to analyze the data as described in the
# associated paper.
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
# Load the Tres Marias dataset
# -----------------------------------------------------------------------------

# Load the Tres Marias dataset containing monthly time series data
# The dataset has two columns:
#   - time: Month names with year (e.g., "Jan 1999")
#   - y: Numeric values representing the time series measurements
tres_marias_data <- read.csv(here::here("data", "tres_marias.csv"))
tres_marias_df <- data.frame(tres_marias_data)

# Generate monthly dates using yearmon format (YYYY-MM)
# This avoids locale issues with month names in different languages
# Creates 301 monthly periods starting from January 1999
tres_marias_df$time <- seq(from = as.yearmon("Jan 1999"),
                         by = 1/12, length.out = 301)

# Create the time series
tres_marias_ts <- ts(tres_marias_df$y,
                   start = c(1999, 01),
                   frequency = 12)
# -----------------------------------------------------------------------------
# Descriptive statistics
# -----------------------------------------------------------------------------
# Some plots of useful volume from Tres Marias reservoir
plot(stl(tres_marias_ts, s.window = 12))

# Seasonality plots
monthplot(tres_marias_ts, ylab = " ")

# Time series plot
plot(tres_marias_ts,
     ylim = c(0,1),
     ylab = "", xlab = "Year",
     xaxt = "n"
)
axis(1, seq(1999, 2024, 2))

# Autocorrelation function
acf(tres_marias_ts)

# Partial autocorrelation function
pacf(tres_marias_ts)

# Descriptive statistics: Useful volume of the tres_marias reservoir
descriptive_df <- data.frame(
  min = min(tres_marias_ts),
  max = max(tres_marias_ts),
  median = median(tres_marias_ts),
  mean = mean(tres_marias_ts),
  sd = sd(tres_marias_ts),
  skewness = moments::skewness(tres_marias_ts),
  excess_kurtosis = moments::kurtosis(tres_marias_ts
  )
)
descriptive_df <- round(descriptive_df, 4)
descriptive_df

# -----------------------------------------------------------------------------
# Split the data into training and test sets
# -----------------------------------------------------------------------------
stop_train_date <- as.yearmon("2023-01")

y_train_loc <- tres_marias_df$time <= stop_train_date
y_test_loc <- tres_marias_df$time > stop_train_date
y_train_df <- tres_marias_df[y_train_loc, ]
y_test_df <- tres_marias_df[y_test_loc, ]

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
start_dry_period_B = as.yearmon("2013-03")
end_dry_period_B  = as.yearmon("2018-10")

dry_period_time_A <- (tres_marias_df$time >= start_dry_period_A &
                        tres_marias_df$time <= end_dry_period_A)

dry_period_time_B <- (tres_marias_df$time >= start_dry_period_B &
                        tres_marias_df$time <= end_dry_period_B)

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
d03_train <- as.numeric(ifelse(cycle(y_train_ts) == 3, 1, 0))
d03_test <- as.numeric(ifelse(cycle(y_test_ts) == 3, 1, 0))


# Construct the regressor matrices for the training and test sets
X_train <- cbind(hs = hs_train,
                 hc = hc_train,
                 dummy_dry = dummy_dry_train,
                 d03 = d03_train)

X_test <- cbind(hs = hs_test,
                hc = hc_test,
                dummy_dry = dummy_dry_test,
                d03 = d03_test)

# -----------------------------------------------------------------------------
# Fit the \beta ARMA model with AR(1, 23) MA(12, 2, 23) and the regressors
# -----------------------------------------------------------------------------
ar_vec <- c(1, 23)
ma_vec <- c(1, 2, 23)
fit_BARMA <- barma(y_train_ts,
                   ar = ar_vec,
                   ma = ma_vec,
                   link = "logit",
                   X = X_train,
                   X_hat = X_test)

# Print the fitted model
fit_BARMA$model

# ---------------------------------------------------------------------------
# Fit SARIMAX with the training time series
# ---------------------------------------------------------------------------
fit_SARIMA <- forecast::Arima(y = y_train_ts,
                              order = c(2, 0, 0),
                              seasonal = c(2, 0, 0))

# package ‘forecast’ version 8.22.0
# > forecast::auto.arima(y = y_train_ts,
# +                                    d = 0, D = 0,
# +                                    allowdrift = FALSE)
# Series: y_train_ts
# ARIMA(2,0,0)(2,0,0)[12] with non-zero mean
#
# Coefficients:
#          ar1      ar2    sar1    sar2    mean
#       1.3884  -0.4510  0.2741  0.3886  0.5789
# s.e.  0.0574   0.0583  0.0553  0.0627  0.1271
#
# sigma^2 = 0.002771:  log likelihood = 438.52
# AIC=-865.05   AICc=-864.75   BIC=-843.05
# ---------------------------------------------------------------------------
# Fit SARIMAX with the training time series
# ---------------------------------------------------------------------------
fit_SARIMAX <- forecast::Arima(y = y_train_ts,
                               order = c(3, 0, 1),
                               seasonal = c(2, 0, 0),
                               xreg =  cbind(
                                 dummy_dry = dummy_dry_train,
                                 d03 = d03_train))

# package ‘forecast’ version 8.22.0
# > forecast::auto.arima(y = y_train_ts,
# +                                     d = 0, D = 0,
# +                                     allowdrift = FALSE,
# +                                     xreg =  cbind(
# +                                       dummy_dry = dummy_dry_train,
# +                                       d03 = d03_train))
# Series: y_train_ts
# Regression with ARIMA(3,0,1)(2,0,0)[12] errors
#
# Coefficients:
#   ar1     ar2      ar3     ma1    sar1    sar2  intercept  dummy_dry
# 0.8077  0.3713  -0.3015  0.5556  0.2332  0.4127     0.6112    -0.0875
# s.e.  0.3721  0.5034   0.1694  0.3753  0.0559  0.0627     0.0953     0.0268
# d03
# 0.0302
# s.e.  0.0152
#
# sigma^2 = 0.00267:  log likelihood = 446.09
# AIC=-872.17   AICc=-871.38   BIC=-835.51

# ---------------------------------------------------------------------------
# Fit ETS with the training time series
# ---------------------------------------------------------------------------
fit_ETS <- forecast::ets(y_train_ts, model = "ANA")

# package ‘forecast’ version 8.22.0
# > forecast::ets(y_train_ts)
# ETS(A,N,A)
#
# Call:
#   forecast::ets(y = y_train_ts)
#
# Smoothing parameters:
#   alpha = 0.895
#   gamma = 0.1031
#
# Initial states:
#   l = 0.5442
#   s = -0.1068 -0.1095 -0.1046 -0.0739 -8e-04 0.0506
#         0.104 0.098 0.1406 0.0826 -1e-04 -0.0801
#
# sigma:  0.0577
#
# AIC      AICc       BIC
# 4.350044  6.108286 59.346444

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
                                           d03 = d03_test),
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
fit_BARMA_loglog <- barma(y_train_ts,
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
fit_BARMA$model
fit_BARMA_loglog$model
fit_BARMA_cloglog$model

