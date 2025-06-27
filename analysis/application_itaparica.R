# =============================================================================
# Application Itaparica data set
# =============================================================================
# PAPER: Dynamic beta modeling of seasonal hydro-environmental time series
# Test inferences and link function selection in dynamic beta modeling of
# seasonal hydro-environmental time series with temporary abnormal regimes
#
# DESCRIPTION: The Itaparica data set represents the useful volume of Itaparica
# reservoirs. This code is used to analyze the data as described in the
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
# Load the Itaparica dataset
# -----------------------------------------------------------------------------

# Load the Itaparica dataset containing monthly time series data
# The dataset has two columns:
#   - time: Month names with year (e.g., "Jan 1999")
#   - y: Numeric values representing the time series measurements
itaparica_data <- read.csv(here::here("data", "itaparica.csv"))
itaparica_df <- data.frame(itaparica_data)

# Generate monthly dates using yearmon format (YYYY-MM)
# This avoids locale issues with month names in different languages
# Creates 301 monthly periods starting from January 1999
itaparica_df$time <- seq(from = as.yearmon("Jan 1999"),
                         by = 1/12, length.out = 301)

# Create the time series
itaparica_ts <- ts(itaparica_df$y,
                     start = c(1999, 01),
                     frequency = 12)

# -----------------------------------------------------------------------------
# Descriptive statistics
# -----------------------------------------------------------------------------
# Some plots of useful volume from itaparica reservoir
plot(stl(itaparica_ts, s.window = 12))

# Seasonality plots
monthplot(itaparica_ts, ylab = " ")

# Time series plot
plot(itaparica_ts,
     ylim = c(0,1),
     ylab = "", xlab = "Year",
     xaxt = "n"
)
axis(1, seq(1999, 2024, 2))

# Autocorrelation function
acf(itaparica_ts)

# Partial autocorrelation function
pacf(itaparica_ts)

# Descriptive statistics: Useful volume of the Itaparica reservoir
descriptive_df <- data.frame(
  min = min(itaparica_ts),
  max = max(itaparica_ts),
  median = median(itaparica_ts),
  mean = mean(itaparica_ts),
  sd = sd(itaparica_ts),
  skewness = moments::skewness(itaparica_ts),
  excess_kurtosis = moments::kurtosis(itaparica_ts
  )
)
descriptive_df <- round(descriptive_df, 2)
descriptive_df

# -----------------------------------------------------------------------------
# Split the data into training and test sets
# -----------------------------------------------------------------------------
stop_train_date <- as.yearmon("2023-01")

y_train_loc <- itaparica_df$time <= stop_train_date
y_test_loc <- itaparica_df$time > stop_train_date
y_train_df <- itaparica_df[y_train_loc, ]
y_test_df <- itaparica_df[y_test_loc, ]

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

# Create a dummy variable for the dry period
start_dry_period <- as.yearmon("2012-10")
end_dry_period <- as.yearmon("2020-05")

dry_period_time <- (itaparica_df$time >= start_dry_period &
                      itaparica_df$time <= end_dry_period)

dummy_dry <- as.numeric(dry_period_time)
dummy_dry_train <- dummy_dry[vec_train]
dummy_dry_test <- dummy_dry[vec_test]
y_train_df$dummy_dry <- dummy_dry_train

# Create harmonic terms and their interactions with the dry period dummy
hs_train <- sin(2 * pi * vec_train / 12)
hc_train <- cos(2 * pi * vec_train / 12)
hs1dummy_dry_train <- hs_train * (1 - dummy_dry_train)
hc1dummy_dry_train <- hc_train * (1 - dummy_dry_train)
hs_test <- sin(2 * pi * vec_test / 12)
hc_test <- cos(2 * pi * vec_test / 12)
hs1dummy_dry_test <- hs_test * (1 - dummy_dry_test)
hc1dummy_dry_test <- hc_test * (1 - dummy_dry_test)

# Create monthly dummy variables and their interactions with the dry period
d03_train <- as.numeric(ifelse(cycle(y_train_ts) == 3, 1, 0))
d10_train <- as.numeric(ifelse(cycle(y_train_ts) == 10, 1, 0))
d03_1dummy_dry_train <- d03_train * (1 - dummy_dry_train)
d10_1dummy_dry_train <- d10_train * (1 - dummy_dry_train)
d03_test <- as.numeric(ifelse(cycle(y_test_ts) == 3, 1, 0))
d10_test <- as.numeric(ifelse(cycle(y_test_ts) == 10, 1, 0))
d03_1dummy_dry_test <- d03_test * (1 - dummy_dry_test)
d10_1dummy_dry_test <- d10_test * (1 - dummy_dry_test)

# Construct the regressor matrices for the training and test sets
X_train <- cbind(hs1dummy_dry = hs1dummy_dry_train,
                 hc1dummy_dry = hc1dummy_dry_train,
                 dummy_dry = dummy_dry_train,
                 d03_1dummy_dry = d03_1dummy_dry_train,
                 d10_1dummy_dry = d10_1dummy_dry_train)

X_test <- cbind(hs1dummy_dry_hat = hs1dummy_dry_test,
                hc1dummy_dry_hat = hc1dummy_dry_test,
                dummy_dry_hat = dummy_dry_test,
                d03_1dummy_dry_hat = d03_1dummy_dry_test,
                d10_1dummy_dry_hat = d10_1dummy_dry_test)


# -----------------------------------------------------------------------------
# Fit the \beta ARMA model with AR(1, 3) MA(12, 19, 24) and the regressors
# -----------------------------------------------------------------------------
ar_vec <- c(1, 3)
ma_vec <- c(12, 19, 24)
fit_BARMA <- barma(y_train_ts,
                   ar = ar_vec,
                   ma = ma_vec,
                   link = "loglog",
                   X = X_train,
                   X_hat = X_test)

# Print the fitted model
fit_BARMA$model

# ---------------------------------------------------------------------------
# Fit SARIMAX with the training time series
# ---------------------------------------------------------------------------
fit_SARIMA <- forecast::Arima(y = y_train_ts,
                              order = c(1, 0, 0),
                              seasonal = c(1, 0, 1))

# package ‘forecast’ version 8.22.0
# > forecast::auto.arima(y = y_train_ts,
# +                      d = 0, D = 0,
# +                      allowdrift = FALSE)
# Series: y_train_ts
# ARIMA(1,0,0)(1,0,1)[12] with non-zero mean
#
# Coefficients:
#          ar1    sar1     sma1    mean
#       0.8254  0.8868  -0.6424  0.5874
# s.e.  0.0354  0.0503   0.0895  0.1088
#
# sigma^2 = 0.01624:  log likelihood = 184.4
# AIC=-358.79   AICc=-358.58   BIC=-340.46
#
# ---------------------------------------------------------------------------
# Fit SARIMAX with the training time series
# ---------------------------------------------------------------------------
fit_SARIMAX <- forecast::Arima(y = y_train_ts,
                               order = c(2, 0, 2),
                               seasonal = c(2, 0, 1),
                               xreg =  cbind(
                                 dummy_dry = dummy_dry_train,
                                 d03_1dummy_dry = d03_1dummy_dry_train,
                                 d10_1dummy_dry = d10_1dummy_dry_train))

# package ‘forecast’ version 8.22.0
# > forecast::auto.arima(y = y_train_ts, d = 0, D = 0, allowdrift = FALSE,
# +                      xreg =  cbind(dummy_dry = dummy_dry_train,
# +                                    d03_1dummy_dry = d03_1dummy_dry_train,
# +                                    d10_1dummy_dry = d10_1dummy_dry_train))
# Series: y_train_ts
# Regression with ARIMA(2,0,2)(2,0,1)[12] errors
#
# Coefficients:
#          ar1      ar2      ma1     ma2    sar1    sar2     sma1  intercept
#       1.5304  -0.6869  -0.8883  0.2469  0.6616  0.1428  -0.5380     0.6980
# s.e.  0.1067   0.0865   0.1267  0.0964  0.1747  0.0887   0.1721     0.0356
#       dummy_dry  d03_1dummy_dry  d10_1dummy_dry
#         -0.3985         -0.0604          0.0477
# s.e.     0.0443          0.0406          0.0409
#
# sigma^2 = 0.01474:  log likelihood = 202.91
# AIC=-381.82   AICc=-380.69   BIC=-337.82
#
# ---------------------------------------------------------------------------
# Fit ETS with the training time series
# ---------------------------------------------------------------------------
fit_ETS <- forecast::ets(y_train_ts, model = "ANA")

# package ‘forecast’ version 8.22.0
# > fit_ETS
# ETS(A,N,A)
#
# Call:
#   forecast::ets(y = y_train_ts)
#
# Smoothing parameters:
#   alpha = 0.7555
# gamma = 1e-04
#
# Initial states:
#   l = 0.5501
#   s = -0.1737 -0.1394 -0.0034 0.0461 0.0978 0.1185
#         0.14 0.1406 0.0807 -0.0864 -0.1228 -0.0979
#
# sigma:  0.1332
#
# AIC     AICc      BIC
# 488.1942 489.9524 543.1906
#

# =============================================================================
# Forecast
# =============================================================================

# \beta ARMA forecast values
forecast_BARMA <- fit_BARMA$forecast

# SARIMA forecast values
forecast_SARIMA <- forecast(fit_SARIMA, h = len_y_test)$mean

# SARIMAx forecast values
forecast_SARIMAX <- forecast(fit_SARIMAX,
                             xreg =  cbind(
                               dummy_dry = dummy_dry_test,
                               d03_1dummy_dry = d03_1dummy_dry_test,
                               d10_1dummy_dry = d10_1dummy_dry_test
                             ),
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
  "topleft",
  bty = "n"
)

# -----------------------------------------------------------------------------
# Fit the \beta ARMA model with AR(1, 3) MA(12, 19, 24), regressors using
# logit link function
# -----------------------------------------------------------------------------
fit_BARMA_logit <- barma(y_train_ts,
                         ar = ar_vec,
                         ma = ma_vec,
                         link = "logit",
                         X = X_train,
                         X_hat = X_test)

# -----------------------------------------------------------------------------
# Fit the \beta ARMA model with AR(1, 3) MA(12, 19, 24), regressors using
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
