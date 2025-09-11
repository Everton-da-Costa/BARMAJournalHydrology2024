#' Data preparation script for Itaparica reservoir data
#'
#' This script processes raw data for the Itaparica reservoir and creates
#' both data.frame and time series objects for use in the package.

# Load necessary libraries
library(here)
library(zoo)
library(usethis)

# Read and validate raw data
itaparica_raw_df <- tryCatch({
  read.csv(here::here("data-raw", "itaparica.csv"))
}, error = function(e) {
  stop("Error reading itaparica.csv: ", e$message)
})

# Create sequence of months
first_date <- as.yearmon("jan 1999", format = "%b %Y")
n_rows <- nrow(itaparica_raw_df)
all_dates <- seq(first_date, by = 1/12, length.out = n_rows)

# Create the data frame with complete time sequence
itaparica_df <- data.frame(
  time = all_dates,
  y = itaparica_raw_df$y
)

# Handle values larger than 1 (with warning)
values_larger_than_one <- itaparica_df$y >= 1
if (any(values_larger_than_one)) {
  warning(sprintf(
    "Found %d values > 1. These will be set to 0.9999",
    sum(values_larger_than_one)
  ))
  print(itaparica_df[values_larger_than_one, ])
  itaparica_df$y[values_larger_than_one] <- 0.9999
}

# Create time series object
itaparica_ts <- ts(
  itaparica_df$y,
  start = c(1999, 1),
  frequency = 12
)

# Add metadata attributes
attr(itaparica_df, "source") <- "Itaparica Reservoir Data"
attr(itaparica_df, "description") <- "Monthly useful volume of Itaparica reservoir"
attr(itaparica_df, "date_processed") <- Sys.Date()

# Save data objects
usethis::use_data(itaparica_df, overwrite = TRUE)
usethis::use_data(itaparica_ts, overwrite = TRUE)

# Document data cleaning steps
cat("Data processing completed:\n",
    "- Loaded raw data from itaparica.csv\n",
    "- Created time sequence from Jan 1999\n",
    "- Handled values > 1\n",
    "- Created data.frame and time series objects\n",
    "- Added metadata\n",
    "- Saved processed data\n",
    file = "data-raw/processing_log.txt",
    append = TRUE)
