#' Monthly Useful Volume of the Itaparica Reservoir
#'
#' @description
#' A time series dataset containing the monthly useful volume, expressed as a
#' proportion (ranging from 0 to 1), of the Itaparica hydroelectric power
#' plant reservoir. The reservoir is part of the São Francisco River Basin (SFRB)
#' in Brazil. The data covers a period marked by seasonal fluctuations and
#' prolonged droughts.
#'
#' @format A data frame with 2 variables:
#' \describe{
#'   \item{time}{Month and year of observation (yearmon format)}
#'   \item{y}{Useful volume as a proportion (between 0 and 1)}
#' }
#' @source Brazilian National Electric System Operator (ONS)
"itaparica_df"

#' Itaparica Reservoir Monthly Useful Volume Time Series
#'
#' Monthly useful volume data from the Itaparica reservoir as a time series object,
#' starting from January 1999.
#'
#' @format A time series object (`ts`) with 301 monthly observations from
#' January 1999 to January 2024.
#'
#' @source Brazilian National Electric System Operator (ONS)
"itaparica_ts"

