#' Simulated data from the trial and the observational study
#'
#' The simulated dataset for combined data, including the trial data and the
#' observational study data.
#'
#' @format ## `dat`
#' A data frame with 2500 rows and 6 columns, including the trial data with the
#' sample size as 500 and the observational study data with the sample size as
#' 1000. The columns contain the following variables.
#' \describe{
#'   \item{x1}{The continuous baseline covariates}
#'   \item{x2}{The binary baseline covariate}
#'   \item{a}{The binary treatment assignment}
#'   \item{t}{The observed event time}
#'   \item{delta}{The event indicator, where `delta = 1` indicates the occurrence
#'   of the event and `delta = 0` indicates censoring}
#'   \item{s}{The data source indicator, where `s = 1` indicates the trial
#'   sample and `s = 0` indicates the observational study sample.}
#' }
"dat"
