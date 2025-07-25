#' A simulation with three groups under the Proportional Hazards hypothesis
#'
#' This dataframe represents a simulation of a study under the Proportional Hazards
#' hypothesis. All three survival curves follows an exponential distribution with
#' different parameters.
#'
#' @format ## `data_under_PH`
#' A data frame with 600 rows and 3 columns:
#' \describe{
#'   \item{time}{Time of events in months with truncation at 60 months}
#'   \item{status}{Indicator of censorship. 1 denotes an event, 0 denotes a censor}
#'   \item{arm}{Integer from 0 to 2. Indicates the group the patient belongs to}
#' }
"data_under_PH"



#' A simulation with three groups without the Proportional Hazards hypothesis
#'
#' This dataframe represents a simulation of a study where the three survival curves
#' cross each other.
#'
#' @format ## `data_not_PH`
#' A data frame with 600 rows and 3 columns:
#' \describe{
#'   \item{time}{Time of events in months with truncation at 60 months}
#'   \item{status}{Indicator of censorship. 1 denotes an event, 0 denotes a censor}
#'   \item{arm}{Integer from 0 to 2. Indicates the group the patient belongs to}
#' }
"data_not_PH"

