#' Toy data set with butterfly count for x species across y sites
#' @format data.table object with butterfly count in long format
#' \describe{
#'   \item{SITE_ID}{number used as site id, the can be and integer or a string}
#'   \item{DATE}{date when the observation/count was recorded, format YEAR-MM-DD}
#'   \item{SPECIES}{number used as species id, this can be an integer or a string (e.g. "Aglais io")}
#'   \item{DAY}{integer, day within the month}
#'   \item{MONTH}{integer, month within the year}
#'   \item{YEAR}{integer, 4 digit year}
#'   \item{COUNT}{integer, acctual butterfly count}
#' }

"m_count"

#' Toy data set with the date when the sites have been visited for monitoring
#' @format data.table object with visit date long format
#' \describe{
#'   \item{SITE_ID}{number used as site id, the can be and integer or a string}
#'   \item{DATE}{date when the observation/count was recorded, format YEAR-MM-DD}
#' }

"m_visit"
