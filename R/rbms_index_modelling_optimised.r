##==========================================
## Deprecated optimised versions of flight_curve() and get_nm()
##
## They gave the same results as flight_curve() and get_nm(); the option
## to fit the years in parallel is now in flight_curve(parallel = TRUE).
##
##      Date:   17.10.2026
##
##==========================================


#' flight_curve_optimised
#' Deprecated, use \link{flight_curve}, which gives the same flight curves and fits the years in parallel with
#' \code{parallel = TRUE}. The arguments are passed to \link{flight_curve}; it will be removed in a future version.
#' @param ... Arguments passed to \link{flight_curve}.
#' @return See \link{flight_curve}.
#' @keywords internal
#' @seealso \link{flight_curve}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @export flight_curve_optimised
#'

flight_curve_optimised <- function(...){

  .Deprecated("flight_curve", package = "rbms", msg = "flight_curve_optimised() is deprecated, use flight_curve() (same results, parallel = TRUE to fit the years in parallel).")
  flight_curve(...)
}


#' get_nm_optimised
#' Deprecated, use \link{get_nm}. The arguments are passed to \link{get_nm}; it will be removed in a future version.
#' @param ... Arguments passed to \link{get_nm}.
#' @return See \link{get_nm}.
#' @keywords internal
#' @seealso \link{get_nm}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @export get_nm_optimised
#'

get_nm_optimised <- function(...){

  .Deprecated("get_nm", package = "rbms", msg = "get_nm_optimised() is deprecated, use get_nm().")
  get_nm(...)
}
