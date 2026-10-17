##==========================================
## Multi-year flight curve: moved to develop/flight_curve_multi/ in the rbms repository
## (https://github.com/RetoSchmucki/rbms) for further development.
##
##      Date:   17.10.2026
##
##==========================================


multi_moved <- function(f){
  stop(paste0(f, "() has been moved out of the rbms package for further development. Use flight_curve(), or source the ",
              "code from develop/flight_curve_multi/fit_bam_multi.r in https://github.com/RetoSchmucki/rbms"), call. = FALSE)
}


#' flight_curve_multi
#' Deprecated, moved to develop/flight_curve_multi/fit_bam_multi.r in the rbms GitHub repository for further development;
#' use \link{flight_curve}. Calling it returns an error with this information.
#' @param ... not used.
#' @keywords internal
#' @seealso \link{flight_curve}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @export flight_curve_multi
#'

flight_curve_multi <- function(...) multi_moved("flight_curve_multi")


#' fit_bam_multi
#' Deprecated, moved to develop/flight_curve_multi/fit_bam_multi.r in the rbms GitHub repository for further development;
#' use \link{flight_curve}. Calling it returns an error with this information.
#' @param ... not used.
#' @keywords internal
#' @seealso \link{flight_curve}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @export fit_bam_multi
#'

fit_bam_multi <- function(...) multi_moved("fit_bam_multi")


#' plot_flight_curve_multi
#' Deprecated, moved to develop/flight_curve_multi/fit_bam_multi.r in the rbms GitHub repository for further development;
#' use \link{plot.pheno_curve}. Calling it returns an error with this information.
#' @param ... not used.
#' @keywords internal
#' @seealso \link{plot.pheno_curve}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @export plot_flight_curve_multi
#'

plot_flight_curve_multi <- function(...) multi_moved("plot_flight_curve_multi")
