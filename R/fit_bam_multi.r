##==========================================
## Multi-year flight curve: removed from the package for further development
## (the previous code is available in rbms 1.3.0).
##
##      Date:   17.10.2026
##
##==========================================


multi_moved <- function(f){
  stop(paste0(f, "() has been removed from the rbms package for further development. Use flight_curve(); ",
              "the previous code is available in rbms 1.3.0."), call. = FALSE)
}


#' flight_curve_multi
#' Deprecated, removed from the package for further development (the previous code is available in rbms 1.3.0);
#' use \link{flight_curve}. Calling it returns an error with this information.
#' @param ... not used.
#' @keywords internal
#' @seealso \link{flight_curve}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @export flight_curve_multi
#'

flight_curve_multi <- function(...) multi_moved("flight_curve_multi")


#' fit_bam_multi
#' Deprecated, removed from the package for further development (the previous code is available in rbms 1.3.0);
#' use \link{flight_curve}. Calling it returns an error with this information.
#' @param ... not used.
#' @keywords internal
#' @seealso \link{flight_curve}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @export fit_bam_multi
#'

fit_bam_multi <- function(...) multi_moved("fit_bam_multi")


#' plot_flight_curve_multi
#' Deprecated, removed from the package for further development (the previous code is available in rbms 1.3.0);
#' use \link{plot.pheno_curve}. Calling it returns an error with this information.
#' @param ... not used.
#' @keywords internal
#' @seealso \link{plot.pheno_curve}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @export plot_flight_curve_multi
#'

plot_flight_curve_multi <- function(...) multi_moved("plot_flight_curve_multi")
