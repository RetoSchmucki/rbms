#' n_components
#' Number of groups of years linked through the sites with an observation. Two years are linked when a site has
#' an observation in both, directly or through other years and sites. With one group, the level of all years
#' can be related in the collated index model; with more, the site and year effects of the groups cannot be
#' related and the collated index of the bootstrap sample can run away. Years without any observation are not counted.
#' @param site vector Site id of each observation (e.g. SINDEX > 0).
#' @param year vector Year of each observation, same length as \code{site}.
#' @return integer The number of groups of linked years, \code{Inf} if there is no observation.
#' @keywords internal bootstrap
#' @seealso \link{boot_sample}, \link{collated_index}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#'

n_components <- function(site, year){

  if(length(site) == 0) return(Inf)
  s <- as.integer(factor(site))
  y <- as.integer(factor(year))
  comp <- seq_len(max(y))

  repeat {
    # a site takes the lowest label of its years, a year the lowest label of its sites
    new <- as.integer(tapply(tapply(comp[y], s, min)[s], y, min))
    if(all(new == comp)) break
    comp <- new
  }

  length(unique(comp))
}


#' draw_strata
#' Draw one stratified bootstrap sample of site indices, sampling with replacement within each stratum the number
#' of sites in that stratum. Uses \code{sample.int()} on positions so that a stratum with a single site is not
#' expanded to \code{1:x} as \code{sample(x)} would do.
#' @param strata_ind list Integer vectors of site indices, one per stratum.
#' @return integer vector The site indices of the bootstrap sample.
#' @keywords internal bootstrap
#' @seealso \link{boot_sample}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#'

draw_strata <- function(strata_ind){

  unlist(lapply(strata_ind, function(s) s[sample.int(length(s), length(s), replace = TRUE)]), use.names = FALSE)
}
