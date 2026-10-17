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


#' year_components
#' Group of linked years for each year, see \link{n_components}.
#' @param site vector Site id of each observation (e.g. count > 0).
#' @param year vector Year of each observation, same length as \code{site}.
#' @return integer vector Group label of each year, named by year (sorted as \code{levels(factor(year))}).
#' @keywords internal bootstrap
#' @seealso \link{n_components}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#'

year_components <- function(site, year){

  s <- as.integer(factor(site))
  y <- factor(year)
  comp <- seq_len(nlevels(y))
  y <- as.integer(y)

  repeat {
    new <- as.integer(tapply(tapply(comp[y], s, min)[s], y, min))
    if(all(new == comp)) break
    comp <- new
  }

  stats::setNames(as.integer(factor(comp)), levels(factor(year)))
}


#' sum_by
#' Sum of x by integer group g in 1:n, groups without value are 0.
#' @keywords internal
#' @noRd

sum_by <- function(x, g, n){

  r <- numeric(n)
  s <- rowsum(x, g, reorder = TRUE)
  r[as.integer(rownames(s))] <- s[, 1]
  r
}


#' site_year_level
#' Expected count per unit of flight curve (level) for each site and year, from the Poisson GLM
#' COUNT ~ factor(SITE_ID) + factor(M_YEAR) + offset(log(NM)), fitted without building the model matrix (no limit on
#' the number of sites). The counts are first summed per site-year, with the flight curve covered (sum of NM).
#' With method = "newton", the site effects are concentrated out of the likelihood (closed form given the year effects)
#' and the concentrated likelihood, a multinomial per site, is maximised by Newton steps on the year effects only. One
#' year per group of linked years (see \link{n_components}) is the reference and years without count have a level of 0.
#' If Newton fails, iterative proportional fitting is used. With method = "ipf", iterative proportional fitting matches
#' the fitted total of each site and each year to the observed total. Both give the maximum likelihood estimates of glm().
#' @param site vector Site id of each count.
#' @param year vector Year of each count.
#' @param count numeric Observed counts.
#' @param nm numeric Proportion of the flight curve (NM) of each count, > 0.
#' @param method character Fitting method, "newton" (default) or "ipf".
#' @param tol numeric Tolerance for convergence, default 1e-8.
#' @param maxit integer Maximum number of iterations, default 100 for Newton and 10000 for IPF.
#' @return numeric vector The level of the site-year of each count, the fitted count being \code{level * nm}. The attribute
#' "effects" is a list with the site (site) and year (year) multiplicative effects, named by site and year, level = site x year
#' (the scale between them is arbitrary, their product is not). With Newton,
#' the attribute "year_effect" is a data.table with, per year, the group of linked years (GROUP), the log year effect
#' relative to the first year of the group (LOG_EFFECT, -Inf for years without count), its standard error (SE) scaled
#' by the quasi-Poisson dispersion and conditional on the site totals, and the dispersion (DISPERSION).
#' @keywords internal site index
#' @seealso \link{impute_count}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#'

site_year_level <- function(site, year, count, nm, method = c("newton", "ipf"), tol = 1e-8, maxit = NULL){

  method <- match.arg(method)

  ## site x year table: total count (C) and flight curve covered (N)
  fs <- factor(site)
  rs <- as.integer(fs)
  fy <- factor(year)
  ry <- as.integer(fy)
  n_s <- max(rs)
  n_y <- nlevels(fy)
  sy <- data.table::data.table(S_IDX = rs, Y_IDX = ry, C_SY = count, N_SY = nm)[, .(C_SY = sum(C_SY), N_SY = sum(N_SY)), by = .(S_IDX, Y_IDX)]
  C_s <- sum_by(sy$C_SY, sy$S_IDX, n_s)
  C_y <- sum_by(sy$C_SY, sy$Y_IDX, n_y)

  b <- NULL
  year_effect <- NULL

  if(method == "newton"){
    fit <- try(newton_year(sy, C_s, C_y, n_s, n_y, tol = tol, maxit = if(is.null(maxit)) 100 else maxit), silent = TRUE)
    if(inherits(fit, "try-error")){
      warning(paste("Newton on the year effects failed (", trimws(as.character(fit)), "), iterative proportional fitting used"), call. = FALSE)
    } else {
      b <- fit$b
      year_effect <- data.table::data.table(M_YEAR = levels(fy), fit$year_effect)
    }
  }

  if(is.null(b)){
    b <- rep(1, n_y)
    level <- rep(1, nrow(sy))
    for(it in seq_len(if(is.null(maxit)) 10000 else maxit)){
      a <- C_s / sum_by(b[sy$Y_IDX] * sy$N_SY, sy$S_IDX, n_s)
      a[!is.finite(a)] <- 0
      b <- C_y / sum_by(a[sy$S_IDX] * sy$N_SY, sy$Y_IDX, n_y)
      new <- a[sy$S_IDX] * b[sy$Y_IDX]
      conv <- max(abs(new - level) / pmax(new, 1e-10)) < tol
      level <- new
      if(conv) break
    }
    if(!conv) warning(paste("site_year_level did not converge after", it, "iterations"), call. = FALSE)
  }

  ## site effects in closed form given the year effects, level of each site-year
  a <- C_s / sum_by(b[sy$Y_IDX] * sy$N_SY, sy$S_IDX, n_s)
  a[!is.finite(a)] <- 0
  level_sy <- numeric(n_s * n_y)
  level_sy[(sy$S_IDX - 1) * n_y + sy$Y_IDX] <- a[sy$S_IDX] * b[sy$Y_IDX]
  level <- level_sy[(rs - 1) * n_y + ry]
  data.table::setattr(level, "effects", list(site = stats::setNames(a, levels(fs)), year = stats::setNames(b, levels(fy))))

  if(!is.null(year_effect)){
    ## quasi-Poisson dispersion from the counts, scales the standard errors
    mu <- level * nm
    df <- length(count) - (sum(C_s > 0) + sum(C_y > 0) - data.table::uniqueN(year_effect[is.finite(LOG_EFFECT), GROUP]))
    phi <- if(df > 0) sum(((count - mu)^2 / mu)[mu > 0]) / df else NA_real_
    year_effect[, SE := SE * sqrt(phi)][, DISPERSION := phi]
    data.table::setattr(level, "year_effect", year_effect[])
  }

  level
}


#' newton_year
#' Newton steps on the concentrated likelihood of the year effects, see \link{site_year_level}.
#' @param sy data.table Site-year table with site (S_IDX), year (Y_IDX), count (C_SY) and flight curve covered (N_SY).
#' @param C_s,C_y numeric Total count per site and per year.
#' @param n_s,n_y integer Number of sites and years.
#' @param tol numeric Tolerance on the gradient, relative to the largest year total.
#' @param maxit integer Maximum number of Newton steps.
#' @return list with the year effects b (exp of the log effects, 0 for years without count) and a data.table with
#' GROUP, LOG_EFFECT and SE (not scaled by the dispersion) per year.
#' @keywords internal site index
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#'

newton_year <- function(sy, C_s, C_y, n_s, n_y, tol = 1e-8, maxit = 100){

  active <- C_y > 0
  d <- sy[active[Y_IDX], ]
  if(nrow(d) == 0) stop("no count")

  ## one reference year (log effect 0) per group of years linked through sites with counts
  group <- rep(NA_integer_, n_y)
  pos <- d[C_SY > 0, ]
  group[active] <- year_components(pos$S_IDX, pos$Y_IDX)[as.character(which(active))]
  ref <- active & !duplicated(group) & !is.na(group)
  free <- active & !ref

  beta <- rep(0, n_y)
  beta[!active] <- -Inf

  loglik <- function(beta){
    D <- sum_by(exp(beta[d$Y_IDX]) * d$N_SY, d$S_IDX, n_s)
    sum(d$C_SY * beta[d$Y_IDX]) - sum((C_s * log(D))[C_s > 0])
  }
  neg_hessian <- function(beta){
    e <- exp(beta[d$Y_IDX]) * d$N_SY
    D <- sum_by(e, d$S_IDX, n_s)
    p <- e / D[d$S_IDX]                     # share of the site total expected in the year
    fit_y <- sum_by(C_s[d$S_IDX] * p, d$Y_IDX, n_y)
    P <- matrix(0, n_s, n_y)
    P[cbind(d$S_IDX, d$Y_IDX)] <- sqrt(C_s[d$S_IDX]) * p
    list(g = C_y - fit_y, H = diag(fit_y, n_y) - crossprod(P))
  }

  conv <- !any(free)
  it <- 0
  while(!conv && it < maxit){
    it <- it + 1
    gh <- neg_hessian(beta)
    if(max(abs(gh$g[free])) < tol * max(C_y)){ conv <- TRUE; break }
    delta <- solve(gh$H[free, free, drop = FALSE], gh$g[free])
    ## step halving to make sure the likelihood increases
    l0 <- loglik(beta)
    t <- 1
    repeat {
      nb <- beta
      nb[free] <- beta[free] + t * delta
      if(is.finite(loglik(nb)) && loglik(nb) >= l0 - 1e-12 * abs(l0)) break
      t <- t / 2
      if(t < 1e-10) stop("no increase of the likelihood")
    }
    beta <- nb
  }
  if(!conv) stop(paste("no convergence after", maxit, "steps"))

  se <- rep(NA_real_, n_y)
  se[ref] <- 0
  if(any(free)) se[free] <- sqrt(diag(solve(neg_hessian(beta)$H[free, free, drop = FALSE])))

  list(b = exp(beta), year_effect = data.table::data.table(GROUP = group, LOG_EFFECT = beta, SE = se))
}


#' index_sample_data
#' Data of a species for a collated index: site-years with NA in \code{sindex_value} removed (warning), the sites of the
#' original data or of the bootstrap sample \code{bootID} (each draw of a site with its own uID), the groups of years linked
#' through sites with observations in the original data, and the number of sites per year. Used by \link{collated_index}
#' and \link{collated_index_mixed}.
#' @param data,s_sp,sindex_value,bootID,boot_ind see \link{collated_index}.
#' @return list with the species data (data), bootID, the sites of the sample (a: uID, SITE_ID), the group of each year
#' (grp), the data of the sample (data_boot) and the number of sites per year (sum_data); or, for a bootstrap sample that
#' failed in \link{boot_sample} (NA row), the output of the collated index with NA (failed).
#' @keywords internal
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#'

index_sample_data <- function(data, s_sp, sindex_value = "SINDEX", bootID = NULL, boot_ind = NULL){

  data <- setDT(data)
  data = data[SPECIES == s_sp, ]

  ## site-years without index (e.g. failed flight curve) are missing, not zero: removed as in site_index()
  n_na <- data[is.na(get(sindex_value)), .N]
  if(n_na > 0){
    warning(paste0(s_sp, ": ", n_na, " site-years with NA in ", sindex_value, " removed before computing the collated index"), call. = FALSE)
    data <- data[!is.na(get(sindex_value)), ]
  }
  if(nrow(data) == 0) stop(paste0("no ", sindex_value, " available for species ", s_sp), call. = FALSE)

  y = as.numeric(as.character(unique(data$M_YEAR)))

  if(is.null(boot_ind)){
    site_id_ <- unique(data[, SITE_ID])
    boot_ind_ <- matrix(NA, nrow = 1, ncol = length(site_id_))
    boot_ind_[1, ] <- seq_len(length(site_id_))
    boot_ind <- list(site_id = site_id_, boot_ind = boot_ind_)
  }

  if(is.null(bootID)){
    bootID <- 0
  }

  ## bootstrap samples drawn per species with boot_sample(s_sp = ) must match the species modelled
  if(!is.null(boot_ind$s_sp) && !identical(as.character(boot_ind$s_sp), as.character(s_sp))){
    stop(paste0("boot_ind was generated for species ", boot_ind$s_sp, ", not for ", s_sp), call. = FALSE)
  }

  ## bootstrap sample that failed in boot_sample() after max_trial trials (NA row): no index computed
  if(bootID > 0 && anyNA(boot_ind$boot_ind[bootID, ])){
    res <- data.table(BOOTi = bootID, M_YEAR = y[order(y)], NSITE = NA_integer_, NSITE_B = NA_integer_,
                      NSITE_OBS = NA_integer_, COL_INDEX = NA_real_, GROUP = NA_integer_)
    setkey(res, M_YEAR)
    return(list(failed = list(col_index = res, site_id = NA_character_, year_effect = NULL)))
  }

  if(bootID == 0){
    a <- data.table(uID = as.character(paste0("s_",seq_along(boot_ind$boot_ind[bootID+1, ]))),
     SITE_ID = as.character(boot_ind$site_id[seq_along(boot_ind$boot_ind[bootID+1, ])]))
  } else {
    a <- data.table(uID = as.character(paste0("s_",seq_along(boot_ind$boot_ind[bootID, ]))),
     SITE_ID = as.character(boot_ind$site_id[boot_ind$boot_ind[bootID, ]]))
  }

  s_y <- data.table(expand.grid(uID = a$uID, M_YEAR = y[order(y)]))
  s_y[ , uID:=as.character(uID)]
  s_y[ , M_YEAR:=as.numeric(as.character(M_YEAR))]
  data[ , SITE_ID:=as.character(SITE_ID)]
  data[ , M_YEAR:=as.numeric(as.character(M_YEAR))]

  ## groups of years linked through sites with observations in the original data (see n_components),
  ## years without observation have no group
  obs_0 <- data[get(sindex_value) > 0, ]
  year_group <- if(nrow(obs_0) > 0) year_components(obs_0$SITE_ID, obs_0$M_YEAR) else integer(0)
  grp <- data.table(M_YEAR = as.numeric(names(year_group)), GROUP = as.integer(year_group))

  setkey(a, uID)
  setkey(s_y, uID)
  s_y <- merge(s_y, a, all.x=TRUE)

  setkey(s_y, SITE_ID, M_YEAR); setkey(data, SITE_ID, M_YEAR)

  if(sindex_value != 'SINDEX'){
    data_boot <- merge(s_y, data[ , .(SITE_ID, M_YEAR, SINDEX, TOTAL_NM)][, (sindex_value) := data[, get(sindex_value)]])
    } else {
    data_boot <- merge(s_y, data[ , .(SITE_ID, M_YEAR, SINDEX, TOTAL_NM)])
    }

  sum_data <- data.table( BOOTi = bootID,
                          M_YEAR = y[order(y)],
                          NSITE = data_boot[, uniqueN(uID), by = M_YEAR][order(M_YEAR), V1],
                          NSITE_B = data_boot[, uniqueN(SITE_ID), by = M_YEAR][order(M_YEAR), V1],
                          NSITE_OBS = data_boot[, sum(get(sindex_value)>0), by = M_YEAR][order(M_YEAR), V1])

  setkey(sum_data, M_YEAR)

  list(data = data, bootID = bootID, a = a, grp = grp, data_boot = data_boot, sum_data = sum_data, failed = NULL)
}


#' site_length_map
#' Length of each site (transect) from a column of the data, checked to be positive, not NA and constant within site.
#' @param data data.table Data of a species with SITE_ID.
#' @param site_length string Name of the column with the length of the site (km).
#' @return data.table with SITE_ID and SITE_LENGTH.
#' @keywords internal
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#'

site_length_map <- function(data, site_length){

  if(!site_length %in% names(data)) stop(paste0("column ", site_length, " not found in data"), call. = FALSE)
  lm_ <- data.table(SITE_ID = as.character(data$SITE_ID), SITE_LENGTH = as.numeric(data[[site_length]]))[
    , .(N_L = uniqueN(SITE_LENGTH), SITE_LENGTH = SITE_LENGTH[1]), by = SITE_ID]
  if(any(lm_$N_L > 1)) stop(paste0(site_length, " must be constant within site"), call. = FALSE)
  if(any(is.na(lm_$SITE_LENGTH) | lm_$SITE_LENGTH <= 0)) stop(paste0(site_length, " must be positive and not NA (length in km)"), call. = FALSE)
  lm_[, .(SITE_ID, SITE_LENGTH)]
}
