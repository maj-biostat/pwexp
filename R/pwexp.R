pwe_interval = c(1E-8, 1e6)

#' Piecewise exponential hazard function
#'
#' @param t times at which to compute hazard
#' @param bins defines bounds of constant hazrad, e.g. c(0, 5, 10)
#' @param rate defines constant hazard in each bin
#'
#' @return
#' @export
hpw <- function(t = 2, bins = c(0,  20,  35), rate = c(1/15, 1/7, 1/20)){
  if(any(t < 0)) stop("Negative values for t not supported")
  y <- sapply(seq_along(t), function(i){
    rate[findInterval(t[i], bins)]
  })
  y
}

#' Piecewise exponential cumulative hazard function
#'
#' @param t times at which to compute hazard
#' @param bins defines bounds of constant hazrad, e.g. c(0, 5, 10)
#' @param rate defines constant hazard in each bin
#'
#' @return
#' @export
Hpw <- function(t = 2, bins = c(0,  20,  35), rate = c(1/15, 1/7, 1/20)){
  if(any(t < 0)) stop("Negative values for t not supported")
  y <- sapply(seq_along(t), function(i){
    idx <- findInterval(t[i], bins)
    if(idx>1){
      HH <- sum(diff(bins)[1:(idx-1)] * rate[1:(idx-1)])
      HH <- HH + (t[i] - bins[idx])*rate[idx]
    } else {
      HH <- (t[i] - bins[idx])*rate[idx]
    }
    HH
  })
  y
}

#' Piecewise exponential survival function
#'
#' @param t times at which to compute hazard
#' @param bins defines bounds of constant hazrad, e.g. c(0, 5, 10)
#' @param rate defines constant hazard in each bin
#'
#' @return
#' @export
Spw <- function(t = 2, bins = c(0,  20,  35), rate = c(1/15, 1/7, 1/20)){
  if(any(t < 0)) stop("Negative values for t not supported")
  exp(-Hpw(t, bins, rate))
}

#' Piecewise exponential cdf
#'
#' @param t times at which to compute hazard
#' @param bins defines bounds of constant hazrad, e.g. c(0, 5, 10)
#' @param rate defines constant hazard in each bin
#'
#' @return
#' @export
ppw <- function(t = 2, bins = c(0,  20,  35), rate = c(1/15, 1/7, 1/20)){
  if(any(t < 0)) stop("Negative values for t not supported")
  1 - Spw(t, bins, rate)
}

#' Piecewise exponential density
#'
#' @param t times at which to compute hazard
#' @param bins defines bounds of constant hazrad, e.g. c(0, 5, 10)
#' @param rate defines constant hazard in each bin
#'
#' @return
#' @export
dpw <- function(t = 2, bins = c(0,  20,  35), rate = c(1/15, 1/7, 1/20)){
  if(any(t < 0)) stop("Negative values for t not supported")
  hpw(t, bins, rate) * exp(-Hpw(t, bins, rate))
}

#' Constrains value to machine precision bounds
#'
#' @param x arbitrary value
#'
#' @return
#' @keywords internal
return_finite <- function(x){
  x <- min(x, .Machine$double.xmax)
  x <- max(x, -.Machine$double.xmax)
  x
}

#' Function for which to find root
#'
#' @param t time point at which to compute surv
#' @param u random uniform
#' @param bins bounds of constant hazard
#' @param rate hazard in each bin
#'
#' @return
#' @keywords internal
rootfn <- function(t, u, bins, rate){
  lambda <- Hpw(t, bins, rate)
  surv <- exp(-lambda);
  return_finite(log(surv) - log(u))
}

#' Piecewise exponential random number generator
#'
#' @param bins bounds of constant hazard
#' @param rate hazard in each bin
#'
#' @return
#' @keywords internal
#' @importFrom stats runif
rand_draw <- function(bins, rate){
  u_i <- runif(1)
  at_limit <- rootfn(pwe_interval[2], u = u_i,
                     bins = bins,
                     rate = rate)
  if(at_limit > 0){
    return(c(pwe_interval[2], 0))
  } else {
    # trying to find the value for t at which S(t) - u = 0
    t_i <- stats::uniroot(rootfn,
                          u = u_i,
                          interval = pwe_interval,
                          check.conv = TRUE,
                          trace = 100,
                          bins = bins,
                          rate = rate)$root
    # our random sample from the surv dist
    return(t_i)
  }
}

#' Piecewise exponential random number generator API
#'
#' @param n number of draws
#' @param bins bounds of constant hazard
#' @param rate hazard in each bin
#'
#' @return
#' @export
#'
rpw <- function(n = 1, bins = c(0,  20,  35), rate = c(1/15, 1/7, 1/20)){
  replicate(n, rand_draw(bins, rate))
}



