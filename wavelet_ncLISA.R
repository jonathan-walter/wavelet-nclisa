## wavelet version of ncLISA from Ottar Bjornstad's 'ncf' package by Jonathan Walter

## wavelet.ncLISA uses the real part of the cross-wavelet transform instead of Pearson correlation
## as the synchrony statistic and is otherwise identical to 'lisa.nc' from 'ncf.' 

## inputs:
##  x: vector of length n representing the x coordinates (or latitude; see latlon)
##  y: vector of length n representing the y coordinates (or longitude; see latlon)
##  z: a matrix of dimensions n x p representing p(>=20) observation at each location
##  neigh: neighborhood size
##  timescales: a range of timescales, given as c(min,max). Defaults to all timescales.
##  resamp: number of resamples under the NULL to generate p-values
##  latlon: if TRUE, coordinates are latitude and longitude
##  quiet: if TRUE, the counter is suppressed during execution

##  Depends on the 'wsyn' package by Reuman and colleagues, available on CRAN

##  for wavelet analysis, z should be detrended and variance standardized, which can be accomplished using
##  wsyn::cleandat. 

library(wsyn)

wavelet.ncLISA<-function(x, y, z, neigh, timescales=c(0,Inf), resamp = 1000, latlon = FALSE, 
          quiet = FALSE) 
{
  
  if (is.null(dim(z))) {
    stop("\n z is univariate. Use ncf::lisa()")
  }
  if (any(!is.finite(unlist(z)))) {
    stop("Missing values in z; deal with these before using wavelet.ncLISA")
  }
  
  n <- dim(z)[1]
  p <- dim(z)[2]
  if (latlon) {
    dmat <- gcdist(x, y)
  }
  else {
    dmat <- sqrt(outer(x, x, "-")^2 + outer(y, y, "-")^2)
  }
  z <- as.matrix(z) + 0
  zx <- synmat(z, times=1:ncol(z), method="ReXWT", tsrange=timescales)
  dkl <- ifelse(dmat > 0 & dmat < neigh, 1, NA)
  nlok <- apply(dkl, 2, sum, na.rm = TRUE)
  dmean <- apply(dmat * dkl, 2, mean, na.rm = TRUE)
  moran <- apply(zx * dkl, 2, mean, na.rm = TRUE)
  p <- NULL
  if (resamp > 0) {
    perm <- matrix(NA, nrow = resamp, ncol = n)
    for (i in 1:resamp) {
      whn = pretty(c(1, resamp), n = 10)
      if (!quiet & any(i == whn)) {
        cat(i, " of ", resamp, "\\r")
        flush.console()
      }
      trekk <- sample(1:n)
      zx2 <- zx[trekk, trekk]
      perm[i, ] <- apply(zx2 * dkl, 2, mean, na.rm = TRUE)
    }
    p <- (apply(moran < t(perm), 1, sum))/(resamp + 1)
    p <- apply(cbind(p, 1 - p), 1, min) + 1/(resamp + 1)
  }
  res <- list(correlation = moran, p = p, n = nlok, dmean = dmean, 
              coord = list(x = x, y = y), call = deparse(match.call()))
  class(res) <- "lisa.nc"
  return(res)
}