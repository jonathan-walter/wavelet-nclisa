## Local wavelet phase synchrony

## Code by Jonathan Walter, borrowing from Ottar Bjornstad's 'ncf' package, with dependencies
## on Dan Reuman's 'wsyn' package and the 'circular' package

## localPhaseSync measures the phase synchrony in local neighborhoods as the mean resultant vector

## this version averages time-localized phase synchronies between a focal cell and all its neighbors

## inputs:
##  x: vector of length n representing the x coordinates (or latitude; see latlon)
##  y: vector of length n representing the y coordinates (or longitude; see latlon)
##  z: a matrix of dimensions n x p representing p(>=20) observation at each location
##  neigh: neighborhood size (a distance)
##  timescales: a range of timescales, given as c(min,max). Defaults to all timescales.
##  resamp: number of resamples under the NULL to generate p-values
##  latlon: if TRUE, coordinates are latitude and longitude
##  quiet: if TRUE, the counter is suppressed during execution

##  Depends on the 'wsyn' package (under development) available on GitHub.
##  wsyn replaces the older "Reumannplatz" package. There isn't anything wrong
##  with the Reumannplatz functions we used, but wsyn is better tested and
##  will receive continued support. Reumannplatz will not.
##  Install it using the following code:
##  library(devtools)
##  install_github('reumandc/wsyn')

##  for wavelet analysis, z should be detrended and variance standardized, which can be accomplished using:
##  z<-Reumannplatz::CleanData(z, normalize=F, detrend=T, rescale=T)$cleandat
##  If normalize=T, a box cox procedure is used to find an optimal transformation to normalize the input
##  timeseries, but this is generally not necessary unless the wavelet quantities are to be significance tested,
##  which is not the case here.
##  OR use:
##  z<-wsyn::cleandat(z, 1:p, clev=2)$cdat
##  See help(cleandat) for the meaning of different values of "clev".

##  Note that I found something glitchy that I have not yet had the chance to fix: "CleanData" overwrites and
##  deletes any object called 'x' in the global environment, so don't call one of the variables 'x'


localPhaseSynch<-function(x, y, z, neigh, timescales=c(0,Inf), latlon = FALSE, 
                          quiet = FALSE) 
{
  
  library(wsyn)
  library(circular)
  
  if (is.null(dim(z))) {
    stop("\n z is univariate.")
  }
  if (any(!is.finite(unlist(z)))) {
    stop("Missing values in z; deal with these before using localPhaseSynch")
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
  
  ## produce a matrix of locations x timesteps containing mean phase angles over relevant timescales
  phi<-matrix(NA, n, p)
  for(nn in 1:n){
    wn<-wsyn::wt(z[nn,], times=1:p)
    phi[nn,]<-Arg(rowMeans(wn$values[, wn$timescales>min(timescales) & wn$timescales<max(timescales)], 
                      na.rm=T))
  }
  
  phasesyn<-matrix(0, n, n)
  for(ii in 2:n){
    for(jj in 1:(ii-1)){
      phasesyn[ii,jj]<-mean(suppressWarnings(apply(cbind(phi[ii,],phi[jj,]), 1, rho.circular, na.rm=TRUE)),na.rm=T)
    }
  }
  phasesyn<-phasesyn+t(phasesyn)
  diag(phasesyn)<-NA
  #phimat<-as.circular(phimat, type="angles", units="radians", modulo="asis", zero=0, rotation="counter")
  dkl <- ifelse(dmat >= 0 & dmat < neigh, 1, NA)
  nlok <- apply(dkl, 2, sum, na.rm = TRUE)
  dmean <- apply(dmat * dkl, 2, mean, na.rm = TRUE)
  phaseSynch <- apply(phasesyn * dkl, 2, mean, na.rm = TRUE)
  
  res <- list(phaseSynch=phaseSynch, n = nlok, dmean = dmean, 
              coord = list(x = x, y = y), call = deparse(match.call()))
  return(res)
}