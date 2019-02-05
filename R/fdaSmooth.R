#' @title Smooth spectra
#'
#' @description Compute a smoothed version of the spectra using splines.
#'
#' @param X_rough A numeric matrix of size \code{N*P}, where \code{N} is
#' the number of observations (spectra) and \code{P} is the number of points
#' (wavelengths) at which spectra are measured.
#'
#' @param wgts A numeric vector of weights for smoothing spectra.
#'
#' @param wvls A numeric vector giving the wavelenghts at which spectra were measured.
#'
#' @return A matrix the same size of \code{X_rough} containing the smoothed spectra.

#' @references
#'
#' @seealso
#'
#' @examples
#'
#' @export
#'
fdaSmooth <- function(X_rough, wgts=NULL, wvls=350:2500){

  NN <- nrow(X_rough)
  PP <- ncol(X_rough)
  ww <- if(is.null(wgts)){ rep(1,PP) }else{ wgts }
  X_smooth <- matrix(0,NN,PP)

  for(i in 1:NN){
    XX <- X_rough[i,]
    aux <- smooth.spline(x = wvls, w = ww, y = X_rough[i,], cv=F)
    X_smooth[i,] <- predict(aux)$y
  }

  return(X_smooth)
}
