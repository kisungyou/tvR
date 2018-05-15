
# 1. TV-L2 Primal Dual Denoising ------------------------------------------
#    under 'RcppCollection_Image.cpp'
#' @keywords internal
#' @noRd
denoise2.TVL2.PrimalDual <- function(image, lambda, niter, normalize){
  nsize  = dim(image)
  if (length(nsize)==2){ # gray-scale
    if (normalize==TRUE){
      imnormal = rcpp_01normalize(image)
    } else if (normalize==FALSE) {
      imnormal = image
    } else {
      stop("")
    }
    output = image_tvl2_primaldual(imnormal, lambda, niter)
  } else {
    nn = nsize[3]
    output = array(0,nsize)
    for (i in 1:nn){
      if (normalize==TRUE){
        imnormal = rcpp_01normalize(as.matrix(image[,,i]))
      } else if (normalize==FALSE){
        imnormal = as.matrix(image[,,i])
      } else {
        stop("")
      }
      output[,,i] = image_tvl2_primaldual(imnormal, lambda, niter)
    }
  }
  return(output)
}



# 2. TV-L2 Finite Difference ----------------------------------------------
#    under 'RcppCollection_Image.cpp'
#' @keywords internal
#' @noRd
denoise2.TVL2.FiniteDifference <- function(image, lambda, niter, normalize){
  nsize  = dim(image)
  if (length(nsize)==2){ # gray-scale
    if (normalize==TRUE){
      imnormal = rcpp_01normalize(image)
    } else if (normalize==FALSE) {
      imnormal = image
    } else {
      stop("")
    }
    output = image_tvl2_FD(imnormal, lambda, niter)
  } else {
    nn = nsize[3]
    output = array(0,nsize)
    for (i in 1:nn){
      if (normalize==TRUE){
        imnormal = rcpp_01normalize(as.matrix(image[,,i]))
      } else if (normalize==FALSE){
        imnormal = as.matrix(image[,,i])
      } else {
        stop("")
      }
      output[,,i] = image_tvl2_FD(imnormal, lambda, niter)
    }
  }
  return(output)
}





# 3. TV-L1 PrimalDual -----------------------------------------------------
#    under 'RcppCollection_Image.cpp'
#' @keywords internal
#' @noRd
denoise2.TVL1.PrimalDual <- function(image, lambda, niter, normalize){
  nsize  = dim(image)
  if (length(nsize)==2){ # gray-scale
    if (normalize==TRUE){
      imnormal = rcpp_01normalize(image)
    } else if (normalize==FALSE) {
      imnormal = image
    } else {
      stop("")
    }
    output = image_tvl1_primaldual(imnormal, lambda, niter)
  } else {
    nn = nsize[3]
    output = array(0,nsize)
    for (i in 1:nn){
      if (normalize==TRUE){
        imnormal = rcpp_01normalize(as.matrix(image[,,i]))
      } else if (normalize==FALSE){
        imnormal = as.matrix(image[,,i])
      } else {
        stop("")
      }
      output[,,i] = image_tvl1_primaldual(imnormal, lambda, niter)
    }
  }
  return(output)
}
