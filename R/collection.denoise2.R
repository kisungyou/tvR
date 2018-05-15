
# 1. TV-L2 Primal Dual Denoising ------------------------------------------
#    under 'RcppCollection_Image.cpp'
#' @keywords internal
#' @noRd
denoise2.TVL2.PrimalDual <- function(image, lambda, niter, iscimg){
  if (dim(image)[4]==1){ ## gray scale case
    imnormal = rcpp_01normalize(image[,,1,1])
    output   = image_tvl2_primaldual(imnormal, lambda, niter)
    if (iscimg){
      oldw <- getOption("warn")
      options(warn = -1)
      output = as.cimg(output)
      options(warn = oldw)
    }
  } else {
    output = array(0,c(dim(image)[1],dim(image)[2],dim(image)[4]))
    for (i in 1:dim(image)[4]){
      imnormal    = rcpp_01normalize(image[,,1,i])
      output[,,i] = image_tvl2_primaldual(imnormal, lambda, niter)
    }
    if (iscimg){
      oldw <- getOption("warn")
      options(warn = -1)
      output = as.cimg(output)
      options(warn = oldw)
    }
  }
  return(output)
}



# 2. TV-L2 Finite Difference ----------------------------------------------
#    under 'RcppCollection_Image.cpp'
#' @keywords internal
#' @noRd
denoise2.TVL2.FiniteDifference <- function(image, lambda, niter, iscimg){
  if (dim(image)[4]==1){ ## gray scale case
    imnormal = rcpp_01normalize(image[,,1,1])
    output   = image_tvl2_FD(imnormal, lambda, niter)
    if (iscimg){
      oldw <- getOption("warn")
      options(warn = -1)
      output = as.cimg(output)
      options(warn = oldw)
    }
  } else {
    output = array(0,c(dim(image)[1],dim(image)[2],dim(image)[4]))
    for (i in 1:dim(image)[4]){
      imnormal    = rcpp_01normalize(image[,,1,i])
      output[,,i] = image_tvl2_FD(imnormal, lambda, niter)
    }
    if (iscimg){
      oldw <- getOption("warn")
      options(warn = -1)
      output = as.cimg(output)
      options(warn = oldw)
    }
  }
  return(output)
}





# 3. TV-L1 PrimalDual -----------------------------------------------------
#    under 'RcppCollection_Image.cpp'
#' @keywords internal
#' @noRd
denoise2.TVL1.PrimalDual <- function(image, lambda, niter, iscimg){
  if (dim(image)[4]==1){ ## gray scale case
    imnormal = rcpp_01normalize(image[,,1,1])
    output   = image_tvl1_primaldual(imnormal, lambda, niter)
    if (iscimg){
      oldw <- getOption("warn")
      options(warn = -1)
      output = as.cimg(output)
      options(warn = oldw)
    }
  } else {
    output = array(0,c(dim(image)[1],dim(image)[2],dim(image)[4]))
    for (i in 1:dim(image)[4]){
      imnormal    = rcpp_01normalize(image[,,1,i])
      output[,,i] = image_tvl1_primaldual(imnormal, lambda, niter)
    }
    if (iscimg){
      oldw <- getOption("warn")
      options(warn = -1)
      output = as.cimg(output)
      options(warn = oldw)
    }
  }
  return(output)
}
