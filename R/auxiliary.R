#' @keywords internal
#' @noRd
check_datamatrix <- function(A){
  cond1 = (is.matrix(A)||(inherits(A, "Matrix")))
  cond2 = (!any(is.infinite(A)))
  cond3 = (!any(is.na(A)))
  if (cond1&&cond2&&cond3){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @keywords internal
#' @noRd
check_lambda <- function(lambda){
  if (length((lambda))==1){
    if (lambda > 0){
      if (!is.na(lambda)){
        if (!is.infinite(lambda)){
          return(TRUE)
        }
      }
    }
  }
  return(FALSE)
}

#' @keywords internal
#' @noRd
check_niter <- function(niter){
  if (length((niter))==1){
    if (niter > 1){
      if ((!is.na(niter))&&(!is.infinite(niter))){
        if (abs(niter-round(niter))<sqrt(.Machine$double.eps)){
          return(TRUE)
        }
      }
    }
  }
  return(FALSE)
}

#' @keywords internal
#' @noRd
check_data_signal <- function(signal){
  cond1 = ((is.vector(signal)) ||(is.matrix(signal)&&(min(nrow(signal),ncol(signal))==1)))
  cond2 = all(!is.na(signal))
  cond3 = all(!is.infinite(signal))
  if (cond1&&cond2&&cond3){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


#  For Non-CIMG type case.
#' @keywords internal
#' @noRd
check_data_image <- function(image){
  ldX = length(dim(image))
  cond1 = (ldX==2)
  cond2 = ((ldX==3)&&(dim(image)[3]==3))
  if (cond1||cond2){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
