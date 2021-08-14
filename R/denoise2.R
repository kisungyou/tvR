#' Total Variation Denoising for Image
#'
#' Given an image \code{f}, it solves an optimization of the form,
#' \deqn{u^* = argmin_u E(u,f)+\lambda V(u)}
#' where \eqn{E(u,f)} is fidelity term and \eqn{V(u)} is total variation regularization term.
#' The naming convention of a parameter \code{method} is \code{<problem type>} + \code{<name of algorithm>}.
#' For more details, see the section below.
#'
#' @section Data format:
#' An input \code{data} can be either (1) 2-dimensional matrix representaing \emph{grayscale} image, or (2) 3-dimensional array
#' for \emph{color} image.
#'
#' @section Algorithms for TV-L1 problem:
#' The cost function for TV-L2 problem is
#' \deqn{min_u |u-f|_1 + \lambda |\nabla u|}
#' where for a given 2-dimensional array, \eqn{|\nabla u| = \sum sqrt(u_x^2 + u_y^2)}
#' Algorithms (in conjunction with model type) for this problems are \describe{
#'   \item{\code{"TVL1.PrimalDual"}}{Primal-Dual algorithm.}
#' }
#' @section Algorithms for TV-L2 problem:
#' The cost function for TV-L2 problem is
#' \deqn{min_u |u-f|_2^2 + \lambda |\nabla u|}
#' and algorithms (in conjunction with model type) for this problems are \describe{
#'   \item{\code{"TVL2.PrimalDual"}}{Primal-Dual algorithm.}
#'   \item{\code{"TVL2.FiniteDifference"}}{Finite Difference scheme with fixed point iteration.}
#' }
#'
#'
#' @param data standard 2d or 3d array.
#' @param lambda regularization parameter (positive real number).
#' @param niter  total number of iterations.
#' @param method indicating problem and algorithm combination.
#' @param normalize a logical; \code{TRUE} to make the range in \eqn{[0,1]}, or \code{FALSE} otherwise.
#'
#' @return denoised array as same size of \code{data}.
#'
#' @examples
#' \dontrun{
#' ## Load grey-scale 'lena' data
#' data(lena128)
#'
#' ## Add white noise
#' sinfo   <- dim(lena128)   # get the size information
#' xnoised <- lena128 + array(rnorm(128*128, sd=10), sinfo)
#'
#' ## apply denoising models
#' xproc1 <- denoise2(xnoised, lambda=10, method="TVL2.FiniteDifference")
#' xproc2 <- denoise2(xnoised, lambda=10, method="TVL1.PrimalDual")
#'
#' ## compare
#' gcol = gray(0:256/256)
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' image(lena128, main="original", col=gcol)
#' image(xnoised, main="noised", col=gcol)
#' image(xproc1, main="TVL2.FiniteDifference", col=gcol)
#' image(xproc2, main="TVL1.PrimalDual", col=gcol)
#' par(opar)
#' }
#'
#' @references
#' \insertRef{rudin_nonlinear_1992}{tvR}
#'
#' \insertRef{chambolle_first-order_2011}{tvR}
#'
#' @export
denoise2 <- function(data, lambda=1.0, niter=100, method=c("TVL1.PrimalDual","TVL2.PrimalDual","TVL2.FiniteDifference"),
                     normalize=FALSE){
  ## Check Data, Lambda, niter
  ##    For image data as cimg, it's fine.
  if (!check_lambda(lambda)){
    stop("* denoise2 : 'lambda' should be a positive real number.")
  }
  if (!check_niter(niter)){
    stop("* denoise2 : 'niter' should be a positive integer larger than 1.")
  }
  if (!is.array(data)){
    stop("* denoise2 : 'data' should be either a matrix or an array.")
  }
  ndim = length(dim(data))
  if (!((ndim==2)||(ndim==3))){
    stop("* denoise2 : 'data' should be either 2d- or 3d-array.")
  }

  ## Method Argument
  if (missing(method)){
    method = "TVL2.PrimalDual"
  } else {
    method = match.arg(method)
  }

  ## Main Computation Part
  output = switch(method,
                  TVL2.PrimalDual       = denoise2.TVL2.PrimalDual(data, 1.0/lambda, niter, normalize),
                  TVL2.FiniteDifference = denoise2.TVL2.FiniteDifference(data, 1.0/lambda, niter, normalize),
                  TVL1.PrimalDual       = denoise2.TVL1.PrimalDual(data, 1.0/lambda, niter, normalize)
                  )

  ## return output
  return(output)
}

## NOTE that for TVL1 and TVL2, sum(sqrt(Ix^2 + Iy^2)) + lambda*||I - g||_2^2 is default used.
#' ROF : fixed point iteration (www.math.ucla.edu/~lvese/285j.1.05s/ROFScheme.pdf)
