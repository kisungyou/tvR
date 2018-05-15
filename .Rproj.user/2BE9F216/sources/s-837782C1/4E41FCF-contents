#' Total Variation Denoising for Image
#'
#' Given an image \code{f}, it solves an optimization of the form,
#' \deqn{u^* = argmin_u E(u,f)+\lambda V(u)}
#' where \eqn{E(u,f)} is fidelity term and \eqn{V(u)} is total variation regularization term.
#' The naming convention of a parameter \code{method} is \code{<problem type>} + \code{<name of algorithm>}.
#' For more details, see the section below.
#'
#' @section Data format:
#' An input \code{image} can have one of three following types: 2-dimensional matrix representaing \emph{grayscale} image, 3-dimensional array
#' for \emph{color} image whose 3rd dimension has size of 3, or \pkg{cimg} object represented as 4-d array.
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
#' @param image standard 2d or 3d array, or \code{cimg} object.
#' @param lambda regularization parameter (positive real number).
#' @param niter  total number of iterations.
#' @param method indicating problem and algorithm combination.
#'
#' @return denoised image of same data type as input \code{image}.
#'
#' @examples
#' ## Load data from 'imager' library and rescale it
#' library(imager)
#' x <- imager::imresize(imager::boats, scale=0.5)
#'
#' ## Add white noise
#' xnoised <- x + imager::imnoise(dim=dim(x), sd=0.1); plot(xnoised)
#'
#' ## apply denoising models
#' xproc1 <- denoise2(xnoised, method="TVL2.FiniteDifference")
#' xproc2 <- denoise2(xnoised, method="TVL1.PrimalDual")
#'
#' ## compare
#' par(mfrow=c(2,2))
#' plot(x, main="original")
#' plot(xnoised, main="noised")
#' plot(xproc1, main="TVL2.FiniteDifference")
#' plot(xproc2, main="TVL1.PrimalDual")
#'
#' @references
#' \insertRef{rudin_nonlinear_1992}{tvR}
#'
#' \insertRef{chambolle_first-order_2011}{tvR}
#'
#' @export
denoise2 <- function(image, lambda=1.0, niter=100, method=c("TVL1.PrimalDual","TVL2.PrimalDual","TVL2.FiniteDifference")){
  ## Check Data, Lambda, niter
  ##    For image data as cimg, it's fine.
  if (!check_lambda(lambda)){
    stop("* denoise2 : 'lambda' should be a positive real number.")
  }
  if (!check_niter(niter)){
    stop("* denoise2 : 'niter' should be a positive integer larger than 1.")
  }
  checkcimg = is.cimg(image)
  if (!checkcimg){
    if (!check_data_image(image)){
      stop("* denoise2 : an 'image' array should be 2d or 3d with 3rd dimension having size of 3.")
    }
    oldw  = getOption("warn")
    options(warn = -1)
    image = as.cimg(image)
    options(warn=oldw)
  }

  ## Method Argument
  if (missing(method)){
    method = "TVL2.PrimalDual"
  } else {
    method = match.arg(method)
  }

  ## Main Computation Part
  output = switch(method,
                  TVL2.PrimalDual       = denoise2.TVL2.PrimalDual(image, 1.0/lambda, niter, checkcimg),
                  TVL2.FiniteDifference = denoise2.TVL2.FiniteDifference(image, 1.0/lambda, niter, checkcimg),
                  TVL1.PrimalDual       = denoise2.TVL1.PrimalDual(image, 1.0/lambda, niter, checkcimg)
                  )

  ## return output
  return(output)
}

## NOTE that for TVL1 and TVL2, sum(sqrt(Ix^2 + Iy^2)) + lambda*||I - g||_2^2 is default used.
#' ROF : fixed point iteration (http://www.math.ucla.edu/~lvese/285j.1.05s/ROFScheme.pdf)
