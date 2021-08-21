#' Total Variation Denoising for Signal
#'
#' Given a 1-dimensional signal \code{f}, it solves an optimization of the form,
#' \deqn{u^* = argmin_u E(u,f)+\lambda V(u)}
#' where \eqn{E(u,f)} is fidelity term and \eqn{V(u)} is total variation regularization term.
#' The naming convention of a parameter \code{method} is \code{<problem type>} + \code{<name of algorithm>}.
#' For more details, see the section below.
#'
#' @section Algorithms for TV-L2 problem:
#' The cost function for TV-L2 problem is
#' \deqn{min_u \frac{1}{2} |u-f|_2^2 + \lambda |\nabla u|}
#' where for a given 1-dimensional vector, \eqn{|\nabla u| = \sum |u_{i+1}-u_{i}|}.
#' Algorithms (in conjunction with model type) for this problems are \describe{
#'   \item{\code{"TVL2.IC"}}{Iterative Clipping algorithm.}
#'   \item{\code{"TVL2.MM"}}{Majorization-Minorization algorithm.}
#' }
#' The codes are translated from MATLAB scripts by Ivan Selesnick.
#'
#'
#' @param signal vector of noisy signal.
#' @param lambda regularization parameter (positive real number).
#' @param niter  total number of iterations.
#' @param method indicating problem and algorithm combination.
#'
#' @return a vector of same length as input \code{signal.}
#'
#' @examples
#' \donttest{
#' ## generate a stepped signal
#' x = rep(sample(1:5,10,replace=TRUE), each=50)
#'
#' ## add some additive white noise
#' xnoised = x + rnorm(length(x), sd=0.25)
#'
#' ## apply denoising process
#' xproc1 = denoise1(xnoised, method = "TVL2.IC")
#' xproc2 = denoise1(xnoised, method = "TVL2.MM")
#'
#' ## plot noisy and denoised signals
#' plot(xnoised, pch=19, cex=0.1, main="Noisy signal")
#' lines(xproc1, col="blue", lwd=2)
#' lines(xproc2, col="red", lwd=2)
#' legend("bottomleft",legend=c("Noisy","TVL2.IC","TVL2.MM"),
#' col=c("black","blue","red"),#' lty = c("solid", "solid", "solid"),
#' lwd = c(0, 2, 2), pch = c(19, NA, NA),
#' pt.cex = c(1, NA, NA), inset = 0.05)
#' }
#'
#' @references
#' \insertRef{rudin_nonlinear_1992}{tvR}
#'
#' \insertRef{selesnick_convex_2015}{tvR}
#'
#' @export
denoise1 <- function(signal, lambda=1.0, niter=100, method=c("TVL2.IC","TVL2.MM")){
  ## Check Data, Lambda, niter
  if (!check_data_signal(signal)){
    stop("* denoise1 : input 'signal' should be a vector with no NA or Inf values allowed.")
  }
  signal = as.vector(signal)
  if (!check_lambda(lambda)){
    stop("* denoise1 : 'lambda' should be a positive real number.")
  }
  if (!check_niter(niter)){
    stop("* denoise1 : 'niter' should be a positive integer larger than 1.")
  }
  ## Method Argument
  if (missing(method)){
    method = "TVL2.IC"
  } else {
    method = match.arg(method)
  }

  ## Main Computation
  output = switch(method,
                  TVL2.IC = denoise1.TVL2.IC(signal, lambda, niter),
                  TVL2.MM = denoise1.TVL2.MM(signal, lambda, niter)
                  )
  return(output)
}
