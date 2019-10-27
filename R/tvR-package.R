#' tvR : Total Variation Regularization
#'
#' \bold{tvR} provides tools for denoising noisy signal and images via
#' Total Variation Regularization. Reducing the total variation of the given signal is known to remove spurious detail while preserving
#' essential structural details. For now, we provide tools for \emph{denoising} only on 1-dimensional signals or 2-dimensional images,
#' where the latter be represented as 2d or 3d array.
#'
#' @docType package
#' @name tvR-package
#' @import Matrix
#' @import Rdpack
#' @importFrom utils packageVersion
#' @importFrom Rcpp evalCpp
#' @useDynLib tvR
NULL


