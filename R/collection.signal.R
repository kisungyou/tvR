# 1. TV-L2 Iterative Clipping by Selesnick --------------------------------
#    under 'RcppCollection_Signal.cpp'
#' @keywords internal
#' @noRd
denoise1.TVL2.IC <- function(signal, lambda, niter){
  output = signal_tvl2_IC(signal, lambda, niter);
  return(as.vector(output));
}


# 2. TV-L2 Majorization-Minorization by Selesnick -------------------------
#    under 'RcppCollection_Signal.cpp'
#' @keywords internal
#' @noRd
denoise1.TVL2.MM <- function(signal, lambda, niter){
  output = signal_tvl2_MM(signal, lambda, niter);
  return(as.vector(output));
}
