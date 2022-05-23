#' Sample Generator for ZINAR(p)
#'
#' This function generates a realization of a ZINAR(p) process.
#'
#' @param alpha The p-dimensional vector (in which p is the process order) of alpha values, the probabilities of an element remaining in the process.All alpha elements must be in [0,1] and their sum must be smaller than 1.
#' @param lambda The Poisson rate parameter. Must be greater than zero.
#' @param n The length of the simulated chain.
#' @param pii The probability of a structural zero (i.e., ignoring the Poisson distribution) under ZIP innovation sequences. Defaults to 0, following a standard Poisson.
#'
#' @return Returns a numeric vector representing a realization of an INAR/ZINAR(p) process.
#'
#' @references  Garay, Aldo M., Francyelle L. Medina, Celso RB Cabral, and Tsung-I. Lin. "Bayesian analysis of the p-order integer-valued AR process with zero-inflated Poisson innovations." Journal of Statistical Computation and Simulation 90, no. 11 (2020): 1943-1964.
#'
#'Garay, Aldo M. ; Medina, Francyelle L. ; Jales, Isaac C. ; Bertail, Patrice. First-order integer valued AR processes with zero-inflated innovations. Cyclostationarity: Theory and Methods, Springer Verlag - 2021, v. 1, p. 19-40.
#'
#' @export

simul_zinarp <- function(n, alpha, lambda, pii = 0){
  if(n%%1!=0|n<=0) stop("n must me a positive integer")
  if(sum(alpha<=0|alpha>=1) > 0) stop("all alpha values must be in (0,1)")
  if(sum(alpha)>=1) stop("sum of alpha values must be less than 1")
  if(lambda<=0) stop("lambda must be greater than zero")
  if(pii>=1|pii<0) stop("pii must be in [0,1)")
  if(length(pii)>1) stop("pii must have length 1")
  if(length(lambda)>1) stop("lambda must have length 1")
  if(length(n)>1) stop("n must have length 1")
  if(pii==0){
    result <- genINARp(alpha = alpha,lambda = lambda, n = n)
  }else{
    result <- genZINARp(alpha = alpha, pii = pii,lambda = lambda,n = n)
  }
  return(result)
}
