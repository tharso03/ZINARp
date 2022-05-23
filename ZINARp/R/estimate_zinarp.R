#' Parameter estimation for ZINARp models
#'
#' This function uses MCMC algorithms (Metropolis-Hastings and Gibbs Sampler) to generate a chain of INAR/ZINAR(p) parameter estimators.
#'
#' @param x A vector containing a discrete non-negative time series dataset.
#' @param p The order of the INAR/ZINAR process.
#' @param iter The number of iterations to be considered. Defaults to 5000.
#' @param thin Lag for posterior sample. Defaults to 2.
#' @param burn Burn-in for posterior sample. Defaults to 0.1. Must be in (0,1).
#' @param innovation Distribution to be used for the innovation : "Poisson" or "ZIP". Defaults to Poisson.
#' @examples
#' test <- simul_zinarp(alpha = 0.1, lambda = 1, n = 100)
#' e.test <- estimate_zinarp(x = test, p = 1, iter = 800, innovation= "Poisson")
#' alpha_hat <- mean(e.test$alpha)
#' lambda_hat <- mean(e.test$lambda)
#'
#' data(slesions)
#' e.slesions <- estimate_zinarp(slesions$y, p = 1, iter = 800, innovation = 'ZIP')
#' alpha_hat_slesions <- mean(e.slesions$alpha)
#' lambda_hat_slesions <- mean(e.slesions$lambda)
#' rho_hat_slesions <- mean(e.slesions$rho)
#' @return Returns a list containing a posteriori samples for the specified model parameters.
#'
#'
#' @references  Garay, Aldo M., Francyelle L. Medina, Celso RB Cabral, and Tsung-I. Lin. "Bayesian analysis of the p-order integer-valued AR process with zero-inflated Poisson innovations." Journal of Statistical Computation and Simulation 90, no. 11 (2020): 1943-1964.
#'
#' Garay, Aldo M., Francyelle L. Medina, Isaac Jales CS, and Patrice Bertail. "First-Order Integer Valued AR Processes with Zero-Inflated Innovations." In Workshop on Nonstationary Systems and Their Applications, pp. 19-40. Springer, Cham, 2021.

#' @export


estimate_zinarp <- function(x,p,iter=5000,thin=2,burn=0.1,innovation="Poisson")
{
    if(length(x)<length(p)) stop("length of x must be greather than length of p")
    if(sum(x<0)>0) stop("all x values must be non-negative")
    if(sum(x%%1!=0)>0) stop("all x values must be integers")
    if(iter%%1!=0|iter<1) stop("number of iterations must be a positive integer")
    if(burn<0 | burn >=1) stop("burn-in must be in [0,1)")
    if(thin%%1!=0|thin<1) stop("thin must be a positive integer")
    if(thin>iter*burn) stop("thin must be smaller than iter*burn")
    if(!(innovation %in% c("Poisson","ZIP"))) stop("innovation must be set to 'Poisson' or 'ZIP'")
    x <- as.vector(x)
    n <- length(x)

    if(innovation=='Poisson'){
      results <- GibbsINARp(x,p,iter,thin,burn)
    }else{
      results <- GibbsZINARp(x,p,iter,thin,burn)
    }

    return(results)
}

