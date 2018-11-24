#' @title Estimate power (probability of success separation) in two population model
#'
#' @param m1 Fraction of marker genes, 0.1 <= \code{m1} <= 10
#' @param pi1 Percentage of the novel type of cells, 5 <= \code{pi1} <= 50
#' @param foldchange Expression fold change in the novel subpopulation, 1 <= \code{foldchange} <= 32
#' @param dropout Overall dropout rate, 0.6 <= \code{dropout} <= 0.999
#' @param p Total number of genes in the simulation, 10000 <= \code{p}
#' @param n Vector of total number of cells, 50 <= min(n), max(n) <= 1000000
#' @param mu Either \code{p} dimensional mean gene expression vector or mean of mean gene expressions (hyperparameter). If supplied by a \code{p} vector, the mean and standard deviation of mean gene expressions are estimated. If supplied by a number, it is regarded as the mean of mean expression vectors. defualt is 1.
#' @param sigma Standard deviation of mean gene expressions (hyperparameter), default 5
#' @param seed Random seed, default 123
#' @param B Number of replicates to estiamte the probability of success separation, default 1000, 100 <= \code{B}
#' @param ncore Number of cores for use in the simulation using \code{\link{parallel}} package. default is 1, 1 <= \code{ncore} <= \code{parallel::detectCores()} 
#' @export
#' @return Vector of estiamted probabilities of same length as \code{n}
#' @examples \dontrun{
#' ncells(m1=1.6, pi1=5, foldchange=4, dropout=.6, p=20000, mu=1, sigma=5)
#' }
ncells <- function(m1, pi1, foldchange, dropout, p=26616, n=seq(100,1000,by=100), mu=1, sigma=5, seed=123, B=1000, ncore=1)
{
    stopifnot(5000 <= p)
    stopifnot(0.1 <= m1 && m1 <= 10)
    stopifnot(50 <= min(n) && max(n) <= 1000000)
    stopifnot(5 <= pi1 && pi1 <= 50)
    stopifnot(1 <= foldchange && foldchange <= 32)
    stopifnot(0.6 <= dropout && dropout <= 0.999)
    stopifnot(0 <= sigma)
    stopifnot(100 <= B)
    stopifnot(1 <= ncore && ncore <= parallel::detectCores())

    pfrac <- m1/100
    nfrac <- pi1/100
    logfc <- log2(foldchange)
    alpha <- dropout/(1-dropout)
    
    set.seed(seed)
    
    model <- list()
    if (length(mu) == p) {
        ## Mean and SD of prior aere estimated from the specified mu vector
        ## Permute the mean vector
        model$muvec <- sample(mu)
        model$marginalmean <- mean(mu)
        model$marginalsd <- sqrt(var(mu) + 1)
    } else if (length(mu) == 1) {
        stopifnot(length(sigma) == 1)
        ## Mean expression vector is simulated
        ## Mean and SD of prior are given
        model$muvec <- rnorm(p)*sigma + mu 
        model$marginalmean <- mu 
        model$marginalsd <- sqrt(sigma^2 + 1) 
    } else {
        stop("length of 'mu' should equal to 'p' or 1")
    }

    set.seed(m1+pi1+foldchange+dropout+p+mu+sigma)
    
    mat <- simulate(pfrac, nfrac, logfc, alpha, p, n, model, B, ncore)

    return(mat)
}






