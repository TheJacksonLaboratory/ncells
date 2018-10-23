#' @title Estimate the number of principal components in two population model
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
#' @export
#' @return Vector of estiamted probabilities of same length as \code{n}
#' @examples \dontrun{
#' ncells(m1=1.6, pi1=5, foldchange=4, dropout=.6, p=20000, mu=1, sigma=5)
#' }
ncomps <- function(m1, pi1, foldchange, dropout, p=26616, mu=1, sigma=5, seed=123)
{
    stopifnot(5000 <= p)
    stopifnot(0.1 <= m1 && m1 <= 10)
    stopifnot(5 <= pi1 && pi1 <= 50)
    stopifnot(1 <= foldchange && foldchange <= 32)
    stopifnot(0.6 <= dropout && dropout <= 0.999)
    stopifnot(0 <= mu)
    stopifnot(0 <= sigma)

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

    F <- function(y) (1-nfrac*pfrac)*pnorm(y,mean=mu,sd=sqrt(sigma^2+1)) + nfrac*pfrac*pnorm(y,mean=mu+logfc,sd=sqrt(sigma^2+1))
    f <- function(x) return(function(y) dnorm(y,mean=x,sd=1))
    g <- function(x) return(function(y) y*F(y)^alpha*f(x)(y))
    h <- function(x) return(function(y) y^2*F(y)^alpha*f(x)(y))

    mu0 <- model$muvec
    d0 <- numeric(p)
    v0 <- numeric(p)
    for (i in 1:p) {
        d0[i] <- integrate(g(mu0[i]),-Inf,Inf)$value
        v0[i] <- integrate(h(mu0[i]),-Inf,Inf)$value - d0[i]^2
    }

    mu1 <- mu0
    d1 <- d0
    v1 <- v0
    pos <- 1:round(p*pfrac)
    mu1[pos] <- mu1[pos] + logfc
    for (i in pos) {
        d1[i] <- integrate(g(mu1[i]),-Inf,Inf)$value
        v1[i] <- integrate(h(mu1[i]),-Inf,Inf)$value - d1[i]^2
    }

    lambda1 <- svd((1-nfrac)*nfrac*tcrossprod(d1[pos]-d0[pos])+diag((1-nfrac)*v0[pos]+nfrac*v1[pos]))$d[1]
    lambda2 <- (1-nfrac)*v0[-pos] + nfrac*v1[-pos]
    ncomp <- sum(lambda2 > lambda1) + 1

    return(ncomp)
}






