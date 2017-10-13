#' @title Estimate the probability of success separation in two population model
#'
#' @param p1 Fraction of marker genes among total number genes in the novel subpopulation, 0.1 <= \code{p1} <= 1.6
#' @param n1 Fraction of cells in the novel subpopulation, 5 <= \code{n1} <= 10
#' @param foldchange Expression fold change in the novel subpopulation, 1.5 <= \code{foldchange} <= 8
#' @param dropout Overall dropout rate, 0.6 <= \code{dropout} <= 0.95
#' @param p Total number of genes in the simulation, 10000 <= \code{p}
#' @param n Vector of total number of cells, 100 <= min(n), max(n) <= 1000
#' @param mu Vector of mean expressions for \code{p} genes,
#' @param tau Standard deviation of mean expressions for \code{p} genes, 0 <= \code{tau}
#' @param perm Permute \code{mu}?
#' @param dfactor Multiplicative factor for separation criterion of two populations using first principal component scores, default is 2.
#' @param seed Random seed,
#' @param B Number of replicates to estiamte the probability of success separation, default is 200, 100 <= \code{B}
#' @param m Number of replicates to generated marginal distribution of expressions given \code{mu}, default is 50.
#' @param progress Print progress of simulations?
#' @param ncore Number of cores for use in the simulation using \code{\link{parallel}} package. default is 1, 1 <= \code{ncore} <= \code{parallel::detectCores()} 
#' @export
#' @return Vector of estiamted probabilities of same length as \code{n}
#' @examples \dontrun{
#' ## if mu is not specified, 
#' ncells(p1=1.6, n1=5, foldchange=4, dropout=.6)
#'
#' ## if mu is specified (for example as below)
#' p <- 26616
#' tau <- 5
#' mu <- rnorm(p)*tau
#' ncells(p1=1.6, n1=5, foldchange=4, dropout=.6, mu=mu)
#' }
ncells <- function(p1, n1, foldchange, dropout, p=26616, n=seq(100,1000,by=100), mu=NULL, tau=5, perm=TRUE, dfactor=2, seed=123, B=200, m=50, progress=TRUE, ncore=1)
{
    stopifnot(10000 <= p)
    stopifnot(0.1 <= p1 && p1 <= 1.6)
    stopifnot(100 <= min(n) && max(n) <= 1000)
    stopifnot(5 <= n1 && n1 <= 10)
    stopifnot(1.5 <= foldchange && foldchange <= 8)
    stopifnot(0.6 <= dropout && dropout <= 0.95)
    stopifnot(0 <= tau)
    stopifnot(100 <= B)
    stopifnot(1 <= ncore && ncore <= parallel::detectCores())

    pfrac <- p1/100
    nfrac <- n1/100
    logfc <- log2(foldchange)
    alpha <- dropout/(1-dropout)
    
    set.seed(seed)
    
    model <- list()
    if (is.null(mu)) { 
        model$parametric <- TRUE
        model$mu <- rnorm(p)*tau
        model$mean <- 0
        model$sd <- sqrt(tau^2 + 1)
        model$cdf <- NULL
        model$ncdf <- NULL
    } else {
        stopifnot(p == length(mu))
        model$parametric <- FALSE
        if (perm) {
            mu <- sample(mu) ## remove initial order
        }
        model$mu <- mu
        model$mean <- NULL
        model$sd <- NULL
        x <- matrix(rnorm(p*m), p, m) + mu
        marker <- 1:p <= p*p1/100
        cell <- 1:m <= m*n1/100
        x[marker, cell] <- x[marker, cell] + logfc
        model$cdf <- sort(x)
        model$ncdf <- length(x)
    }
    model$progress <- progress

    mat <- simulate(pfrac, nfrac, logfc, alpha, p, n, model, B, ncore)

    dpooled <- function(s0, s1, n0, n1) {
        return(ifelse(is.na(s1), s0, sqrt((s0^2*(n0 - 1) + s1^2*(n1 - 1))/(n0 + n1 - 2))))
    }
    res <- apply(mat, 3, function(x) mean(abs(x[,"m0"]-x[,"m1"]) > dfactor*dpooled(x[,"s0"], x[,"s1"], n-n1/100*n, n1/100*n)))

    return(res)
}






