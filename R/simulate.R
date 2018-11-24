simulate <- function(pfrac, nfrac, logfc, alpha, p, n, model, B, ncore)
{
    estmat <- array(, dim=c(B, 5, length(n)))
    dimnames(estmat) <- list(NULL, c("m0", "m1", "s0", "s1", "detect"), paste("n=", n))

    for (i in 1:length(n)) {
        pos <- list()
        marker <- 1:p <= pfrac*p
        cell <- 1:n[i] <= nfrac*n[i]

        loopfunc <- function(b)
        {
            x <- rnorm(n[i]*p) ## random N(0,1) error
            keep <- integer(p)

            .Call("dropout_with_centering", as.double(x), as.double(model$muvec), as.integer(p), as.integer(sum(marker)), as.integer(n[i]), as.integer(sum(cell)), as.double(alpha), as.double(logfc), as.double(pfrac*nfrac), as.double(model$marginalsd), as.double(model$marginalmean), as.integer(keep))
            
            xmat <- t(matrix(x, n[i], p)[,keep > 0])
            l <- irlba::irlba(xmat, 1, right_only=TRUE, maxit=10000)
            s <- l$v * l$d[1]
            x0 <- s[!cell]
            x1 <- s[cell]

            return(c(mean(x0), mean(x1), sd(x0), sd(x1), sum(keep > 0)))
        }
        res <- parallel::mclapply(1:B, loopfunc, mc.cores=ncore)
        estmat[,,i] <- matrix(unlist(res), ncol=5, byrow=TRUE)
    }
    return(estmat)
}




