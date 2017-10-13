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
            ## x <- RcppZiggurat::zrnormLZLLV(n[i]*p)
            x <- nrorm(n[i]*p)
            keep <- integer(p)

            if (model$parametric) {
                .Call("dropout_with_centering", as.double(x), as.double(model$mu), as.integer(p), as.integer(sum(marker)), as.integer(n[i]), as.integer(sum(cell)), as.double(alpha), as.double(logfc), as.double(pfrac*nfrac), as.double(model$sd), as.integer(keep))
            } else {
                .Call("dropout_with_centering_cdf", as.double(x), as.double(model$mu), as.integer(p), as.integer(sum(marker)), as.integer(n[i]), as.integer(sum(cell)), as.double(alpha), as.double(logfc), as.double(model$cdf), as.integer(model$ncdf), as.integer(keep))
            }
            
            xmat <- matrix(x, n[i], p)[,keep > 0]

            l <- irlba::irlba(xmat, nu=1, nv=0)
            s <- l$u * l$d[1]
            x0 <- s[!cell]
            x1 <- s[cell]

            if (model$progress) {
                rb <- round(b/B*100)
                if (rb %in% seq(10,100,by=10)) {
                    str <- paste0(rb, "% is completed among B=", B, " iterations for n=", n[i])
                    print(str)
                }
            }

            return(c(mean(x0), mean(x1), sd(x0), sd(x1), ncol(xmat)))
        }

        res <- parallel::mclapply(1:B, loopfunc, mc.cores=ncore)

        estmat[,,i] <- matrix(unlist(res), ncol=5, byrow=TRUE)
    }

    return(estmat)
}




