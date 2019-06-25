library(pcalg)

## add noise options to the pcalg data generating function:

rmvDAG_2 <- function (n, dag, errDist = c("normal", "cauchy", "t4", "mix",
    "mixt3", "mixN100"), normpars = c(0,1), mix = 0.1, errMat = NULL, back.compatible = FALSE,
    use.node.names = !back.compatible) {
    stopifnot(is(dag, "graph"), (p <- length(nodes(dag))) >=
        2)
    weightMatrix <- if (back.compatible)
        wgtMatrix.0(dag)
    else wgtMatrix(dag)
    nonZeros <- which(weightMatrix != 0, arr.ind = TRUE)
    if (nrow(nonZeros) > 0) {
        if (any(nonZeros[, 1] - nonZeros[, 2] < 0) || any(diag(weightMatrix) !=
            0))
            stop("Input DAG must be topologically ordered!")
    }
    errDist <- match.arg(errDist)
    if (grepl("^mix", errDist))
        eMat <- function(outs) {
            X <- c(rnorm(n * p - length(outs)), outs)
            matrix(sample(X), nrow = n)
        }
    if (is.null(errMat)) {
        errMat <- switch(errDist, normal = matrix(rnorm(n * p, normpars[1], normpars[2]),
            nrow = n), cauchy = matrix(rcauchy(n * p), nrow = n),
            t4 = matrix(rt(n * p, df = 4), nrow = n), mix = eMat(rcauchy(round(mix *
                n * p))), mixt3 = eMat(rt(round(mix * n * p),
                df = 3)), mixN100 = eMat(rnorm(round(mix * n *
                p), sd = 10)))
    }
    else {
        stopifnot(!is.null(dim.eM <- dim(errMat)), dim.eM ==
            c(n, p), is.numeric(errMat))
    }
    if (use.node.names)
        colnames(errMat) <- nodes(dag)
    if (sum(weightMatrix) > 0) {
        X <- errMat
        for (j in 2:p) {
            ij <- 1:(j - 1)
            X[, j] <- X[, j] + X[, ij, drop = FALSE] %*% weightMatrix[j,
                ij]
        }
        X
    }
    else errMat
}

## put noise (here resample) on the weights of a dag (for tumor samples):

newWeights <- function(g, lB = -1, uB = 1) {
    n <- length(g@nodes)
    w <- g@edgeData@data
    for (i in 1:length(w)) {
        w[[i]]$weight <- runif(1, lB, uB)
    }
    g@edgeData@data <- w
    return(g)
}

## parameters for simulations:

runs <- 100
p <- 0.2 # edge prob
## uniform limits:
lB <- -1
uB <- 1
##
n <- 10 # nodes
m <- 100 # samples
sd <- 1 # standard deviation

acc <- matrix(0, runs, 4)

for (i in 1:runs) {

    normal <- randomDAG(n, p, lB, uB)

    dn <- rmvDAG_2(m, normal, normpars = c(0,sd))

    tumor <- newWeights(normal) # resample edge weights

    dt <- rmvDAG_2(m, tumor, normpars = c(0,sd))

    cn <- trueCov(normal)

    ct <- trueCov(tumor)

    gm <- as(normal, "matrix")
    gm[which(gm != 0)] <- 1

    gtc <- mnem:::mytc(gm) # transitively closed graph as matrix
    diag(gtc) <- 0

    dce <- (cn - ct)*gtc # gtn for differential causal effects

    dcei <- dce*0 # inferred

    for (j in 1:n) {

        dcei[j, ] <- idaFast(j, 1:n, cov(dn), normal) - idaFast(j, 1:n, cov(dt), tumor)

    }

    dcei <- dcei*gtc

    ## what as accuracy? Acount for random chance.

    acc[i, 1] <- cor(as.vector(dce), as.vector(dcei), method = "s")

    acc[i, 2] <- dist(matrix(c(dce, dcei), 2))

}

## comparative plot:

efreq <- round(dce[which(dce != 0)], 2)
efreqi <- round(dcei[which(dcei != 0)], 2)
efreq[abs(efreq) > 2] <- 2 # even though the gtn does not allow for dce > 2, noise can pass that limit
efreqi[abs(efreqi) > 2] <- 2

par(mfrow=c(1,2))
dnf <- mnem:::adj2dnf(gtc)
dnf <- dnf[grep("=", dnf)] # this preprocessing will be implemented in plotDnf (bug)
mnem::plotDnf(dnf, labels = efreq, edgecol = rgb(abs(efreq)/2,0,(2-abs(efreq))/2))
mnem::plotDnf(dnf, labels = efreqi, edgecol = rgb(abs(efreqi)/2,0,(2-abs(efreqi))/2))
