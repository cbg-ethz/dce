#' Simulation study
#'
#' This function takes several parameters to define a simulation study.
#' @param nodes numebr of genes
#' @param samples vector of length two for sample number of both conditions
#' @param simruns simulations runs
#' @param mu mean expression
#' @param sd standard deviation
#' @param effRange vector of length four with lowest, second lowest, second
#'     highest and highest possible effect the effects are uniformly drawn from
#' @param truePos the fraction of differential effects; truePos=1 means all
#' are potentially differential (see effRange); if truePos=0.5, 50% of
#' differential effects are set to 0.
#' @param perturb positive or negative frantion of edges to be added or
#' removed, respectively
#' @param cormeth method for the correlation accuracy (see ?cor)
#' @param prob edge probability; either probability or "runif" to draw a
#' probability in each run
#' @param bootstrap can be either "none" (default) or include one or more
#' of the following: "basic", "full"
#' @param verbose verbose output, if TRUE
#' @param ... additional parameters for compute_differential_causal_effects
#' @author Martin Pirkl
#' @return accuracy for several different methods and statistics
#' for the ground truth as a list of two arrays
#' @export
#' @importFrom nem transitive.reduction
#' @importFrom Matrix rankMatrix
#' @examples
#' 1
simDce <- function(
    nodes=5, samples=c(10,10),simruns=10,mu=0,sd=1,
    effRange=c(-1,0,0,1),truePos=1,perturb=0,cormeth="p",
    prob="runif",bootstrap="none",verbose=FALSE, ...
    ) {
    cutoff <- 0.5
    getRates <- function(t, i, cutoff = 0.5) { # maybe AUC
        tp <- sum(abs(t) > cutoff & abs(i) > cutoff & sign(t) == sign(i))
        fn <- sum(abs(t) > cutoff & abs(i) <= cutoff)
        fp <- sum(abs(t) <= cutoff & abs(i) > cutoff)
        tn <- sum(abs(t) <= cutoff & abs(i) <= cutoff)
        return(c(tp, fp, tn, fn))
    }
    acc <- array(
        0, c(simruns,6,6),
        dimnames = list(
            runs = paste0("run_", seq_len(simruns)),
            methods = c(
                "random",
                "basic",
                "full",
                "simple correlation",
                "basic bootstrap",
                "full boostrap"
            ),
            metrics = c("correlation", "time", "tp", "fp", "tn", "fn")
        )
    )
    gtnfeat <- array(
        0, c(simruns, 4),
        dimnames = list(
            runs = paste0("run_", seq_len(simruns)),
            features = c(
                "avg children/parents",
                "max children",
                "max parents",
                "density"
            )
        )
    )
    n <- nodes
    m <- samples
    lB <- effRange[seq_len(2)]
    uB <- effRange[3:4]
    truepos <- truePos
    for (run in seq_len(simruns)) {
        if (prob %in% "runif") {
            p2 <- runif(1)
        } else {
            p2 <- prob
        }
        normal <- create_random_DAG(n, p2, lB, uB)
        tumor <- resample_edge_weights(normal, lB, uB, truepos)
        dn <- simulate_data(normal, m[2], normpars = c(mu,sd))
        dt <- simulate_data(tumor, m[1], normpars = c(mu,sd))
        gm <- as(normal, "matrix")
        gm[which(gm != 0)] <- 1
        cn <- trueEffects(normal)
        ct <- trueEffects(tumor)
        gtc <- nem::transitive.closure(gm, mat=TRUE)
        ## save features of gtn, which might correlate with accuracy:
        gtnfeat[run, 1] <- mean(apply(gm, 1, sum))
        gtnfeat[run, 2] <- max(apply(gm, 1, sum))
        gtnfeat[run, 3] <- max(apply(gm, 2, sum))
        gtnfeat[run, 4] <- sum(gm)
        ## ground truth:
        dcet <- (cn - ct)*gtc # gtn for differential causal effects
        dcegtn <- list(dce = dcet, graph = normal, dcefull = dcet)
        class(dcegtn) <- "dce"
        ## perturb network:
        if (perturb < 0) {
            adjn <- as(normal, "matrix")
            remedge <- sample(
                which(adjn != 0),
                floor(sum(adjn != 0)*abs(perturb))
            )
            adjn[remedge] <- 0
            adjn[which(adjn != 0)] <- 1
            normal <- tumor <- as(adjn, "graphNEL")
        }
        if (perturb > 0) {
            adjn <- as(normal, "matrix")
            addedge <- sample(
                which(adjn == 0 & upper.tri(adjn)),
                floor(sum(adjn == 0 & upper.tri(adjn))*perturb)
            )
            adjn[addedge] <- 1
            adjn[which(adjn != 0)] <- 1
            normal <- tumor <- as(adjn, "graphNEL")
        }
        ## full linear model:
        start <- as.numeric(Sys.time())
        dcei <- compute_differential_causal_effects(
            normal, dn,
            tumor, dt, method = "full", ...
        )
        acc[run, 3, 2] <- as.numeric(Sys.time()) - start
        dceifl <- dcei
        dcei <- dcei$dce
        acc[run, 3, 1] <- cor(
            as.vector(dcet), as.vector(dcei),
            method = cormeth, use = "complete.obs"
        )
        acc[run, 3, 3:6] <- getRates(dcet, dcei, cutoff = cutoff)
        if ("full" %in% bootstrap) {
            start <- as.numeric(Sys.time())
            dcei <- compute_differential_causal_effects(
                normal, dn,
                tumor, dt, method = "full",
                bootstrap = TRUE, ...
            )
            acc[run, 6, 2] <- as.numeric(Sys.time()) - start
            dceifl <- dcei
            dcei <- dcei$dce
            acc[run, 6, 1] <- cor(
                as.vector(dcet), as.vector(dcei),
                method = cormeth, use = "complete.obs"
            )
            acc[run, 6, 3:6] <- getRates(dcet, dcei, cutoff = cutoff)
        }
        ## normal
        Cn <- cov(dn)
        Ct <- cov(dt)
        if (
            Matrix::rankMatrix(Cn) == nrow(Cn) &
            Matrix::rankMatrix(Ct) == nrow(Ct)
        ) {
            start <- as.numeric(Sys.time())
            dcei <- compute_differential_causal_effects(
                normal, dn,
                tumor, dt, method = "normal", ...
            )
            acc[run, 2, 2] <- as.numeric(Sys.time()) - start
            dcein <- dcei
            dcei <- dcei$dce
            acc[run, 2, 1] <- cor(
                as.vector(dcet), as.vector(dcei),
                method = cormeth, use = "complete.obs"
            )
            acc[run, 2, 3:6] <- getRates(dcet, dcei, cutoff = cutoff)
            if ("basic" %in% bootstrap) {
                start <- as.numeric(Sys.time())
                dcei <- compute_differential_causal_effects(
                    normal, dn,
                    tumor, dt, method = "normal",
                    bootstrap = TRUE, ...
                )
                acc[run, 5, 2] <- as.numeric(Sys.time()) - start
                dcein <- dcei
                dcei <- dcei$dce
                acc[run, 5, 1] <- cor(
                    as.vector(dcet), as.vector(dcei),
                    method = cormeth, use = "complete.obs"
                )
                acc[run, 5, 3:6] <- getRates(dcet, dcei, cutoff = cutoff)
            }
        }
        ## simple correlation:
        start <- as.numeric(Sys.time())
        dcei <- (cor(dn) - cor(dt))*gtc
        acc[run, 4, 2] <- as.numeric(Sys.time()) - start
        dcec <- list(dce = dcei, graph = as(gtc, "graphNEL"), dcefull = dcei)
        dceic <- dcec
        class(dceic) <- "dce"
        dcec <- dcec$dce
        acc[run, 4, 1] <- cor(
            as.vector(dcet), as.vector(dcec),
            method = cormeth, use = "complete.obs"
        )
        acc[run, 4, 3:6] <- getRates(dcet, dcei, cutoff = cutoff)
        ## random base line:
        dcei <- dceicn <- dceict <- dcet
        start <- as.numeric(Sys.time())
        dcei[which(gtc != 0)] <-
            runif(sum(gtc != 0), lB[1], uB[2]) -
            runif(sum(gtc != 0), lB[1], uB[2])
        acc[run, 1, 2] <- as.numeric(Sys.time()) - start
        dcer <- list(dce = dcei, graph = as(gtc, "graphNEL"), dcefull = dcei)
        dceir <- dcer
        class(dceir) <- "dce"
        dcer <- dcer$dce
        coridx <- which(dcet != 0 | dcer != 0)
        acc[run, 1, 1] <- cor(
            as.vector(dcet), as.vector(dcer),
            method = cormeth, use = "complete.obs"
        )
        acc[run, 1, 3:6] <- getRates(dcet, dcei, cutoff = cutoff)
        if (verbose) {
            cat(paste0(run, "."))
        }
    }
    dceSim <- list(acc=acc,gtnFeat=gtnfeat)
    class(dceSim) <- "dceSim"
    return(dceSim)
}
#' Plot simulation study
#'
#' Takes a dceSim object and produces a figure.
#' @param x dceSim object
#' @param col `col` argument passed to `boxplot`
#' @param showMeth which method results to plot
#' @param showFeat which features to plot
#' @param methNames method plot labels
#' @param ... additional arguments for boxplot
#' @author Martin Pirkl
#' @method plot dceSim
#' @return plot
#' @export
#' @importFrom nem transitive.reduction
#' @importFrom graphics par boxplot axis
plot.dceSim <- function(
    x, col = seq_len(4),
    showMeth = seq_len(4), showFeat = 1, methNames = NULL,
    ...
) {
    runs <- dim(x$acc)[1]
    if (is.null(methNames)) {
        methNames <- dimnames(x$acc)[[2]][showMeth]
    }
    par(mfrow=c(1,length(showFeat)))
    if (1 %in% showFeat) {
        boxplot(
            x$acc[seq_len(runs), showMeth, 1], col = col,
            main="Correlation", , xaxt = "n", ...
        )
        axis(1, seq_len(length(showMeth)), labels = methNames)
    }
    if (4 %in% showFeat) {
        tmp <- x$acc[seq_len(runs), showMeth, 3]/
            (x$acc[seq_len(runs), showMeth, 3]+
             x$acc[seq_len(runs), showMeth, 4])
        boxplot(
            tmp, col = col,
            main="PPV", , xaxt = "n", ...
        )
        axis(1, seq_len(length(showMeth)), labels = methNames)
    }
    if (5 %in% showFeat) {
        tmp <- x$acc[seq_len(runs), showMeth, 5]/
            (x$acc[seq_len(runs), showMeth, 5]+
             x$acc[seq_len(runs), showMeth, 6])
        boxplot(
            tmp, col = col,
            main="NPV", , xaxt = "n", ...
        )
        axis(1, seq_len(length(showMeth)), labels = methNames)
    }
    if (2 %in% showFeat) {
        boxplot(
            x$acc[seq_len(runs), showMeth, 2], col = col,
            main="Time", , xaxt = "n", ...
        )
        axis(1, seq_len(length(showMeth)), labels = methNames)
    }
    if (3 %in% showFeat) {
        boxplot(
            x$gtnFeat[seq_len(runs), seq_len(4)], col = col,
            main="Ground truth Features", , xaxt = "n", ...
        )
        axis(
            1, seq_len(4),
            labels = dimnames(x$gtnFeat)[[2]][seq_len(4)]
        )
    }
}
