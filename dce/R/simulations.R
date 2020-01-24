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
#' @param bootstrap can be either 0 (default) or have the first argument as
#' the number of runs and include one or more
#' of the following: "basic", "full"
#' @param verbose verbose output, if TRUE
#' @param test if greater than 0, tests the pathway for enrichment with
#' test permutation runs
#' @param enriched fraction of runs with enriched pathways
#' @param theta
#' @param partial
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
    prob="runif",bootstrap=0,verbose=FALSE,test=0,
    enriched=1,theta=NULL,partial=TRUE,...
    ) {
    bootruns <- as.numeric(bootstrap[1])
    if (bootruns == 0) {
        bootstrap <- 0
    }
    cutoff <- 0.75
    acc <- array(
        0, c(simruns,4,8),
        dimnames = list(
            runs = paste0("run_", seq_len(simruns)),
            methods = c(
                "random",
                "NB",
                "simple correlation",
                "NB bootstrap"
            ),
            metrics = c("correlation", "time",
                        "tp", "fp", "tn", "fn",
                        "p-value (perm)", "AUC")
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
        if (run <= floor(simruns*enriched)) {
            truepos <- 1
        } else {
            truepos <- 0
        }
        tumor <- resample_edge_weights(normal, lB, uB, truepos)
        dn <- simulate_data(normal, m[2], mu, sd)
        dt <- simulate_data(tumor, m[1], mu, sd)
        gm <- as(normal, "matrix")
        gm[which(gm != 0)] <- 1
        cn <- trueEffects(normal, partial = partial)
        ct <- trueEffects(tumor, partial = partial)
        gtc <- nem::transitive.closure(gm, mat=TRUE)
        ## save features of gtn, which might correlate with accuracy:
        gtnfeat[run, 1] <- mean(apply(gm, 1, sum))
        gtnfeat[run, 2] <- max(apply(gm, 1, sum))
        gtnfeat[run, 3] <- max(apply(gm, 2, sum))
        gtnfeat[run, 4] <- sum(gm)
        ## ground truth:
        if (partial) {
            dcet <- (cn - ct)*gm
        } else {
            dcet <- (cn - ct)*gtc
        }
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
        ## negative binomial:
        coracc <- function(x, y, method = cormeth) {
            x <- as.vector(x)
            y <- as.vector(y)
            idx <- which(x != 0 | y != 0)
            a <- cor(x[idx], y[idx], method = cormeth)
            return(a)
        }
        computeAUC <- function(x, y, errDist = errDist) {
            x <- abs(x)
            y <- abs(y)
            cl <- 100
            ppv <- rec <- numeric(cl)
            cs <- seq(max(max(x), max(y)), 0, length.out = cl)
            for (c in seq_len(cl)) {
                tp <- sum(y > cs[c] & x > cs[c])
                fp <- sum(y > cs[c] & x <= cs[c])
                fn <- sum(y <= cs[c] & x > cs[c])
                ppv[c] <- tp/(tp+fp)
                rec[c] <- tp/(tp+fn)
            }
            ppv[is.na(ppv)] <- 0
            rec[is.na(rec)] <- 0
            ppv <- ppv[order(rec)]
            rec <- rec[order(rec)]
            auc <- function(a,b) {
                n <- length(a)
                c <- sum((a-c(0,a[-n]))*(b+c(0,b[-n]))/c(1,rep(2,n-1)))
                return(c)
            }
            AUC <- auc(rec, ppv)
            return(AUC)
        }
        ## nbinom:
        start <- as.numeric(Sys.time())
        dcei <- compute_differential_causal_effects(
            normal, dn,
            tumor, dt, method = "full",
            partial = partial, ...
        )
        acc[run, 2, 2] <- as.numeric(Sys.time()) - start
        dceifl <- dcei
        dcei <- dcei$dce
        acc[run, 2, 1] <- coracc(dcet, dcei)
        acc[run, 2, 8] <- computeAUC(dcet, dcei)
        acc[run, 2, 3:6] <- as.numeric(get_prediction_counts(dcet, dcei, cutoff = cutoff))
        if (test > 0) {
            p.tmp <- compute_enrichment(normal, dn, dt, permutation_count = test, partial = partial)
            acc[run, 2, 7] <- p.tmp[[1]]
        }
        if ("NB" %in% bootstrap) {
            start <- as.numeric(Sys.time())
            dcei <- compute_differential_causal_effects(
                normal, dn,
                tumor, dt, method = "full",
                bootstrap = TRUE, runs = bootruns,
                errDist = errDist, partial = partial, ...
            )
            acc[run, 4, 2] <- as.numeric(Sys.time()) - start
            dceifl <- dcei
            dcei <- dcei$dce
            acc[run, 4, 1] <- coracc(dcet, dcei)
            acc[run, 4, 8] <- computeAUC(dcet, dcei)
            acc[run, 4, 3:6] <- as.numeric(get_prediction_counts(dcet, dcei, cutoff = cutoff))
        }
        ## simple correlation:
        start <- as.numeric(Sys.time())
        if (partial) {
            dcei <- (cor(dn) - cor(dt))*gtc
        } else {
            dcei <- (cor(dn) - cor(dt))*gtc
        }
        acc[run, 3, 2] <- as.numeric(Sys.time()) - start
        dcec <- list(dce = dcei, graph = as(gtc, "graphNEL"), dcefull = dcei)
        dceic <- dcec
        class(dceic) <- "dce"
        dcec <- dcec$dce
        acc[run, 3, 1] <- coracc(dcet, dcec)
        acc[run, 3, 8] <- computeAUC(dcet, dcec)
        acc[run, 3, 3:6] <- as.numeric(get_prediction_counts(dcet, dcei, cutoff = cutoff))
        ## random base line:
        dcei <- dceicn <- dceict <- dcet
        start <- as.numeric(Sys.time())
        if (partial) {
            dcei[which(gm != 0)] <-
                runif(sum(gm != 0), lB[1], uB[2]) -
                runif(sum(gm != 0), lB[1], uB[2])
        } else {
            dcei[which(gtc != 0)] <-
                runif(sum(gtc != 0), lB[1], uB[2]) -
                runif(sum(gtc != 0), lB[1], uB[2])
        }
        acc[run, 1, 2] <- as.numeric(Sys.time()) - start
        dcer <- list(dce = dcei, graph = as(gtc, "graphNEL"), dcefull = dcei)
        dceir <- dcer
        class(dceir) <- "dce"
        dcer <- dcer$dce
        coridx <- which(dcet != 0 | dcer != 0)
        acc[run, 1, 1] <- coracc(dcet, dcei)
        acc[run, 1, 8] <- computeAUC(dcet, dcei)
        acc[run, 1, 3:6] <- as.numeric(get_prediction_counts(dcet, dcei, cutoff = cutoff))
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
#' @param border `border` argument passed to `boxplot`
#' @param showMeth which method results to plot
#' @param showFeat which features to plot
#' @param methNames method plot labels
#' @param ... additional arguments for boxplot
#' @author Martin Pirkl
#' @method plot dceSim
#' @return plot
#' @export
#' @importFrom graphics par boxplot axis
plot.dceSim <- function(
    x, col = seq_len(4), border = col,
    showMeth = seq_len(4), showFeat = 1, methNames = NULL,
    ...
    ) {
    if (length(col) < length(showMeth)) {
        col <- c(col, col)
        border <- c(border, border)
    }
    runs <- dim(x$acc)[1]
    if (is.null(methNames)) {
        methNames <- dimnames(x$acc)[[2]][showMeth]
    }
    par(mfrow=c(1,length(showFeat)))
    if (1 %in% showFeat) {
        boxplot(
            x$acc[seq_len(runs), showMeth, 1], col = col, border = border,
            main="Correlation", xaxt = "n", ...
        )
        axis(1, seq_len(length(showMeth)), labels = methNames)
    }
    if (4 %in% showFeat) {
        tmp <- x$acc[seq_len(runs), showMeth, 3]/
            (x$acc[seq_len(runs), showMeth, 3]+
             x$acc[seq_len(runs), showMeth, 4])
        boxplot(
            tmp, col = col, border = border,
            main="PPV", xaxt = "n", ...
        )
        axis(1, seq_len(length(showMeth)), labels = methNames)
    }
    if (5 %in% showFeat) {
        tmp <- x$acc[seq_len(runs), showMeth, 5]/
            (x$acc[seq_len(runs), showMeth, 5]+
             x$acc[seq_len(runs), showMeth, 6])
        boxplot(
            tmp, col = col, border = border,
            main="NPV", xaxt = "n", ...
        )
        axis(1, seq_len(length(showMeth)), labels = methNames)
    }
    if (2 %in% showFeat) {
        boxplot(
            x$acc[seq_len(runs), showMeth, 2], col = col, border = border,
            main="Time", xaxt = "n", ...
        )
        axis(1, seq_len(length(showMeth)), labels = methNames)
    }
    if (3 %in% showFeat) {
        boxplot(
            x$gtnFeat[seq_len(runs), seq_len(4)], col = col, border = border,
            main="Ground truth Features", xaxt = "n", ...
        )
        axis(
            1, seq_len(4),
            labels = dimnames(x$gtnFeat)[[2]][seq_len(4)]
        )
    }
    if (6 %in% showFeat) {
        showMeth2 <- showMeth[which(showMeth %in% c(3,4))]
        boxplot(
            x$acc[seq_len(runs), showMeth2, 7],
            col = col[which(showMeth %in% c(3,4))],
            border = border[which(showMeth %in% c(3,4))],
            main="empirical p-value", xaxt = "n", ylim = c(0,1), ...
        )
        abline(h=c(0.1,0.05,0.01), col = rgb(0.5,0.5,0.5,0.75), lty = 2)
        axis(1, seq_len(length(showMeth2)),
             labels = methNames[which(showMeth %in% c(3,4))])
    }
    if (7 %in% showFeat) {
        boxplot(
            x$acc[seq_len(runs), showMeth, 8], col = col, border = border,
            main="AUC", xaxt = "n", ...
        )
        axis(1, seq_len(length(showMeth)), labels = methNames)
    }
}
