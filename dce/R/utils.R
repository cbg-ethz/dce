#' Modified theta.ml function from package MASS
#'
#' fixes a bug, if theta estimation breaks
#' see ?MASS::theta.ml for argument values
theta.ml.rob <- function (y, mu, n = sum(weights), weights, limit = 10, eps = .Machine$double.eps^0.25,
    trace = FALSE)
{
    score <- function(n, th, mu, y, w) sum(w * (digamma(th +
        y) - digamma(th) + log(th) + 1 - log(th + mu) - (y +
        th)/(mu + th)))
    info <- function(n, th, mu, y, w) sum(w * (-trigamma(th +
        y) + trigamma(th) - 1/th + 2/(mu + th) - (y + th)/(mu +
        th)^2))
    if (inherits(y, "lm")) {
        mu <- y$fitted.values
        y <- if (is.null(y$y))
            mu + residuals(y)
        else y$y
    }
    if (missing(weights))
        weights <- rep(1, length(y))
    t0 <- n/sum(weights * (y/mu - 1)^2)
    it <- 0
    del <- 1
    if (trace)
        message(sprintf("theta.ml.rob: iter %d 'theta = %f'", it,
            signif(t0)), domain = NA)
    while ((it <- it + 1) < limit && abs(del) > eps) {
        t0 <- abs(t0)
        del <- score(n, t0, mu, y, weights)/(i <- info(n, t0,
            mu, y, weights))
        t0 <- t0 + del
        if (trace)
            message("theta.ml.rob: iter", it, " theta =", signif(t0))
        if (is.nan(del) | is.infinite(del)) {
            warning("NaNs produced. Resetting theta to 0.")
            t0 <- 0
            del <- 0
        }
    }
    if (t0 < 0) {
        t0 <- 0
        warning("estimate truncated at zero")
        attr(t0, "warn") <- gettext("estimate truncated at zero")
    }
    if (it == limit) {
        warning("iteration limit reached")
        attr(t0, "warn") <- gettext("iteration limit reached")
    }
    attr(t0, "SE") <- sqrt(1/i)
    t0
}
#' Modified glm.nb function from package MASS
#'
#' fixes a bug, if theta estimation breaks
#' see ?MASS::glm.nb for argument values
glm.nb.rob <- function (formula, data, weights, subset, na.action, start = NULL,
    etastart, mustart, control = glm.control(...), method = "glm.fit",
    model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, ...,
    init.theta, link = log)
{
    loglik <- function(n, th, mu, y, w) sum(w * (lgamma(th +
        y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y *
        log(mu + (y == 0)) - (th + y) * log(th + mu)))
    link <- substitute(link)
    fam0 <- if (missing(init.theta))
        do.call("poisson", list(link = link))
    else do.call("negative.binomial", list(theta = init.theta,
        link = link))
    mf <- Call <- match.call()
    m <- match(c("formula", "data", "subset", "weights", "na.action",
        "etastart", "mustart", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval.parent(mf)
    Terms <- attr(mf, "terms")
    if (method == "model.frame")
        return(mf)
    Y <- model.response(mf, "numeric")
    X <- if (!is.empty.model(Terms))
        model.matrix(Terms, mf, contrasts)
    else matrix(, NROW(Y), 0)
    w <- model.weights(mf)
    if (!length(w))
        w <- rep(1, nrow(mf))
    else if (any(w < 0))
        stop("negative weights not allowed")
    offset <- model.offset(mf)
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    n <- length(Y)
    if (!missing(method)) {
        if (!exists(method, mode = "function"))
            stop(gettextf("unimplemented method: %s", sQuote(method)),
                domain = NA)
        glm.fitter <- get(method)
    }
    else {
        method <- "glm.fit"
        glm.fitter <- stats::glm.fit
    }
    if (control$trace > 1)
        message("Initial fit:")
    fit <- glm.fitter(x = X, y = Y, w = w, start = start, etastart = etastart,
        mustart = mustart, offset = offset, family = fam0, control = list(maxit = control$maxit,
            epsilon = control$epsilon, trace = control$trace >
                1), intercept = attr(Terms, "intercept") > 0)
    class(fit) <- c("glm", "lm")
    mu <- fit$fitted.values
    th <- as.vector(theta.ml.rob(Y, mu, sum(w), w, limit = control$maxit,
                                   trace = control$trace > 2))
    if (th < 1) { th <- 1 }
    if (control$trace > 1)
        message(gettextf("Initial value for 'theta': %f", signif(th)),
              domain = NA)
    negative.binomial <- MASS::negative.binomial
    fam <- do.call("negative.binomial", list(theta = th, link = link))
    iter <- 0
    d1 <- sqrt(2 * max(1, fit$df.residual))
    d2 <- del <- 1
    g <- fam$linkfun
    Lm <- loglik(n, th, mu, Y, w)
    Lm0 <- Lm + 2 * d1
    while ((iter <- iter + 1) <= control$maxit && (abs(Lm0 -
        Lm)/d1 + abs(del)/d2) > control$epsilon) {
        eta <- g(mu)
        fit <- glm.fitter(x = X, y = Y, w = w, etastart = eta,
            offset = offset, family = fam, control = list(maxit = control$maxit,
                epsilon = control$epsilon, trace = control$trace >
                  1), intercept = attr(Terms, "intercept") >
                0)
        t0 <- th
        th <- theta.ml.rob(Y, mu, sum(w), w, limit = control$maxit,
            trace = control$trace > 2)
        fam <- do.call("negative.binomial", list(theta = th,
            link = link))
        mu <- fit$fitted.values
        del <- t0 - th
        Lm0 <- Lm
        Lm <- loglik(n, th, mu, Y, w)
        if (control$trace) {
            Ls <- loglik(n, th, Y, Y, w)
            Dev <- 2 * (Ls - Lm)
            message(sprintf("Theta(%d) = %f, 2(Ls - Lm) = %f",
                iter, signif(th), signif(Dev)), domain = NA)
        }
        if (is.nan(Lm)) {
            th <- t0
            Lm <- Lm0
            warning("no convergence: try larger number of iterations or larger epsilon")
            break()
        }
    }
    if (!is.null(attr(th, "warn")))
        fit$th.warn <- attr(th, "warn")
    if (iter > control$maxit) {
        warning("alternation limit reached")
        fit$th.warn <- gettext("alternation limit reached")
    }
    if (length(offset) && attr(Terms, "intercept")) {
        null.deviance <- if (length(Terms))
            glm.fitter(X[, "(Intercept)", drop = FALSE], Y, w,
                offset = offset, family = fam, control = list(maxit = control$maxit,
                  epsilon = control$epsilon, trace = control$trace >
                    1), intercept = TRUE)$deviance
        else fit$deviance
        fit$null.deviance <- null.deviance
    }
    class(fit) <- c("negbin", "glm", "lm")
    fit$terms <- Terms
    fit$formula <- as.vector(attr(Terms, "formula"))
    Call$init.theta <- signif(as.vector(th), 10)
    Call$link <- link
    fit$call <- Call
    if (model)
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x)
        fit$x <- X
    if (!y)
        fit$y <- NULL
    fit$theta <- as.vector(th)
    fit$SE.theta <- attr(th, "SE")
    fit$twologlik <- as.vector(2 * Lm)
    fit$aic <- -fit$twologlik + 2 * fit$rank + 2
    fit$contrasts <- attr(X, "contrasts")
    fit$xlevels <- .getXlevels(Terms, mf)
    fit$method <- method
    fit$control <- control
    fit$offset <- offset
    fit
}
#' permutation test for (partial) correlation on non-Gaussian data
#' @param x wild type data set
#' @param y mutant data set
#' @param i number of iterations (permutations)
#' @param fun function to compute the statistic, e.g., cor or pcor
pcor_perm <- function(x, y, iter = 1000, fun = pcor, method = "spearman") {
    z <- fun(y) - fun(x)
    p <- z*0
    for (i in seq_len(iter)) {
        xp <- x[, sample(seq_len(ncol(x)), ncol(x))]
        yp <- y[, sample(seq_len(ncol(y)), ncol(y))]
        zp <- fun(yp) - fun(xp)
        idx <- which(abs(zp) >= abs(z))
        p[idx] <- p[idx] + 1
    }
    p <- p/iter
    return(p)
}
#' graphNEL with 0 edge weights to proper adjacency matrix
#'
#' @param g graphNEL object
#' @export
as.adjmat <- function(g) {
    adj <- as(g, "matrix")
    for (p in names(g@edgeData@data)) {
        a <- gsub("\\|.*", "", p)
        b <- gsub(".*\\|", "", p)
        adj[a, b] <- 1
    }
    return(adj)
}

#' robust partial correlation of column variables
#'
#' @param x matrix
#' @importFrom ppcor pcor
#' @export
pcor <- function(x, method = "spearman") {
    rho <- try(ppcor::pcor(x, method = method), silent = TRUE)
    if (length(grep("Error", rho)) > 0) {
        warning(paste0("Moore-Penrose generalized matrix invers in ",
                       "function ppcor::pcor crashed. Using ",
                       "matlib::Ginv instead."))
        omega <- cor(x, method = method)
        p <- Gsolve(omega)
        pdiag <- diag(p)%*%t(diag(p))
        rho <- -p/(pdiag^0.5)
        diag(rho) <- 1
        rho[which(is.na(rho) | is.infinite(rho))] <- 0
    } else {
        rho <- rho$estimate
    }
    return(rho)
}
#' Graph to DAG
#'
#' Converts a general graph to a dag with minimum distance to the
#' original graph. The general idea is to transitively close the graph
#' to detect cycles and remove them based on the rule "the more outgoing
#' edges a node has, the more likely it is that incoming edges from a
#' cycle will be deleted, and vice versa. However, this is too rigorous
#' and deletes too many edges, which do not lead to a cycle. These
#' edges are added back in the final step.
#' @param g graph as adjacency matrix
#' @param tc if TRUE computes the transitive closure
#' @author Ken Adams
#' @return dag as adjacency matrix
#' @export
#' @examples
#' g <- matrix(c(1,0,1,0, 1,1,0,0, 0,1,1,0, 1,1,0,1), 4, 4)
#' rownames(g) <- colnames(g) <- LETTERS[1:4]
#' dag <- g2dag(g)
g2dag <- function(g, tc = FALSE) {
    ord <- order(apply(g, 1, sum) - apply(g, 2, sum), decreasing = 1)
    g <- g[ord, ord]
    cyc <- intersect(which(g+t(g) > 1), which(lower.tri(g) == TRUE))
    g[cyc] <- 0
    diag(g) <- 1
    g2 <- g
    g0 <- g*0
    epiNEM::HeatmapOP(g, Colv = 0, Rowv = 0)
    for (i in seq_len(nrow(g))) {
        g2 <- g2%*%g
        g2[which(g2 > 0)] <- 1
        ord <- order(apply(g2, 1, sum) - apply(g2, 2, sum), decreasing = 1)
        g2 <- g2[ord, ord]
        g <- g[ord, ord]
        cyc <- intersect(which(g2+t(g2) > 1), which(lower.tri(g2) == TRUE))
        g2[cyc] <- 0
        if (all(g0 == g2)) { break() }
        g0 <- g2
        epiNEM::HeatmapOP(g2, Colv = 0, Rowv = 0)
    }
    if (!tc) {
        print(g)
        g3 <- g2*0
        g3[which(g2 == 1 & g == 1)] <- 1
        for (e in which(g == 1 & g3 == 0)) {
            for (a in 1:ncol(g)) {
                b <- e - a*nrow(g)
                if (b < 0) {
                    b <- e - (a-1)*nrow(g)
                    break()
                }
            }
            print(e)
            print(b)
            print(a)
            if (sum(g[b, ]) >= sum(g[a, ])) {
                g3[b, a] <- 1
            }
            g4 <- mytc(g3)
            ord <- order(apply(g4, 1, sum) - apply(g4, 2, sum), decreasing = 1)
            g4 <- g4[ord, ord]
            if (any(g2+t(g2) > 1)) {
                g3[a, b] <- 0
            }
        }
        g2 <- g3
    }
    return(g2)
}
#' @importFrom matlib Ginv
#' @noRd
Gsolve <- function(a, b=NULL, ...) {
    if (!is.null(b)) {
        x <- t(solve(diag(1, nrow(a)), matlib::Ginv(a, ...)%*%b))
    } else {
        x <- matlib::Ginv(a, ...)
    }
    return(x)
}
#' Compute true positive/... counts
#'
#' Useful for performance evaluations
#' @param truth Ground truth
#' @param inferred Computed results
#' @param cutoff Threshold for classification
#' @author Hans Wurst
#' @return data.frame
#' @export
#' @examples
#' get_prediction_counts(c(1,0), c(1,1))
get_prediction_counts <- function(truth, inferred, cutoff = 0.5) {
    tp <- sum(
        abs(truth) > cutoff &
        abs(inferred) > cutoff &
        sign(truth) == sign(inferred) & inferred != 0,
        na.rm = TRUE
    )
    fn <- sum(
        abs(truth) > cutoff &
        abs(inferred) <= cutoff &
        inferred != 0,
        na.rm = TRUE
    )
    fp <- sum(
        abs(truth) <= cutoff &
        abs(inferred) > cutoff |
        (
            abs(truth) > cutoff &
            abs(inferred) > cutoff &
            sign(truth) != sign(inferred)
        ) &
        inferred != 0,
        na.rm = TRUE
    )
    tn <- sum(
        abs(truth) <= cutoff &
        abs(inferred) <= cutoff &
        inferred != 0,
        na.rm = TRUE
    )

    return(data.frame(
        true.positive=tp,
        false.positive=fp,
        true.negative=tn,
        false.negative=fn
    ))
}
#' Create random DAG (topologically ordered)
#'
#' Creates a DAG according to given parameters.
#' @param n Number of nodes
#' @param prob Probability of creating an edge
#' @param lB Lower bound for edge weights
#' @param uB Upper bound for edge weights
#' @param node.labels Node labels
#' @author Martin Pirkl
#' @return graph
#' @export
#' @importFrom methods new
#' @examples
#' dag <- create_random_DAG(30, 0.2)
create_random_DAG <- function (
    n, prob,
    lB = -1, uB = 1,
    node.labels = paste0("n", as.character(seq_len(n)))
) {
    stopifnot(
        n >= 2, is.numeric(prob), length(prob) == 1, 0 <= prob, prob <= 1,
        is.numeric(lB), is.numeric(uB), lB <= uB
    )
    if (length(uB) == 1) { uB <- c(0, uB) }
    if (length(lB) == 1) { lB <- c(lB, 0) }
    if (lB[1] - lB[2] == 0) {
        lB <- uB
    }
    if (uB[1] - uB[2] == 0) {
        uB <- lB
    }
    edL <- vector("list", n)
    nmbEdges <- 0L
    for (i in seq_len(n - 2)) {
        listSize <- rbinom(1, n - i, prob)
        nmbEdges <- nmbEdges + listSize
        edgeList <- sample(seq(i + 1, n), size = listSize)
        posedges <- sample(
            seq_len(length(edgeList)), ceiling(length(edgeList)/2)
        )
        negedges <- seq_len(length(edgeList))[-posedges]
        weightList <- numeric(length(edgeList))
        weightList[posedges] <- runif(
            length(posedges), min = uB[1], max = uB[2]
        )
        if (length(negedges) > 0) {
            weightList[negedges] <- runif(
                length(negedges), min = lB[1], max = lB[2]
            )
        }
        edL[[i]] <- list(edges = edgeList, weights = weightList)
    }
    listSize <- rbinom(1, 1, prob)
    if (listSize > 0) {
        nmbEdges <- nmbEdges + 1
        edgeList <- n
        weightList <- sample(c(runif(1, min = lB[1], max = lB[2]), runif(1, min = uB[1], max = uB[2])), 1)
    }
    else {
        edgeList <- integer(0)
        weightList <- numeric(0)
    }
    edL[[n - 1]] <- list(edges = edgeList, weights = weightList)
    if (nmbEdges > 0) {
        edL[[n]] <- list(edges = integer(0), weights = numeric(0))
        names(edL) <- node.labels
        new("graphNEL", nodes = node.labels, edgeL = edL, edgemode = "directed")
    }
    else new("graphNEL", nodes = node.labels, edgemode = "directed")
}
#' Resample network edge weights
#'
#' Takes a graph and modifies edge weights.
#' @param g original graph
#' @param tp fraction of edge weights which will be modified
#' @param mineff minimal differential effect size
#' @param maxeff maximum effect effect size or standard deviation,
#' if method is "gauss"
#' @param method method for drawing the differential for the causal effects. can be "runif", "exp" or "gauss".
#' @author Martin Pirkl
#' @return graph with new edge weights
#' @export
#' @examples
#' graph.wt <- as(matrix(c(0,0,0,1,0,0,0,1,0), 3), "graphNEL")
#' graph.mt <- resample_edge_weights(graph.wt)
resample_edge_weights <- function(g, tp = 1,
                                  mineff = 1, maxeff = 2,
                                  method = "runif") {
    gold <- g
    g <- as(g, "matrix")
    changes <- floor((1-tp)*sum(g != 0))
    gbin <- g
    gbin[which(gbin != 0)] <- 1
    gtr <- gtr2 <- nem::transitive.reduction(gbin)
    diag(gtr) <- 1
    noedges <- 1
    start <- TRUE
    keep <- NULL
    while(changes > 0 & sum(gtr2) != 0 | start) {
        start <- FALSE
        if (changes == 1 & sum(gtr2) == 1) {
            keep <- c(keep, which(gtr2 == 1))
        } else {
            keep <- c(keep,
                      sample(which(gtr2 == 1),
                             min(changes, sum(gtr2))))
        }
        changes <- changes - sum(gtr2)
        diag(gtr2) <- 1
        gtr2[noedges] <- 1
        gtr2[keep] <- 1
        gtr2 <- gtr2%*%gtr
        gtr2[which(gtr2 > 0)] <- 1
        noedges <- which(gbin == 0 & gtr2 == 1)
        gtr2[noedges] <- 0
        gtr2[keep] <- 0
    }
    g2 <- g
    change <- which(g != 0)
    change <- unlist(lapply(change, function(x) {
        eff <- g2[x]
        die <- sample(c(0,1), 1)
        if (die) {
            if (method == "runif") {
                effnew <- runif(1, eff+mineff, eff+maxeff)
            } else if (method == "gauss") {
                effnew <- eff+mineff+abs(rnorm(1, 0, maxeff))
            } else if (method == "exp") {
                effnew <- eff+mineff+rexp(1, maxeff)
            }
        } else {
            if (method == "runif") {
                effnew <- runif(1, eff-maxeff, eff-mineff)
            } else if (method == "gauss") {
                effnew <- eff-mineff-abs(rnorm(1, 0, maxeff))
            } else if (method == "exp") {
                effnew <- eff-mineff-rexp(1, maxeff)
            }
        }
        return(effnew)
    }))
    g[which(g != 0)] <- change
    g[keep] <- g2[keep]
    edges <- which(gbin == 1, arr.ind = TRUE)
    for (i in seq(nrow(edges))) {
        edge <- paste0(rownames(gbin)[edges[i, 1]], "|", rownames(gbin)[edges[i, 2]])
        gold@edgeData@data[[edge]]$weight <- g[edges[i, 1],edges[i, 2]]
    }
    return(gold)
}
#' @noRd
#' @author modified code from the package 'bayesm'
#' @param y numeric vector of response
#' @param X numeric matrix with variables as columns
#' @param theta inverse dispersion parameter
#' @param link link function as character: "identity" or "log"
#' @param intercept logical to model intercept or not
#' @param family character of distribution: "nbinom"
glm.mle <- function (formula, data=NULL, theta=NULL, link="identity",
                     intercept=TRUE, family = "nbinom") {
    if (link %in% "identity") {
        linkfun <- meanfun <- function(x) return(x)
    } else if (link %in% "log") {
        linkfun <- function(x) return(log(x))
        meanfun <- function(x) return(exp(x))
    } else {
        linkfun <- function(x) return(log(x, link))
        meanfun <- function(x) return(link^x)
    }
    if (is(formula, "formula")) {
        D <- model.frame(formula)
        D <- as(D, "matrix")
        terms <- terms(formula)
        factors <- apply(attr(terms, "factors"), c(1,2), as.numeric)
        y <- D[, 1]
        X <- Xcn <- NULL
        for (i in seq_len(ncol(factors))) {
            vars <- rownames(factors)[which(factors[, i] != 0)]
            vars2 <- gsub("\\[", "\\\\[", vars)
            vars2 <- gsub("\\]", "\\\\]", vars2)
            if (sum(factors[, i]) > 1) {
                vars3 <- rownames(factors)[which(factors[, i] != 0)]
                vars3 <- gsub("\\[", "\\\\[", vars3)
                vars3 <- gsub("\\]", "\\\\]", vars3)
                tmp <- paste0(c(paste0("^", vars3[1], "$"),
                                paste0("^", vars3[1], "\\.")), collapse = "|")
                A <- D[, grep(tmp, colnames(D)), drop = FALSE]
                tmp <- paste0(c(paste0("^", vars3[2], "$"),
                                paste0("^", vars3[2], "\\.")), collapse = "|")
                B <- D[, grep(tmp, colnames(D)), drop = FALSE]
                Ap <- A[, rep(seq_len(ncol(A)), each = ncol(B)), drop = FALSE]
                X <- cbind(X, Ap*B)
                Xcn <- c(Xcn, paste0(colnames(A), ":", colnames(B)))
            } else {
                tmp <- paste0(c(paste0("^", vars2, "$"),
                                paste0("^", vars2, "\\.")), collapse = "|")
                X <- cbind(X, D[ , grep(tmp, colnames(D)), drop = FALSE])
                Xcn <- c(Xcn, colnames(D)[grep(tmp, colnames(D))])
            }
        }
        colnames(X) <- Xcn
    } else {
        y <- formula
        X <- data
    }
    if (is.numeric(intercept)) {
        int.fixed <- intercept
        intercept <- FALSE
    } else {
        int.fixed <- 0
    }
    llnegbin = function(par, X, y, nvar) {
        beta = par[1:nvar]
        if (intercept) {
            int <- par[nvar+1]
        } else {
            int <- 0
        }
        if (link %in% "identity") {
            int <- int + int.fixed
            mu <- meanfun(X %*% beta + int - min(X %*% beta) + mean(y))
        } else {
            mu <- meanfun(X %*% beta)
        }
        if (is.null(theta)) {
            alpha <- par[nvar+2]
        } else {
            alpha <- theta
        }
        out = dnbinom(y, size = alpha, mu = mu, log = TRUE)
        return(sum(out))
    }
    nvar = ncol(X)
    nobs = length(y)
    par = c(rep(1, nvar+1), 1)
    mle = optim(par, llnegbin, X = X, y = y, nvar = nvar,
                method = "BFGS",
                hessian = TRUE, control = list(fnscale = -1, maxit = 1000))
    if (!intercept) { mle$par[nvar+1] <- int.fixed }
    beta = mle$par[1:nvar]
    intercept <- mle$par[nvar+1]
    coefficients <- c(intercept, beta)
    names(coefficients) <- c("intercept", colnames(X))
    hessian <- mle$hessian[c(nvar+1, seq_len(nvar)), c(nvar+1, seq_len(nvar))]
    colnames(hessian) <- rownames(hessian) <- names(coefficients)
    if (is.null(theta)) {
        alpha <- mle$par[nvar+2]
    } else {
        alpha <- theta
    }
    result <- list(coefficients = coefficients, hessian = hessian,
                   theta=alpha)
    class(result) <- "glmmle"
    return(result)
}
#' Estimate p-values
#'
#' Function to compute p-values based on regression coefficients
#' @param x object of class "glmmle"
#' @author Martin Pirkl
#' @return data.frame with coefficients and p-values
#' @export
#' @importFrom aod wald.test
#' @method summary glmmle
summary.glmmle <- function(x) {
    coef <- x$coefficients
    hess <- x$hessian
    df <- seq_len(length(coef))
    hess <- hess[df, df]
    hessGlob <- hess
    coef <- coef[df]
    cov <- Gsolve(-hess)
    var.cf <- diag(cov)
    s.err <- sqrt(var.cf)
    tvalue <- coef/s.err
    ptvalue <- 2*pt(-abs(tvalue), ncol(hess))
    pzvalue <- 2*pnorm(-abs(tvalue))
    y <- list()
    y$coefficients <- cbind(coef, s.err, tvalue, ptvalue, pzvalue)
    colnames(y$coefficients) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "P(>|z|)")
    return(y)
}
#' @importFrom edgeR DGEList calcNormFactors estimateDisp
#' @noRd
estimateTheta <- function(data, ...) {
    if (is.null(dim(data))) {
        data <- matrix(data, length(data))
    }
    if (ncol(data) < 5) {
        mus <- apply(data, 2, mean)
        sigmas <- apply(data, 2, sd)
        thetas <- mus^2/(sigmas^2 - mus)
        theta <- mean(thetas)
    } else {
        y <- edgeR::DGEList(counts=t(data))
        y <- edgeR::calcNormFactors(y)
        y <- edgeR::estimateDisp(y, ...)
        theta <- 1/y$common.dispersion
    }
    return(theta)
}
#' @export
make.log.link <- function(base=exp(1)) {
    structure(list(
        linkfun=function(mu) { log(mu, base) },
        linkinv=function(eta) { base**(eta) },
        mu.eta=function(eta) { base**(eta) },
        valideta=function(eta) { TRUE }
    ), class="link-glm")
}
#' Identity robust glm fit.
#'
#' A modified glm.fit function which prevents negative mean
#' values in case of the identity link function. Convergence
#' is not guaranteed!
#' @export
glm.dce.fit <- function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
    mustart = NULL, offset = rep(0, nobs), family = gaussian(),
    control = list(), intercept = TRUE, singular.ok = TRUE)
{
    control <- do.call("glm.control", control)
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y))
        rownames(y)
    else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object",
            call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null) if (is.null(x))
        if.null
    else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta))
            stop("invalid linear predictor values in empty model",
                call. = FALSE)
        mu <- linkinv(eta)
        if (!validmu(mu))
            stop("invalid fitted means in empty model", call. = FALSE)
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric()
        iter <- 0L
    }
    else {
        coefold <- NULL
        eta <- if (!is.null(etastart))
            etastart
        else if (!is.null(start))
            if (length(start) != nvars)
                stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                  nvars, paste(deparse(xnames), collapse = ", ")),
                  domain = NA)
            else {
                coefold <- start
                offset + as.vector(if (NCOL(x) == 1L)
                  x * start
                else x %*% start)
            }
        else family$linkfun(mustart)
        offset0 <- 1
        mu <- linkinv(eta)
        if (any(mu <= 0)) {
            offset2 <- - min(eta) + offset0
            mu <- linkinv(eta <- eta + offset2)
        }
        if (!(validmu(mu) && valideta(eta)))
            stop("cannot find valid starting values: please specify some",
                call. = FALSE)
        devold <- sum(dev.resids(y, mu, weights))
        boundary <- conv <- FALSE
        for (iter in 1L:control$maxit) {
            good <- weights > 0
            varmu <- variance(mu)[good]
            if (any(is.na(varmu)))
                stop("NAs in V(mu)")
            if (any(varmu == 0))
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good])))
                stop("NAs in d(mu)/d(eta)")
            good <- (weights > 0) & (mu.eta.val != 0)
            if (all(!good)) {
                conv <- FALSE
                warning("no observations informative at iteration ",
                  iter)
                break
            }
            z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
            ngoodobs <- as.integer(nobs - sum(!good))
            fit <- lm.fit(x = x[good, , drop = FALSE] * w, y = z *
                          w, singular.ok = TRUE, tol = min(1e-07, control$epsilon/1000))
            fit$coefficients[is.na(fit$coefficients)] <- 0
            if (any(!is.finite(fit$coefficients))) {
                conv <- FALSE
                warning(gettextf("non-finite coefficients at iteration %d",
                  iter), domain = NA)
                break
            }
            if (nobs < fit$rank)
                stop(gettextf("X matrix has rank %d, but only %d observations",
                  fit$rank, nobs), domain = NA)
            start <- fit$coefficients
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
            if (any(mu <= 0)) {
                offset2 <- - min(eta) + offset0
                mu <- linkinv(eta <- eta + offset2)
            }
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace)
                cat("Deviance =", dev, "Iterations -", iter,
                  "\n")
            boundary <- FALSE
            if (!is.finite(dev)) {
                if (is.null(coefold)) {
                    stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
                }
                warning("step size truncated due to divergence", call. = FALSE)
                ii <- 1
                while (!is.finite(dev)) {
                  if (ii > control$maxit)
                    stop("inner loop 1; cannot correct step size",
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  if (any(mu <= 0)) {
                      offset2 <- - min(eta) + offset0
                      mu <- linkinv(eta <- eta + offset2)
                  }
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace)
                  cat("Step halved: new deviance =", dev, "\n")
            }
            if (!(valideta(eta) && validmu(mu))) {
                if (is.null(coefold)) {
                    stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
                }
                warning("step size truncated: out of bounds", call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxit)
                    stop("inner loop 2; cannot correct step size",
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  if (any(mu <= 0)) {
                      offset2 <- - min(eta) + offset0
                      mu <- linkinv(eta <- eta + offset2)
                  }
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace)
                  cat("Step halved: new deviance =", dev, "\n")
            }
            if (((dev - devold)/(0.1 + abs(dev)) >= control$epsilon) &
                (iter > 1)) {
                if (is.null(coefold)) {
                    stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
                }
                warning("step size truncated due to increasing deviance", call. = FALSE)
                ii <- 1
                while ((dev - devold)/(0.1 + abs(dev)) > -control$epsilon) {
                  if (ii > control$maxit)
                    break
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  if (any(mu <= 0)) {
                      offset2 <- - min(eta) + offset0
                      mu <- linkinv(eta <- eta + offset2)
                  }
                  dev <- sum(dev.resids(y, mu, weights))
                }
                if (ii > control$maxit)
                    warning("inner loop 3; cannot correct step size")
                else if (control$trace)
                    cat("Step halved: new deviance =", dev, "\n")
            }
            if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
                conv <- TRUE
                coef <- start
                break
            }
            else {
                devold <- dev
                coef <- coefold <- start
            }
        }
        if (!conv)
            warning("glm.dce.fit: algorithm did not converge. Try increasing the maximum iterations", call. = FALSE)
        if (boundary)
            warning("glm.dce.fit: algorithm stopped at boundary value",
                call. = FALSE)
        eps <- 10 * .Machine$double.eps
        if (family$family == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps))
                warning("glm.dce.fit: fitted probabilities numerically 0 or 1 occurred",
                  call. = FALSE)
        }
        if (family$family == "poisson") {
            if (any(mu < eps))
                warning("glm.dce.fit: fitted rates numerically 0 occurred",
                  call. = FALSE)
        }
        if (fit$rank < nvars) {
            if (!singular.ok)
                stop("singular fit encountered")
            coef[fit$qr$pivot][seq.int(fit$rank + 1, nvars)] <- NA
        }
        xxnames <- xnames[fit$qr$pivot]
        residuals <- (y - mu)/mu.eta(eta)
        fit$qr$qr <- as.matrix(fit$qr$qr)
        nr <- min(sum(good), nvars)
        if (nr < nvars) {
            Rmat <- diag(nvars)
            Rmat[1L:nr, 1L:nvars] <- fit$qr$qr[1L:nr, 1L:nvars]
        }
        else Rmat <- fit$qr$qr[1L:nvars, 1L:nvars]
        Rmat <- as.matrix(Rmat)
        Rmat[row(Rmat) > col(Rmat)] <- 0
        names(coef) <- xnames
        colnames(fit$qr$qr) <- xxnames
        dimnames(Rmat) <- list(xxnames, xxnames)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    if (!EMPTY)
        names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("",
            sum(good) - fit$rank))
    wtdmu <- if (intercept)
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY)
        0
    else fit$rank
    resdf <- n.ok - rank
    aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
    list(coefficients = coef, residuals = residuals, fitted.values = mu,
        effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat,
        rank = rank, qr = if (!EMPTY) structure(fit$qr[c("qr",
            "rank", "qraux", "pivot", "tol")], class = "qr"),
        family = family, linear.predictors = eta, deviance = dev,
        aic = aic.model, null.deviance = nulldev, iter = iter,
        weights = wt, prior.weights = weights, df.residual = resdf,
        df.null = nulldf, y = y, converged = conv, boundary = boundary)
}
#' Plot dce
#'
#' Wrapper function for mnem::plotDnf and epiNEM::HeatmapOP
#' @param x object of class dce
#' @param type "graph" or "heatmap"
#' @param log if TRUE, applies log to dces
#' @param scalefac scalar for scaling the dces
#' @param edgecol optional
#' @param edgewidth otpional
#' @param genelabels named vector with names like the original names
#' used in the dce object and values as subtitute names
#' @param sort Colv and Rowv have to be set to FALSE;
#' default sorting is by "name"; "strength" does sort for
#' effects strength and "topo" does a topological sorting,
#' alternatively a clusterx argument can be set to sort the rows/columns
#' by the clustering of clusterx
#' @param ... more arguments for mnem::plotDnf or epiNEM::HeatmapOP
#' @importFrom mnem plotDnf
#' @importFrom epiNEM HeatmapOP
#' @export
plotDce <- function(x, type = "graph", log = FALSE, scalefac = NULL,
                    edgecol = NULL, edgewidth = NULL, genelabels = NULL,
                    sort = "name", ...) {
    if (type == "graph") {
        graph <- x$graph
        graph <- as(x$graph, "matrix")
        graph1 <- which(graph == 1, arr.ind = TRUE)
        graph <- apply(graph1, 1, function(x) {
            y <- paste0(rownames(graph)[x[1]], "=", colnames(graph)[x[2]])
            return(y)
        })
        tmp <- x$dce
        tmp <- as.vector(tmp[which(tmp != 0)])
        if (log) {
            geq0 <- which(tmp > 0)
            leq0 <- which(tmp < 0)
            tmp[geq0] <- log(tmp[geq0] + 1)
            tmp[leq0] <- -log(-tmp[leq0] + 1)
        }
        if (is.null(scalefac)) {
            scalefac <- max(abs(tmp))
        }
        tmpBlue <- tmp/scalefac
        tmpRed <- -tmp/scalefac
        tmpBlue[which(tmpBlue < 0)] <- 0
        tmpRed[which(tmpRed < 0)] <- 0
        if (is.null(edgecol)) {
            if (max(tmpRed) > 1) {
                tmpRed2 <- tmpRed/max(tmpRed)
            } else {
                tmpRed2 <- tmpRed
            }
            if (max(tmpBlue) > 1) {
                tmpBlue2 <- tmpBlue/max(tmpBlue)
            } else {
                tmpBlue2 <- tmpBlue
            }
            if (log) {
                lowerbound <- 0
            } else {
                lowerbound <- 0.1
            }
            edgecol <- rgb(tmpRed2,0,tmpBlue2,apply(cbind(tmpBlue2+tmpRed2,
                                                          lowerbound), 1, max))
        }
        if (is.null(edgewidth)) {
            edgewidth <- apply(cbind(tmpBlue, tmpRed), 1, max)
        }
        mnem::plotDnf(graph, edgecol = edgecol,
                      edgewidth = (edgewidth*9)+1, ...)
    } else if (type == "heatmap") {
        graph <- as(x$graph, "matrix")
        tmp <- x$dce
        NAidx <- which(is.na(tmp) == TRUE)
        tmp[NAidx] <- 0
        if (log) {
            geq0 <- which(tmp > 0)
            leq0 <- which(tmp < 0)
            tmp[geq0] <- log(tmp[geq0] + 1)
            tmp[leq0] <- -log(-tmp[leq0] + 1)
        }
        tmp <- tmp/scalefac
        tmp[NAidx] <- NA
        tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]
        graph <- nem::transitive.closure(graph[order(rownames(graph)),
                                               order(colnames(graph))],
                                         mat = TRUE)
        if (sort == "topo") {
            topoord <- order(apply(graph, 1, sum) -
                             apply(graph, 2, sum), decreasing = TRUE)
            tmp <- tmp[topoord, topoord]
        } else if (sort == "strength") {
            strord <- order(apply(abs(tmp), 1, sum, na.rm = TRUE) -
                            apply(abs(tmp), 2, sum, na.rm = TRUE),
                            decreasing = TRUE)
            tmp <- tmp[strord, strord]
        }
        if (!is.null(genelabels)) {
            rownames(tmp) <-
                unlist(lapply(rownames(tmp), function(x) {
                    y <- genelabels[which(names(genelabels) == x)]
                }))
            colnames(tmp) <-
                unlist(lapply(colnames(tmp), function(x) {
                    y <- genelabels[which(names(genelabels) == x)]
                }))
        }
        if (sort == "names") {
            tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]
        }
        epiNEM::HeatmapOP(tmp, ...)
    }
}
