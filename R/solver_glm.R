#' @export
negative.binomial.special <- function(
    theta = NULL,
    linkinv = function(mu, offset = 1) { mu - min(mu) + offset }
) {
    funcs <- list(
        family = "negative.binomial.special",
        link = "special",

        linkfun = linkinv, # well...
        linkinv = linkinv
    )

    if (!is.null(theta)) {
        tmp <- MASS::negative.binomial(theta, link = "identity")

        funcs <- c(
            funcs,
            mu.eta = tmp$mu.eta,
            variance = tmp$variance,
            dev.resids = tmp$dev.resids,
            aic = tmp$aic
        )
    }

    structure(funcs, class = c("dce.nb.family", "family"))
}


#' Extra robust model fitting
#'
#' `glm.nb` uses two families: poisson (initial theta fit) and negative binomial (final fit).
#' The solver can fail if mu becomes negative
#' @export
glm.dce.nb.fit <- function(...) {
    args <- list(...)

    # if (startsWith(args$family$family, "Negative Binomial")) {
    #   theta <- as.numeric(str_match(args$family$family, "\\((.*?)\\)")[[2]]) # this is horrible
    #   args$family <- dce::negative.binomial.special(theta = theta)
    # }
    args$family$linkinv <- function(eta, offset = 1) {
        mu <- eta # identity

        if (any(mu <= 0)) {
            mu <- mu - min(mu) + offset
        }

        mu
    }

    # do.call(glm.fit, args)
    do.call(glm2::glm.fit2, args)
}


#' Modified glm.nb function from package MASS
#'
#' fixes a bug, if theta estimation breaks
#' see ?MASS::glm.nb for argument values
#' @export
glm.nb.rob <- function (
    formula, data, weights, subset, na.action, start = NULL,
    etastart, mustart, control = glm.control(...), method = "glm.fit",
    model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, ...,
    init.theta, link = log
)
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
