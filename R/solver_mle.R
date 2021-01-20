#' Modified theta.ml function from package MASS
#'
#' fixes a bug, if theta estimation breaks
#' see ?MASS::theta.ml for argument values
theta.ml.rob <- function(
    y, mu, n = sum(weights), weights,
    limit = 10, eps = .Machine$double.eps^0.25,
    trace = FALSE
) {
    score <- function(n, th, mu, y, w) sum(w * (digamma(th +
        y) - digamma(th) + log(th) + 1 - log(th + mu) - (y +
        th) / (mu + th)))
    info <- function(n, th, mu, y, w) sum(w * (-trigamma(th +
        y) + trigamma(th) - 1 / th + 2 / (mu + th) - (y + th) / (mu +
        th)^2))
    if (inherits(y, "lm")) {
        mu <- y$fitted.values
        y <- if (is.null(y$y))
            mu + residuals(y)
        else y$y
    }
    if (missing(weights))
        weights <- rep(1, length(y))
    t0 <- n / sum(weights * (y / mu - 1)^2)
    it <- 0
    del <- 1
    if (trace)
        message(sprintf("theta.ml.rob: iter %d 'theta = %f'", it,
            signif(t0)), domain = NA)
    while ((it <- it + 1) < limit && abs(del) > eps) {
        t0 <- abs(t0)
        del <- score(n, t0, mu, y, weights) / (i <- info(n, t0,
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
    attr(t0, "SE") <- sqrt(1 / i)
    t0
}


#' @noRd
#' @author modified code from the package 'bayesm'
#' @param y numeric vector of response
#' @param X numeric matrix with variables as columns
#' @param theta inverse dispersion parameter
#' @param link link function as character: "identity" or "log"
#' @param intercept logical to model intercept or not
#' @param family character of distribution: "nbinom"
glm.mle <- function(
    formula, data=NULL, theta=NULL, link="identity",
    intercept=TRUE, family = "nbinom"
) {
    if (link %in% "identity") {
        linkfun <- meanfun <- function(x) return(x)  # nolint
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
        factors <- apply(attr(terms, "factors"), c(1, 2), as.numeric)
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
                X <- cbind(X, Ap * B)
                Xcn <- c(Xcn, paste0(colnames(A), ":", colnames(B)))
            } else {
                tmp <- paste0(c(paste0("^", vars2, "$"),
                                paste0("^", vars2, "\\.")), collapse = "|")
                X <- cbind(X, D[, grep(tmp, colnames(D)), drop = FALSE])
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
    llnegbin <- function(par, X, y, nvar) {
        beta <- par[1:nvar]
        if (intercept) {
            int <- par[nvar + 1]
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
            alpha <- par[nvar + 2]
        } else {
            alpha <- theta
        }
        out <- dnbinom(y, size = alpha, mu = mu, log = TRUE)
        return(sum(out))
    }
    nvar <- ncol(X)
    nobs <- length(y)  # nolint
    par <- c(rep(1, nvar + 1), 1)
    mle <- optim(
        par, llnegbin, X = X, y = y, nvar = nvar,
        method = "BFGS", hessian = TRUE,
        control = list(fnscale = -1, maxit = 1000)
    )
    if (!intercept) {
        mle$par[nvar + 1] <- int.fixed
    }
    beta <- mle$par[1:nvar]
    intercept <- mle$par[nvar + 1]
    coefficients <- c(intercept, beta)
    names(coefficients) <- c("intercept", colnames(X))
    hessian <- mle$hessian[
        c(nvar + 1, seq_len(nvar)),
        c(nvar + 1, seq_len(nvar))
    ]
    colnames(hessian) <- rownames(hessian) <- names(coefficients)
    if (is.null(theta)) {
        alpha <- mle$par[nvar + 2]
    } else {
        alpha <- theta
    }
    result <- list(
        coefficients = coefficients,
        hessian = hessian,
        theta = alpha
    )
    class(result) <- "glmmle"
    return(result)
}


#' @export
glm.mle.new <- function(form, family, data, control=list()) {
    form <- as.formula(form)

    # massage data
    model.X <- model.matrix(form, data = data)

    df.model.frame <- model.frame(form, data = data)
    model.terms <- attributes(terms(df.model.frame))

    # solve (minimization)
     # terms, intercept and theta
    params <- rep(1, length(model.terms$term.labels) + 2)

    fit <- optim(
        params, loglikeli.func,
        X = model.X,
        Y = df.model.frame[, model.terms$response],
        family = family, method = "BFGS", hessian = TRUE,
        control = control
    )

    # return formatted result
    coef <- setNames(fit$par, c("Intercept", model.terms$term.labels, "Theta"))

    structure(
        list(
            coefficients = coef,
            log.likelihood = -fit$value,
            hessian = fit$hessian
        ),
        class = "glm.mle"
    )
}


#' @export
loglikeli.func <- function(params, X, Y, family) {
    beta.vec <- params[-length(params)]
    theta <- params[length(params)]

    mu <- family$linkinv(X %*% beta.vec)

    if (any(mu <= 0)) {
        return(NA)
    }

    -sum(dnbinom(Y, size = theta, mu = mu, log = TRUE))
}


#' @export
#' @method summary glm.mle
summary.glm.mle <- function(object, ...) {
    # get coefficients
    coef <- object$coefficients %>%
        as.data.frame %>%
        dplyr::rename(Estimate = ".")

    # compute p-values for $H_0: \beta_i = 0$
    cov.mat <- MASS::ginv(object$hessian)
    var.vec <- diag(cov.mat)
    sd.vec <- sqrt(var.vec)

    t.value <- object$coefficients / sd.vec
    pt.value <- 2 * pt(-abs(t.value), ncol(object$hessian))

    coef[, "Pr(>|t|)"] <- pt.value

    # prepare output
    list(
        coefficients = coef
    )
}
