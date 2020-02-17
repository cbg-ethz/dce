#' @export
#' @importFrom graph graphNEL
#' @importFrom igraph igraph
#' @importFrom Matrix sparseMatrix
setGeneric(
    "dce",
    function(
        graph, df.expr.wt, df.expr.mt,
        family, solver = "glm2",
        verbose = FALSE
    ) {
        standardGeneric("dce")
    },
    package = "dce"
)


setOldClass("igraph") # "igraph" is not a formal S4 class, make it compatible with `signature` call
setMethod(
    "dce",
    signature = signature(graph = "igraph"),
    function(
        graph, df.expr.wt, df.expr.mt,
        family, solver = "glm2",
        verbose = FALSE
    ) {
        dce(
            as(igraph::as_adjacency_matrix(graph, attr=NULL), "matrix"),
            df.expr.wt, df.expr.mt,
            family, solver,
            verbose
        )
    }
)


setMethod(
    "dce",
    signature = signature(graph = "graphNEL"),
    function(
        graph, df.expr.wt, df.expr.mt,
        family, solver = "glm2",
        verbose = FALSE
    ) {
        dce(
            as(graph, "matrix"),
            df.expr.wt, df.expr.mt,
            family, solver,
            verbose
        )
    }
)


setMethod(
    "dce",
    signature = signature(graph = "matrix"),
    function(
        graph, df.expr.wt, df.expr.mt,
        family, solver = "glm2",
        verbose = FALSE
    ) {
        # preparations
        graph[graph != 0] <- 1 # ignore edge weights

        # validate input
        nonzero.idx <- which(graph != 0, arr.ind = TRUE)
        if (
            (nrow(nonzero.idx) > 0) &&
            (
                any(nonzero.idx[, 2] - nonzero.idx[, 1] < 0) ||
                any(diag(graph) != 0)
            )
        ) {
            stop("Input DAG must be topologically ordered!")
        }

        # fit model
        .dce(
            graph, df.expr.wt, df.expr.mt,
            family, solver,
            verbose
        )
    }
)


#' @noRd
.dce <- function(
    graph, df.expr.wt, df.expr.mt,
    family, solver,
    verbose
) {
    # compute DCEs
    res <- purrr::pmap_dfr(
        as.data.frame(which(graph != 0, arr.ind = TRUE)),
        function (row, col) {
            if (row == col) {
                stop("oh no")
            }

            # concatenate data
            df.data <- data.frame(
                X = c(df.expr.wt[, row], df.expr.mt[, row]),
                Y = c(df.expr.wt[, col], df.expr.mt[, col]),
                N = c(rep(0, dim(df.expr.wt)[[1]]), rep(1, dim(df.expr.mt)[[1]]))
            )

            # incorporate adjustment set
            valid.adjustment.set <- which(graph[, row] != 0) # parents

            form.adjustment.suffix <- ""
            for (idx in valid.adjustment.set) {
                name <- paste0("Z", idx)
                form.adjustment.suffix <- paste0(
                    form.adjustment.suffix,
                    " + N * ",
                    name
                )
                df.data[, name] <- c(df.expr.wt[, idx], df.expr.mt[, idx])
            }

            # fit model
            form <- paste0("Y ~ N * X", form.adjustment.suffix)

            if (verbose) {
                print(df.data %>% head)
                print(form)
            }

            fit <- glm.solver(
                form = form, df = df.data,
                family = family, solver = solver
            )

            # extract results
            coef.mat <- summary(fit)$coefficients

            data.frame(
                row = row,
                col = col,
                dce = coef.mat["N:X", "Estimate"],
                p.value = coef.mat["N:X", "Pr(>|t|)"]
            )
        }
    )

    # return result (TODO: make not computed p-values NA)
    structure(list(
        graph = graph,
        dce = as.matrix(Matrix::sparseMatrix(res$row, res$col, x = res$dce, dims = dim(graph))),
        dce.pvalue = as.matrix(Matrix::sparseMatrix(res$row, res$col, x = res$p.value, dims = dim(graph)))
    ), class="dce")
}


#' @export
glm.solver <- function(form, df, family, solver) {
    if (solver == "glm2") {
        fit <- glm2::glm2(form, family = family, data = df)
    } else if (solver == "mle") {
        fit <- glm.mle.new(form, family = family, data = df)
    } else {
        stop(paste("Invalid solver", solver))
    }

    fit
}


#' @export
loglikeli.func <- function(params, X, Y, family) {
    beta.vec <- params[-length(params)]
    theta <- params[length(params)]

    mu <- family$linkinv(X %*% beta.vec)

    if (any(mu <= 0)) {
        return(NA)
    }

    -sum(dnbinom(Y, size=theta, mu=mu, log=TRUE))
}
#' @export
glm.mle.new <- function(form, family, data, control=list()) {
    form <- as.formula(form)

    # massage data
    model.X <- model.matrix(form, data=data)

    df.model.frame <- model.frame(form, data=data)
    model.terms <- attributes(terms(df.model.frame))

    # solve (minimization)
    params <- rep(1, length(model.terms$term.labels) + 2) # terms, intercept and theta
    fit <- optim(
        params, loglikeli.func,
        X=model.X, Y=df.model.frame[, model.terms$response], family = family,
        method="BFGS", hessian=TRUE,
        control=control
    )

    # return formatted result
    coef <- setNames(fit$par, c("Intercept", model.terms$term.labels, "Theta"))

    structure(
        list(
            coefficients = coef,
            log.likelihood = -fit$value,
            hessian = fit$hessian
        ),
        class="glm.mle"
    )
}
#' @export
#' @method summary glmmle
summary.glm.mle <- function(x) {
    # get coefficients
    coef <- x$coefficients %>%
        as.data.frame %>%
        dplyr::rename(Estimate=".")

    # compute p-values for $H_0: \beta_i = 0$
    cov.mat <- matlib::Ginv(x$hessian)
    var.vec <- diag(cov.mat)
    sd.vec <- sqrt(var.vec)

    t.value <- x$coefficients / sd.vec
    pt.value <- 2 * pt(-abs(t.value), ncol(x$hessian))

    coef[, "Pr(>|t|)"] <- pt.value

    # prepare output
    list(
        coefficients = coef
    )
}
