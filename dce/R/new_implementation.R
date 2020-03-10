#' dce - main function
#'
#' Main function to compute differential causal effects and its enrichment
#' @param graph valid object defining a directed acyclic graph
#' @param df.expr.wt data frame with wild type expression values
#' @param df.expr.mt data from with mutation type expression values
#' @param solver character with name of solver function
#' @param solver.args additional arguments for the solver function
#' @param adjustment.type character string for the method to define
#' the adjustment set Z for the regression
#' @param p.method character string. "hmp" for harmonic mean
#' or any method from package 'metap', e.g., "meanp" or "sump"
#' @param verbose logical for verbose output
#' @export
#' @importFrom graph graphNEL
#' @importFrom igraph as_adjacency_matrix
#' @importFrom Matrix sparseMatrix
#' @import metap
setGeneric(
    "dce",
    function(
        graph, df.expr.wt, df.expr.mt,
        solver = "glm2", solver.args = list(method = "glm.dce.fit"),
        adjustment.type = "parents",
        p.method = "meanp",
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
        solver = "glm2", solver.args = list(method = "glm.dce.fit"),
        adjustment.type = "parents",
        p.method = "meanp",
        verbose = FALSE
    ) {
        graph <- igraph::igraph.to.graphNEL(graph)
        dce(
            as.adjmat(graph),
            df.expr.wt, df.expr.mt,
            solver, solver.args,
            adjustment.type,
            p.method,
            verbose
        )
    }
)


setMethod(
    "dce",
    signature = signature(graph = "graphNEL"),
    function(
        graph, df.expr.wt, df.expr.mt,
        solver = "glm2", solver.args = list(method = "glm.dce.fit"),
        adjustment.type = "parents",
        p.method = "meanp",
        verbose = FALSE
    ) {
        dce(
            as.adjmat(graph),
            df.expr.wt, df.expr.mt,
            solver, solver.args,
            adjustment.type,
            p.method,
            verbose
        )
    }
)


setMethod(
    "dce",
    signature = signature(graph = "matrix"),
    function(
        graph, df.expr.wt, df.expr.mt,
        solver = "glm2", solver.args = list(method = "glm.dce.fit"),
        adjustment.type = "parents",
        p.method = "meanp",
        verbose = FALSE
    ) {
        # preparations
        graph[graph != 0] <- 1 # ignore edge weights

        # validate input
        if (!all(colnames(graph) %in% colnames(df.expr.wt))) {
            stop("Not all nodes have expression vector in WT data")
        }
        if (!all(colnames(graph) %in% colnames(df.expr.mt))) {
            stop("Not all nodes have expression vector in MT data")
        }

        # fit model
        .dce(
            graph, df.expr.wt, df.expr.mt,
            solver, solver.args,
            adjustment.type,
            p.method,
            verbose
        )
    }
)


#' @noRd
.dce <- function(
    graph, df.expr.wt, df.expr.mt,
    solver, solver.args,
    adjustment.type,
    p.method,
    verbose
) {
    # handle empty graph (no edges)
    if (sum(graph) == 0) {
        return(structure(list(
            graph = graph,
            dce = graph * NA,
            dce.pvalue = graph * NA,
            pathway.pvalue = NA
        ), class="dce"))
    }

    # compute DCEs
    res <- purrr::pmap_dfr(
        as.data.frame(which(graph != 0, arr.ind = TRUE)),
        function (row, col) {
            if (row == col) {
                return(data.frame(
                    row = row,
                    col = col,
                    dce = NA,
                    p.value = NA
                ))
            }

            # concatenate data
            df.data <- data.frame(
                X = c(df.expr.wt[, row], df.expr.mt[, row]),
                Y = c(df.expr.wt[, col], df.expr.mt[, col]),
                N = c(rep(0, dim(df.expr.wt)[[1]]), rep(1, dim(df.expr.mt)[[1]]))
            )

            # incorporate adjustment set
            valid.adjustment.set <- get.adjustment.set(graph, row, col, adjustment.type)

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
                solver = solver, solver.args = solver.args
            )

            # extract results
            coef.mat <- summary(fit)$coefficients

            data.frame(
                row = row,
                col = col,
                dce = coef.mat["N:X", "Estimate"],
                p.value = coef.mat["N:X", if (solver == "glm.nb") "Pr(>|z|)" else "Pr(>|t|)"]
            )
        }
    )

    # process result
    dce.mat <- as.matrix(Matrix::sparseMatrix(
        res$row, res$col, x = res$dce, dims = dim(graph)
    ))
    rownames(dce.mat) <- colnames(dce.mat) <- rownames(graph)

    dce.pvalue.mat <- as.matrix(Matrix::sparseMatrix(
        res$row, res$col, x = res$p.value, dims = dim(graph)
    ))
    rownames(dce.pvalue.mat) <- colnames(dce.pvalue.mat) <- rownames(graph)

    # make uncomputed values NA
    graph.tc <- nem::transitive.closure(graph, mat = TRUE)

    diag(dce.mat) <- NA
    dce.mat[which(graph.tc == 0)] <- NA

    diag(dce.pvalue.mat) <- NA
    dce.pvalue.mat[which(graph.tc == 0)] <- NA

    # compute overall pathway enrichment
    tmp <- dce.pvalue.mat[!is.na(dce.pvalue.mat)]
    tmp[tmp == 0] <- min(tmp[tmp != 0])
    if (p.method == "hmp") {
        pathway.pvalue <- as.numeric(harmonicmeanp::p.hmp(tmp, L = length(tmp)))
    } else {
        require(metap)
        pathway.pvalue <- do.call(p.method, list(p=tmp))$p
    }
    
    # return appropriate object
    structure(list(
        graph = graph,
        dce = dce.mat,
        dce.pvalue = dce.pvalue.mat,
        pathway.pvalue = pathway.pvalue
    ), class="dce")
}


#' @export
dce.nb <- function(
    graph, df.expr.wt, df.expr.mt,
    solver.args = list(method = "glm.dce.nb.fit", link = "identity"),
    adjustment.type = "parents",
    p.method = "meanp",
    verbose = FALSE
) {
    dce(
        graph, df.expr.wt, df.expr.mt,
        solver = "glm.nb", solver.args = solver.args,
        adjustment.type,
        p.method,
        verbose
    )
}

#' @export
get.adjustment.set <- function(graph, x, y, adjustment.type) {
    check_parents <- function(g, x, y) {
        parents <- which(g[, x] == 1)
        gxy <- g
        gxy[parents, x] <- 0
        gxy <- nem::transitive.closure(gxy, mat = TRUE)
        grandparents <- unlist(lapply(parents, function(u) {
            v <- which(gxy[, u] == 1)
            w <- any(gxy[v, y] == 1)
            return(w)
        }))
        set <- which(gxy[parents, y] == 1 | grandparents == TRUE)
        if (length(set) > 0) {
            minset <- parents[set]
        } else {
            minset <- NULL
        }
        return(minset)
    }
    switch(
        adjustment.type,
        parents = {
            names(which(graph[, x] != 0))
        },
        minimal = {
            tmp <- pcalg::adjustment(graph, "dag", x, y, "minimal")

            if (length(tmp) > 0) {
                rownames(graph)[tmp[[1]]]
            } else {
                vector(mode = "character")
            }
        },
        parents_filtered = {
            rownames(graph)[check_parents(graph, x, y)]
        }
    )
}


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


#' @export
glm.solver <- function(form, df, solver, solver.args) {
    solver.func <- switch(
        solver,
        "glm2" = glm2::glm2,
        "glm.nb" = glm.nb,
        "mle" = glm.mle.new
    )
    if (is.null(solver.func)) {
        stop(paste("Invalid solver", solver))
    }

    func.args <- c(list(formula = form, data = df), solver.args)
    do.call(solver.func, func.args)
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
#' @method summary glm.mle
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
