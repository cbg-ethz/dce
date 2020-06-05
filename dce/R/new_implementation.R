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
#' @param p.method character string. "mean", "sum" for standard summary
#' functions, "hmp" for harmonic mean, "test" for the selfcontained
#' test of package 'CombinePValue' or any method from package 'metap',
#' e.g., "meanp" or "sump".
#' @param test either "wald" for testing significance with the
#' wald test or "lr" for using a likelihood ratio test
#' @param lib.size either a numeric vector of the same length as the
#' sum of wild type and mutant samples or a logical. If TRUE, it is
#' recommended that both data sets include not only the genes
#' included in the graph but all genes available in the original data set.
#' @param verbose logical for verbose output
#' @export
#' @importFrom graph graphNEL
#' @importFrom igraph as_adjacency_matrix
#' @importFrom Matrix sparseMatrix
#' @import metap CombinePValue
setGeneric(
    "dce",
    function(
        graph, df.expr.wt, df.expr.mt,
        solver = "glm2", solver.args = list(method = "glm.dce.fit"),
        adjustment.type = "parents",
        p.method = "mean",
        test = "wald",
        lib.size = FALSE,
        conservative = FALSE,
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
        p.method = "mean",
        test = "wald",
        lib.size = FALSE,
        conservative = FALSE,
        verbose = FALSE
    ) {
        graph <- igraph::igraph.to.graphNEL(graph)
        dce(
            as.adjmat(graph),
            df.expr.wt, df.expr.mt,
            solver, solver.args,
            adjustment.type,
            p.method,
            test,
            lib.size,
            conservative,
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
        p.method = "mean",
        test = "wald",
        lib.size = FALSE,
        conservative = FALSE,
        verbose = FALSE
    ) {
        dce(
            as.adjmat(graph),
            df.expr.wt, df.expr.mt,
            solver, solver.args,
            adjustment.type,
            p.method,
            test,
            lib.size,
            conservative,
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
        p.method = "mean",
        test = "wald",
        lib.size = FALSE,
        conservative = FALSE,
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
            test,
            lib.size,
            conservative,
            verbose
        )
    }
)


#' @importFrom naturalsort naturalorder
#' @noRd
.dce <- function(
    graph, df.expr.wt, df.expr.mt,
    solver, solver.args,
    adjustment.type,
    p.method,
    test,
    lib.size,
    conservative,
    verbose
) {
    # handle lib.size + filter data
    if (!is.numeric(lib.size)) {
        if (lib.size) {
            lib.size <- apply(rbind(df.expr.wt, df.expr.mt), 1, sum)
            lib.size <- round(lib.size/(10^min(round(log10(lib.size)))))
        }
    }
    df.expr.wt <- df.expr.wt[, which(colnames(df.expr.wt) %in% colnames(graph))]
    df.expr.mt <- df.expr.mt[, which(colnames(df.expr.mt) %in% colnames(graph))]

    # ensure the data and graph have the same order of genes
    df.expr.wt <- df.expr.wt[, naturalorder(colnames(df.expr.wt))]
    df.expr.mt <- df.expr.mt[, naturalorder(colnames(df.expr.mt))]
    graph <- graph[naturalorder(rownames(graph)), naturalorder(colnames(graph))]

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
                if (conservative) {
                    add <- " + "
                } else {
                    add <- " + N * "
                }
                form.adjustment.suffix <- paste0(
                    form.adjustment.suffix,
                    add,
                    name
                )
                df.data[, name] <- c(df.expr.wt[, idx], df.expr.mt[, idx])
            }

            # fit model
            if (!is.numeric(lib.size)) {
                form <- paste0("Y ~ N * X", form.adjustment.suffix)
            } else {
                df.data <- cbind(df.data, lib.size = factor(lib.size))
                form <- paste0("Y ~ N * X + N*lib.size", form.adjustment.suffix)
            }

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
            coef.xn <- NA
            pval.xn <- NA

            if (test == "lr") {
                if (!is.numeric(lib.size)) {
                    form2 <- paste0("Y ~ N + X", form.adjustment.suffix)
                } else {
                    form2 <- paste0("Y ~ N + X + N*lib.size", form.adjustment.suffix)
                }

                if (verbose) {
                    print(form2)
                }

                fit2 <- glm.solver(
                    form = form2, df = df.data,
                    solver = solver, solver.args = solver.args
                )

                if (length(grep("N:X", rownames(coef.mat))) != 0) {
                    coef.xn <- coef.mat["N:X", "Estimate"]
                    pval.xn <- lmtest::lrtest(fit, fit2)[[5]][2]
                }
            } else if (test == "wald" & length(grep("N:X", rownames(coef.mat))) != 0) {
                coef.xn <- coef.mat["N:X", "Estimate"]
                pval.xn <- coef.mat["N:X", if (solver == "glm.nb") "Pr(>|z|)" else "Pr(>|t|)"]
            }
            data.frame(
                row = row,
                col = col,
                dce = coef.xn,
                p.value = pval.xn
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
    diag(dce.mat) <- NA
    dce.mat[which(graph == 0)] <- NA

    diag(dce.pvalue.mat) <- NA
    dce.pvalue.mat[which(graph == 0)] <- NA

    # compute overall pathway enrichment
    tmp <- dce.pvalue.mat[!is.na(dce.pvalue.mat)]
    tmp[tmp == 0] <- min(tmp[tmp != 0])
    if (p.method == "hmp") {
        pathway.pvalue <- as.numeric(harmonicmeanp::p.hmp(tmp, L = length(tmp)))
    } else if (p.method == "test") {
        pathway.pvalue <- CombinePValue::selfcontained.test(tmp, weight=NA)[[1]]
    } else if (p.method %in% c("mean", "median", "sum", "max", "min")) {
        pathway.pvalue <- do.call(p.method, list(x=tmp))
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
    p.method = "mean",
    test = "wald",
    lib.size = FALSE,
    conservative = FALSE,
    verbose = FALSE
) {
    dce(
        graph, df.expr.wt, df.expr.mt,
        solver = "glm.nb", solver.args = solver.args,
        adjustment.type,
        p.method,
        test,
        lib.size,
        conservative,
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
        "glm.nb" = glm.nb.rob,
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


#' Plot dce object
#'
#' This function takes a differential causal effects object and plots
#' the dag with the dces
#' @param x dce object
#' @param nodename.map node names
#' @param edge.colorscale.limits Limits for scale_edge_color_gradient2 (should contain 0). Useful to make plot comparable to others
#' @param nodesize Node sizes
#' @param labelsize Node label sizes
#' @param show.edge.labels Whether to show edge labels (DCEs)
#' @param use.symlog Scale edge colors using dce::symlog
#' @param ... additional parameters
#' @author Martin Pirkl, Kim Philipp Jablonski
#' @method plot dce
#' @return plot of dag and dces
#' @export
#' @import tidyverse ggraph purrr
#' @importFrom ggplot2 aes theme element_rect arrow unit
#' @importFrom tidygraph as_tbl_graph activate mutate
#' @importFrom rlang .data
plot.dce <- function(
    x,
    nodename.map = NULL, edgescale.limits = NULL,
    nodesize = 17, labelsize = 3,
    show.edge.labels = FALSE, use.symlog = FALSE,
    ...
) {
    coords.dot <- purrr::map_dfr(
        Rgraphviz::agopen(as(x$graph, "graphNEL"), name="foo", layoutType="dot")@AgNode,
        function(node) {
            data.frame(x=node@center@x, y=node@center@y)
        }
    )

    if(is.null(edgescale.limits)) {
        edgescale.limits <- c(-max(abs(x$dce), na.rm=TRUE), max(abs(x$dce), na.rm=TRUE))
    }

    as_tbl_graph(as(x$graph, "graphNEL")) %>%
        activate(nodes) %>%
        mutate(
            label=if(is.null(nodename.map)) .data$name else nodename.map[.data$name],
            nodesize=nodesize
        ) %>%
        activate(edges) %>%
        mutate(
            dce=pmap_dbl(
                list(.data$from, .data$to),
                function (f, t) { x$dce[f, t] }
            ),
            dce.symlog=symlog(dce),
            label=.data$dce %>% round(2) %>% as.character
        ) %>%
    ggraph(layout=coords.dot) + # "sugiyama"
        geom_edge_diagonal(
            aes(
                color=if(use.symlog) .data$dce.symlog else .data$dce,
                alpha=abs(.data$dce),
                width=abs(.data$dce),
                label=if(show.edge.labels) .data$label else NULL,
                # start_cap = circle(.data$node1.nodesize, unit="native"), # uncomment once https://github.com/thomasp85/ggraph/pull/246 has been merged
                end_cap = circle(.data$node2.nodesize, unit="native")
            ),
            strength=0.5,
            arrow=arrow(type="closed", length=unit(3, "mm"))
        ) +
        geom_node_circle(aes(r=.data$nodesize), color="black", fill="white") +
        geom_node_text(aes(label=.data$label), size=labelsize) +
        coord_fixed() +
        scale_edge_color_gradient2(
            low="red", mid="grey", high="blue",
            midpoint=0, limits=if(use.symlog) symlog(edgescale.limits) else edgescale.limits,
            name=if(use.symlog) "DCE (symlog)" else "DCE"
        ) +
        scale_edge_width(range=c(1, 3), limits=c(0, edgescale.limits[[2]]), guide=FALSE) +
        scale_edge_alpha(range=c(.1, 1), limits=c(0, edgescale.limits[[2]]), guide=FALSE) +
        theme(
            panel.background=element_rect(fill="white"),
            # legend.position="none"
        )
}


#' @export
symlog <- function(x, base = 10, threshold = 1) {
    if(is.null(x)) {
        return(NULL)
    }

    ifelse(
        abs(x) < threshold,
        x,
        sign(x) * (threshold + log(sign(x) * x / threshold, base))
    )
}
