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
        effect.type = "total",
        p.method = "mean",
        test = "wald",
        lib.size = FALSE,
        latent = 0,
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
        effect.type = "total",
        p.method = "mean",
        test = "wald",
        lib.size = FALSE,
        latent = 0,
        conservative = FALSE,
        verbose = FALSE
    ) {
        graph <- igraph::igraph.to.graphNEL(graph)
        dce(
            as.adjmat(graph),
            df.expr.wt, df.expr.mt,
            solver, solver.args,
            adjustment.type,
            effect.type,
            p.method,
            test,
            lib.size,
            latent,
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
        effect.type = "total",
        p.method = "mean",
        test = "wald",
        lib.size = FALSE,
        latent = 0,
        conservative = FALSE,
        verbose = FALSE
    ) {
        dce(
            as.adjmat(graph),
            df.expr.wt, df.expr.mt,
            solver, solver.args,
            adjustment.type,
            effect.type,
            p.method,
            test,
            lib.size,
            latent,
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
        effect.type = "total",
        p.method = "mean",
        test = "wald",
        lib.size = FALSE,
        latent = 0,
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
            effect.type,
            p.method,
            test,
            lib.size,
            latent,
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
    effect.type,
    p.method,
    test,
    lib.size,
    latent,
    conservative,
    verbose
) {
    # handle latent variables
    if (latent > 0) {
        pca.wt <- prcomp(df.expr.wt)
        lat.wt <- pca.wt$x[, seq_len(latent), drop = FALSE]
        colnames(lat.wt) <- paste0("H", seq_len(latent))
        pca.mt <- prcomp(df.expr.mt)
        lat.mt <- pca.mt$x[, seq_len(latent), drop = FALSE]
        colnames(lat.mt) <- paste0("H", seq_len(latent))
        lat.data <- rbind(lat.wt,lat.mt)
    }

    # handle lib.size
    if (!is.numeric(lib.size)) {
        if (lib.size) {
            lib.size <- apply(rbind(df.expr.wt, df.expr.mt), 1, sum)
            lib.size <- round(lib.size/(10^min(round(log10(lib.size)))))
        }
    }

    if (length(unique(lib.size)) == 1 & lib.size[1] != FALSE) {
        print("Only single library size detected, disabling correction!")
        lib.size <- FALSE
    }

    # subset expression data to pathway genes
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
            valid.adjustment.set <- get.adjustment.set(graph, row, col,
                                                       adjustment.type,
                                                       effect.type)

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
            if (latent > 0) {
                df.data <- cbind(df.data, lat.data)
                form <- paste0(form, " + ",
                               paste(paste0("N*H", seq_len(latent)),
                                     collapse = " + "))
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
            stderr.xn <- NA
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
                    stderr.xn <- coef.mat["N:X", "Std. Error"]
                    pval.xn <- lmtest::lrtest(fit, fit2)[[5]][2]
                }
            } else if (test == "wald" & length(grep("N:X", rownames(coef.mat))) != 0) {
                coef.xn <- coef.mat["N:X", "Estimate"]
                stderr.xn <- coef.mat["N:X", "Std. Error"]
                pval.xn <- coef.mat["N:X", if (solver == "glm.nb") "Pr(>|z|)" else "Pr(>|t|)"]
            }
            data.frame(
                row = row,
                col = col,
                dce = coef.xn,
                stderr = stderr.xn,
                p.value = pval.xn
            )
        }
    )

    # process result
    dce.mat <- as.matrix(Matrix::sparseMatrix(
        res$row, res$col, x = res$dce, dims = dim(graph)
    ))
    rownames(dce.mat) <- colnames(dce.mat) <- rownames(graph)

    dce.stderr.mat <- as.matrix(Matrix::sparseMatrix(
        res$row, res$col, x = res$stderr, dims = dim(graph)
    ))
    rownames(dce.stderr.mat) <- colnames(dce.stderr.mat) <- rownames(graph)

    dce.pvalue.mat <- as.matrix(Matrix::sparseMatrix(
        res$row, res$col, x = res$p.value, dims = dim(graph)
    ))
    rownames(dce.pvalue.mat) <- colnames(dce.pvalue.mat) <- rownames(graph)

    # make uncomputed values NA
    diag(dce.mat) <- NA
    dce.mat[which(graph == 0)] <- NA

    diag(dce.stderr.mat) <- NA
    dce.stderr.mat[which(graph == 0)] <- NA

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
        dce.stderr = dce.stderr.mat,
        dce.pvalue = dce.pvalue.mat,
        pathway.pvalue = pathway.pvalue
    ), class="dce")
}


#' @export
dce.nb <- function(
    graph, df.expr.wt, df.expr.mt,
    solver.args = list(method = "glm.dce.nb.fit", link = "identity"),
    adjustment.type = "parents",
    effect.type = "total",
    p.method = "mean",
    test = "wald",
    lib.size = FALSE,
    latent = 0,
    conservative = FALSE,
    verbose = FALSE
) {
    dce(
        graph, df.expr.wt, df.expr.mt,
        solver = "glm.nb", solver.args = solver.args,
        adjustment.type,
        effect.type,
        p.method,
        test,
        lib.size,
        latent,
        conservative,
        verbose
    )
}


#' @export
get.adjustment.set <- function(graph, x, y, adjustment.type = "parents",
                               effect.type = "total") {
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
            set <- names(which(graph[, x] != 0))
        },
        minimal = {
            tmp <- pcalg::adjustment(graph, "dag", x, y, "minimal")

            if (length(tmp) > 0) {
                set <- rownames(graph)[tmp[[1]]]
            } else {
                set <- vector(mode = "character")
            }
        },
        parents_filtered = {
            set <- rownames(graph)[check_parents(graph, x, y)]
        }
    )
    if (effect.type == "direct") {
        xkids <- names(which(graph[x, ] != 0))
        xkids <- xkids[which(xkids != colnames(graph)[y])]
        set <- c(set, xkids)
    }
    return(set)
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
#' @importFrom reshape2 melt
as.data.frame.dce <- function(x) {
    x$dce %>%
        melt(.) %>%
        rename(dce = value, source = Var1, target = Var2) %>%
        mutate(dce.stderr = melt(x$dce.stderr)$value) %>%
        mutate(dce.pvalue = melt(x$dce.pvalue)$value)
}
