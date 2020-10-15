#' dce - main function
#'
#' Main function to compute differential causal effects and its enrichment
#' @param graph valid object defining a directed acyclic graph
#' @param df_expr_wt data frame with wild type expression values
#' @param df_expr_mt data from with mutation type expression values
#' @param solver character with name of solver function
#' @param solver_args additional arguments for the solver function
#' @param adjustment_type character string for the method to define
#' the adjustment set Z for the regression
#' @param effect_type method of computing causal effects
#' @param p_method character string. "mean", "sum" for standard summary
#' functions, "hmp" for harmonic mean, "test" for the selfcontained
#' test of package 'CombinePValue' or any method from package 'metap',
#' e.g., "meanp" or "sump".
#' @param test either "wald" for testing significance with the
#' wald test or "lr" for using a likelihood ratio test
#' @param lib_size either a numeric vector of the same length as the
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
        graph, df_expr_wt, df_expr_mt,
        solver = "glm2", solver_args = list(method = "glm.dce.fit"),
        adjustment_type = "parents",
        effect_type = "total",
        p_method = "mean",
        test = "wald",
        lib_size = FALSE,
        latent = 0,
        conservative = FALSE,
        verbose = FALSE
    ) {
        standardGeneric("dce")
    },
    package = "dce"
)


# "igraph" is not a formal S4 class, make it compatible with `signature` call
setOldClass("igraph")
setMethod(
    "dce",
    signature = signature(graph = "igraph"),
    function(
        graph, df_expr_wt, df_expr_mt,
        solver = "glm2", solver_args = list(method = "glm.dce.fit"),
        adjustment_type = "parents",
        effect_type = "total",
        p_method = "mean",
        test = "wald",
        lib_size = FALSE,
        latent = 0,
        conservative = FALSE,
        verbose = FALSE
    ) {
        graph <- igraph::igraph.to.graphNEL(graph)
        dce(
            as_adjmat(graph),
            df_expr_wt, df_expr_mt,
            solver, solver_args,
            adjustment_type,
            effect_type,
            p_method,
            test,
            lib_size,
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
        graph, df_expr_wt, df_expr_mt,
        solver = "glm2", solver_args = list(method = "glm.dce.fit"),
        adjustment_type = "parents",
        effect_type = "total",
        p_method = "mean",
        test = "wald",
        lib_size = FALSE,
        latent = 0,
        conservative = FALSE,
        verbose = FALSE
    ) {
        dce(
            as_adjmat(graph),
            df_expr_wt, df_expr_mt,
            solver, solver_args,
            adjustment_type,
            effect_type,
            p_method,
            test,
            lib_size,
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
        graph, df_expr_wt, df_expr_mt,
        solver = "glm2", solver_args = list(method = "glm.dce.fit"),
        adjustment_type = "parents",
        effect_type = "total",
        p_method = "mean",
        test = "wald",
        lib_size = FALSE,
        latent = 0,
        conservative = FALSE,
        verbose = FALSE
    ) {
        # preparations
        graph[graph != 0] <- 1 # ignore edge weights

        # detect cycles in DAG
        graph <- topologically_ordering(graph)
        if (any(graph[lower.tri(graph)] == 1)) {
            warning("Cycle(s) detected in network")
        }

        # validate input
        if (!all(colnames(graph) %in% colnames(df_expr_wt))) {
            stop("Not all nodes have expression vector in WT data")
        }
        if (!all(colnames(graph) %in% colnames(df_expr_mt))) {
            stop("Not all nodes have expression vector in MT data")
        }

        # fit model
        .dce(
            graph, df_expr_wt, df_expr_mt,
            solver, solver_args,
            adjustment_type,
            effect_type,
            p_method,
            test,
            lib_size,
            latent,
            conservative,
            verbose
        )
    }
)


#' @export
dce_nb <- function(
    graph, df_expr_wt, df_expr_mt,
    solver_args = list(method = "glm.dce.nb.fit", link = "identity"),
    adjustment_type = "parents",
    effect_type = "total",
    p_method = "mean",
    test = "wald",
    lib_size = FALSE,
    latent = 0,
    conservative = FALSE,
    verbose = FALSE
) {
    dce(
        graph, df_expr_wt, df_expr_mt,
        solver = "glm.nb", solver_args = solver_args,
        adjustment_type,
        effect_type,
        p_method,
        test,
        lib_size,
        latent,
        conservative,
        verbose
    )
}


#' @importFrom naturalsort naturalorder
#' @noRd
.dce <- function(
    graph, df_expr_wt, df_expr_mt,
    solver, solver_args,
    adjustment_type,
    effect_type,
    p_method,
    test,
    lib_size,
    latent,
    conservative,
    verbose
) {
    # handle latent variables
    if (latent > 0) {
        pca_wt <- prcomp(scale(df_expr_wt))
        lat_wt <- pca_wt$x[, seq_len(latent), drop = FALSE]
        colnames(lat_wt) <- paste0("H", seq_len(latent))
        pca_mt <- prcomp(scale(df_expr_mt))
        lat_mt <- pca_mt$x[, seq_len(latent), drop = FALSE]
        colnames(lat.mt) <- paste0("H", seq_len(latent))
        lat_data <- rbind(lat_wt, lat_mt)
    }

    # handle lib_size
    if (!is.numeric(lib_size)) {
        if (lib_size) {
            lib_size <- apply(rbind(df_expr_wt, df_expr_mt), 1, sum)
            lib_size <- round(lib_size / (10^min(round(log10(lib_size)))))
        }
    }

    if (length(unique(lib_size)) == 1 & lib_size[1] != FALSE) {
        print("Only single library size detected, disabling correction!")
        lib_size <- FALSE
    }

    # subset expression data to pathway genes
    df_expr_wt <- df_expr_wt[, which(colnames(df_expr_wt) %in% colnames(graph))]
    df_expr_mt <- df_expr_mt[, which(colnames(df_expr_mt) %in% colnames(graph))]

    # ensure the data and graph have the same order of genes
    df_expr_wt <- df_expr_wt[, naturalorder(colnames(df_expr_wt))]
    df_expr_mt <- df_expr_mt[, naturalorder(colnames(df_expr_mt))]
    graph <- graph[naturalorder(rownames(graph)), naturalorder(colnames(graph))]

    # handle empty graph (no edges)
    if (sum(graph) == 0) {
        return(structure(list(
            graph = graph,
            dce = graph * NA,
            dce_pvalue = graph * NA,
            pathway_pvalue = NA
        ), class = "dce"))
    }

    # compute DCEs
    res <- purrr::pmap_dfr(
        as.data.frame(which(graph != 0, arr.ind = TRUE)),
        function(row, col) {
            if (row == col) {
                return(data.frame(
                    row = row,
                    col = col,
                    dce = NA,
                    stderr = NA,
                    p_value = NA
                ))
            }

            # concatenate data
            df_data <- data.frame(
                X = c(
                    df_expr_wt[, row],
                    df_expr_mt[, row]
                ),
                Y = c(
                    df_expr_wt[, col],
                    df_expr_mt[, col]
                ),
                N = c(
                    rep(0, dim(df_expr_wt)[[1]]),
                    rep(1, dim(df_expr_mt)[[1]])
                )
            )

            # incorporate adjustment set
            valid_adjustment_set <- get_adjustment_set(
                graph, row, col,
                adjustment_type,
                effect_type
            )

            form_adjustment_suffix <- ""
            for (idx in valid_adjustment_set) {
                name <- paste0("Z_", idx)
                if (conservative) {
                    add <- " + "
                } else {
                    add <- " + N * "
                }
                form_adjustment_suffix <- paste0(
                    form_adjustment_suffix,
                    add,
                    "`", name, "`"
                )
                df_data[, name] <- c(df_expr_wt[, idx], df_expr_mt[, idx])
            }

            # fit model
            if (!is.numeric(lib_size)) {
                form <- paste0("Y ~ N * X", form_adjustment_suffix)
            } else {
                df_data <- cbind(df_data, lib_size = factor(lib_size))
                form <- paste0("Y ~ N * X + N*lib_size", form_adjustment_suffix)
            }
            if (latent > 0) {
                df_data <- cbind(df_data, lat_data)
                form <- paste0(form, " + ",
                               paste(paste0("N*H", seq_len(latent)),
                                     collapse = " + "))
            }

            if (verbose) {
                print(df_data %>% head)
                print(form)
                cat("\n")
            }

            fit <- glm_solver(
                form = form, df = df_data,
                solver = solver, solver_args = solver_args
            )

            # extract results
            if (is.matrix(fit)) {
                # better support for custom solver functions
                coef_mat <- fit
            } else {
                coef_mat <- summary(fit)$coefficients
            }

            coef_xn <- NA
            stderr_xn <- NA
            pval_xn <- NA

            if (test == "lr") {
                if (!is.numeric(lib_size)) {
                    form2 <- paste0("Y ~ N + X", form_adjustment_suffix)
                } else {
                    form2 <- paste0("Y ~ N + X + N*lib_size", form_adjustment_suffix)
                }

                if (verbose) {
                    print(form2)
                }

                fit2 <- glm_solver(
                    form = form2, df = df_data,
                    solver = solver, solver_args = solver_args
                )

                if (length(grep("N:X", rownames(coef_mat))) != 0) {
                    coef_xn <- coef_mat["N:X", "Estimate"]
                    stderr_xn <- coef_mat["N:X", "Std. Error"]
                    pval_xn <- lmtest::lrtest(fit, fit2)[[5]][2]
                }
            } else if (test == "wald" & length(grep("N:X", rownames(coef_mat))) != 0) {
                coef_xn <- coef_mat["N:X", "Estimate"]
                stderr_xn <- coef_mat["N:X", "Std. Error"]
                pval_xn <- coef_mat["N:X", if (as.character(quote(solver)) == "glm.nb") "Pr(>|z|)" else "Pr(>|t|)"]
            } else if (test == "vcovHC" & length(grep("N:X", rownames(coef_mat))) != 0) {
                coef_xn <- coef_mat["N:X", "Estimate"]
                stderr_xn <- coef_mat["N:X", "Std. Error"]

                robust <- lmtest::coeftest(
                    fit, vcov = sandwich::vcovHC(fit, type = "HC0")
                )
                pval_xn <- robust["N:X", "Pr(>|t|)"]
            }

            data.frame(
                row = row,
                col = col,
                dce = coef_xn,
                stderr = stderr_xn,
                p_value = pval_xn
            )
        }
    )

    # process result
    dce_mat <- as.matrix(Matrix::sparseMatrix(
        res$row, res$col, x = res$dce, dims = dim(graph)
    ))
    rownames(dce_mat) <- colnames(dce_mat) <- rownames(graph)

    dce_stderr_mat <- as.matrix(Matrix::sparseMatrix(
        res$row, res$col, x = res$stderr, dims = dim(graph)
    ))
    rownames(dce_stderr_mat) <- colnames(dce_stderr_mat) <- rownames(graph)

    dce_pvalue_mat <- as.matrix(Matrix::sparseMatrix(
        res$row, res$col, x = res$p_value, dims = dim(graph)
    ))
    rownames(dce_pvalue_mat) <- colnames(dce_pvalue_mat) <- rownames(graph)

    # make uncomputed values NA
    diag(dce_mat) <- NA
    dce_mat[which(graph == 0)] <- NA

    diag(dce_stderr_mat) <- NA
    dce_stderr_mat[which(graph == 0)] <- NA

    diag(dce_pvalue_mat) <- NA
    dce_pvalue_mat[which(graph == 0)] <- NA

    # compute overall pathway enrichment
    tmp <- dce_pvalue_mat[!is.na(dce_pvalue_mat)]
    if ((length(tmp) > 0) && (sum(tmp != 0) == 0)) {
        # there are non-NA p-values, but all are zero
        pathway_pvalue <- 0
    } else {
        tmp[tmp == 0] <- min(tmp[tmp != 0])
        if (p_method == "hmp") {
            pathway_pvalue <- as.numeric(harmonicmeanp::p.hmp(tmp, L = length(tmp)))
        } else if (p_method == "test") {
            pathway_pvalue <- CombinePValue::selfcontained.test(tmp, weight = NA)[[1]]
        } else if (p_method %in% c("mean", "median", "sum", "max", "min")) {
            pathway_pvalue <- do.call(p_method, list(x = tmp))
        } else {
            require(metap)
            pathway_pvalue <- do.call(p_method, list(p = tmp))$p
        }
    }

    # return appropriate object
    structure(list(
        graph = graph,
        dce = dce_mat,
        dce_stderr = dce_stderr_mat,
        dce_pvalue = dce_pvalue_mat,
        pathway_pvalue = pathway_pvalue
    ), class = "dce")
}


#' @importFrom mnem transitive.closure
#' @export
get_adjustment_set <- function(
    graph, x, y,
    adjustment_type = "parents", effect_type = "total"
) {
    check_parents <- function(g, x, y) {
        parents <- which(g[, x] == 1)
        gxy <- g
        gxy[parents, x] <- 0
        gxy <- mnem::transitive.closure(gxy, mat = TRUE)
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
        adjustment_type,
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
    if (effect_type == "direct") {
        xkids <- names(which(graph[x, ] != 0))
        xkids <- xkids[which(xkids != colnames(graph)[y])]
        set <- c(set, xkids)
    }
    return(set)
}


#' @export
glm_solver <- function(form, df, solver, solver_args) {
    # handle general functions
    if (is.function(solver)) {
        return(solver(formula = form, data = df))
    }

    # lm solver
    if (solver == "lm") {
        return(lm(formula = form, data = df))
    }

    # glm solver
    solver_func <- switch(
        solver,
        "glm2" = glm2::glm2,
        "glm.nb" = glm.nb.rob,
        "mle" = glm.mle.new
    )
    if (is.null(solver_func)) {
        stop(paste("Invalid solver", solver))
    }

    func_args <- c(list(formula = form, data = df), solver_args)
    do.call(solver_func, func_args)
}


#' @export
#' @importFrom reshape2 melt
#' @importFrom dplyr rename mutate
as.data.frame.dce <- function(x) {
    x$dce %>%
        melt(.) %>%
        rename(dce = value, source = Var1, target = Var2) %>%
        mutate(dce_stderr = melt(x$dce_stderr)$value) %>%
        mutate(dce_pvalue = melt(x$dce_pvalue)$value)
}
