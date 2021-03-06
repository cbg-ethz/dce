% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_utils.R
\name{topologically_ordering}
\alias{topologically_ordering}
\title{Topological ordering}
\usage{
topologically_ordering(adja_mat, alt = FALSE)
}
\arguments{
\item{adja_mat}{Adjacency matrix of network}

\item{alt}{Use igraph implementation}
}
\value{
topologically ordered matrix
}
\description{
Order rows/columns of a adjacency matrix topologically
}
\examples{
adj <- matrix(c(0,1,0,0,0,1,0,0,0),3,3)
topologically_ordering(adj)
}
