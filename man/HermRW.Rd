% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/different_model.R
\name{HermRW}
\alias{HermRW}
\title{use herm model to cluster the graph}
\usage{
HermRW(k, cg1, D, n)
}
\arguments{
\item{k}{number of clusters}

\item{cg1}{complex matrix 1 generated from function DSBM}

\item{D}{diagonal matrix whose parameter is the degree of nodes}

\item{n}{nodes number}
}
\value{
kmeans result
}
\description{
use herm model to cluster the graph
}
\examples{
 G<- DSBM(4,100,0.8,0.4)
 HermRW(4,G$cg1,G$D,100)

}
