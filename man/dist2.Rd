% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{dist2}
\alias{dist2}
\title{This is the distance function extracted from SNFtool package.
Computes the Euclidean distances between all pairs of data point given.}
\usage{
dist2(X, C)
}
\arguments{
\item{X}{A data matrix where each row is a different data point}

\item{C}{A data matrix where each row is a different data point. If this matrix is the same as X, pairwise distances for all data points are computed.}
}
\value{
Returns an N x M matrix where N is the number of rows in X and M is the number of rows in M. element (n,m) is the squared Euclidean distance between nth data point in X and mth data point in C
}
\description{
This is the distance function extracted from SNFtool package.
Computes the Euclidean distances between all pairs of data point given.
}
\examples{
data(miRNAExp)
Dist1=dist2(t(miRNAExp),t(miRNAExp))

}
