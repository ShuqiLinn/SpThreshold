\name{create_random_W}
\alias{create_random_W}
\title{Generate a Random Connected Adjacency Matrix}

\usage{
create_random_W(n_locs)
}

\arguments{
\item{n_locs}{Integer number of locations, at least 2.}
}

\description{
Randomly orders the locations and connects each location after the
first to a uniformly selected preceding location.
}

\value{
An \code{n_locs} by \code{n_locs} symmetric binary adjacency matrix
with a zero diagonal, representing a connected tree with exactly
\code{n_locs - 1} edges.
}