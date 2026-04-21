\name{create_random_W}
\alias{create_random_W}
\title{Generate a Random Connected Adjacency Matrix}
\description{
Constructs a symmetric binary adjacency matrix for a connected graph on \code{n_locs} nodes using a random spanning-tree process. This is the graph structure used in the simulation study.
}
\usage{
create_random_W(n_locs)
}
\arguments{
\item{n_locs}{Integer number of spatial units. Must be at least 2.}
}
\value{
An \code{n_locs} by \code{n_locs} symmetric matrix with entries in \eqn{\{0, 1\}} and zeros on the diagonal. The resulting graph is guaranteed to be connected with exactly \code{n_locs - 1} edges.
}
\details{
The algorithm draws a uniformly random permutation of the node indices and then, for each node after the first, adds an edge between that node and a randomly selected node drawn from those already placed. The result is a uniformly random labeled spanning tree on \code{n_locs} nodes.
}
\examples{
set.seed(1)
W <- create_random_W(25)
dim(W)
rowSums(W)
}
\seealso{
\code{\link{create_W_from_shapefile}} for generating \code{W} from a real spatial polygon dataset.
}
