\name{create_W_from_shapefile}
\alias{create_W_from_shapefile}
\title{Build an Adjacency Matrix from a Spatial Polygon Object}
\description{
Constructs a binary adjacency matrix from an object of class \code{sf} or \code{SpatialPolygonsDataFrame}, using shared borders or vertices to define neighbors. Wraps \code{spdep::poly2nb()} and \code{spdep::nb2mat()}.
}
\usage{
create_W_from_shapefile(polygons, queen = TRUE, require_connected = TRUE)
}
\arguments{
\item{polygons}{An object of class \code{sf}, \code{sfc}, or \code{SpatialPolygonsDataFrame}.}

\item{queen}{Logical. If \code{TRUE} (default), two areas are neighbors when they share at least one boundary point (queen contiguity). If \code{FALSE}, neighbors must share an edge of positive length (rook contiguity).}

\item{require_connected}{Logical. If \code{TRUE} (default), the function errors when the resulting graph is disconnected. Set to \code{FALSE} to return the disconnected adjacency matrix anyway.}
}
\value{
A symmetric \code{n_locs} by \code{n_locs} numeric matrix with entries in \eqn{\{0, 1\}} and zeros on the diagonal, where \code{n_locs} is the number of polygons.
}
\details{
Requires the \pkg{spdep} package. If \pkg{spdep} is not installed, the function raises an informative error prompting the user to install it.

Islands and disconnected components are common in real shapefiles and can cause issues in CAR-based models. By default the function checks for connectivity and errors if the graph is disconnected, so that users are alerted early.
}
\examples{
\dontrun{
library(sf)
nc <- st_read(system.file("shape/nc.shp", package = "sf"))
W <- create_W_from_shapefile(nc)
dim(W)
}
}
\seealso{
\code{\link{create_random_W}} for simulated adjacency structures.
}
