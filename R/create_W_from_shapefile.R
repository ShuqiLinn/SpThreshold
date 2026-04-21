#' Build an adjacency matrix from a spatial polygon object
#'
#' Constructs a binary adjacency matrix from an object of class \code{sf} or
#' \code{SpatialPolygonsDataFrame}, using shared borders or vertices to define
#' neighbors. Wraps \code{spdep::poly2nb()} and \code{spdep::nb2mat()}.
#'
#' @param polygons An object of class \code{sf}, \code{sfc}, or
#'   \code{SpatialPolygonsDataFrame}.
#' @param queen Logical. If \code{TRUE} (default), two areas are neighbors
#'   when they share at least one boundary point (queen contiguity). If
#'   \code{FALSE}, neighbors must share an edge of positive length (rook
#'   contiguity).
#' @param require_connected Logical. If \code{TRUE} (default), the function
#'   errors when the resulting graph is disconnected. Set to \code{FALSE} to
#'   return the disconnected adjacency matrix anyway.
#'
#' @return A symmetric \code{n_locs} by \code{n_locs} numeric matrix with
#'   entries in \{0, 1\} and zeros on the diagonal, where \code{n_locs} is the
#'   number of polygons.
#'
#' @details
#' Requires the \pkg{spdep} package. If \pkg{spdep} is not installed, the
#' function raises an informative error prompting the user to install it.
#'
#' Islands and disconnected components are common in real shapefiles and can
#' cause issues in CAR-based models. By default the function checks for
#' connectivity and errors if the graph is disconnected, so that users are
#' alerted early.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"))
#' W <- create_W_from_shapefile(nc)
#' dim(W)
#' }
#'
#' @seealso \code{\link{create_random_W}} for simulated adjacency structures.
#'
#' @export
create_W_from_shapefile <- function(polygons, queen = TRUE,
                                    require_connected = TRUE) {

   if (!requireNamespace("spdep", quietly = TRUE)) {
      stop("Package 'spdep' is required for create_W_from_shapefile(). ",
           "Install it with install.packages('spdep').")
   }

   nb <- spdep::poly2nb(polygons, queen = queen)

   ## Islands (regions with zero neighbors) will cause spdep::nb2mat to error
   ## under style = "B" unless zero.policy is set.
   n_islands <- sum(sapply(nb, function(x) length(x) == 1 && x == 0))

   W <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)

   ## Check connectivity if requested
   if (require_connected) {
      comps <- spdep::n.comp.nb(nb)
      if (comps$nc > 1) {
         stop(sprintf(
            "The adjacency graph has %d disconnected components (including %d isolated islands). ",
            comps$nc, n_islands),
            "Set require_connected = FALSE to return the matrix anyway, or resolve the ",
            "disconnection by, for example, linking isolated units to their nearest neighbor ",
            "using spdep::knn2nb().")
      }
   }

   ## Strip dimnames to return a plain numeric matrix
   dimnames(W) <- NULL
   W

}
