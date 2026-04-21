#' Generate a random connected adjacency matrix
#'
#' Constructs a symmetric binary adjacency matrix for a connected graph on
#' \code{n_locs} nodes using a random spanning-tree process. This is the graph
#' structure used in the simulation study of the manuscript.
#'
#' @param n_locs Integer number of spatial units. Must be at least 2.
#'
#' @return An \code{n_locs} by \code{n_locs} symmetric matrix with entries in
#'   \{0, 1\}, and zeros on the diagonal. The resulting graph is guaranteed to
#'   be connected with exactly \code{n_locs - 1} edges.
#'
#' @details
#' The algorithm draws a uniformly random permutation of the node indices and
#' then, for each node after the first, adds an edge between that node and a
#' randomly selected node drawn from those already placed. The result is a
#' uniformly random labeled spanning tree on \code{n_locs} nodes.
#'
#' @examples
#' set.seed(1)
#' W <- create_random_W(25)
#' dim(W)
#' rowSums(W)   ## degree of each node
#'
#' @seealso \code{\link{create_W_from_shapefile}} for generating \code{W} from
#'   a real spatial polygon dataset.
#'
#'
#' @export
create_random_W <- function(n_locs) {

   if (!is.numeric(n_locs) || length(n_locs) != 1 || n_locs < 2 ||
       n_locs != as.integer(n_locs)) {
      stop("n_locs must be a single integer >= 2.")
   }
   n_locs <- as.integer(n_locs)

   W <- matrix(0, nrow = n_locs, ncol = n_locs)
   perm <- sample(n_locs)

   for (i in 2:n_locs) {
      a <- perm[i]
      b <- perm[sample.int(i - 1, 1)]
      W[a, b] <- 1
      W[b, a] <- 1
   }

   W

}
