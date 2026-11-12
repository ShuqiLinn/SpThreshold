#' Generate a random connected adjacency matrix
#'
#' Randomly orders the locations and connects each location after the first
#' to a randomly selected preceding location.
#'
#' @param n_locs Integer number of locations, at least 2.
#' @return A symmetric binary adjacency matrix with a zero diagonal,
#'   representing a connected tree with \code{n_locs - 1} edges.
#' @export
create_random_W <- function(n_locs) {
   if (!is.numeric(n_locs) || length(n_locs) != 1 || n_locs < 2 ||
       n_locs != as.integer(n_locs))
      stop("n_locs must be a single integer >= 2.")
   
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