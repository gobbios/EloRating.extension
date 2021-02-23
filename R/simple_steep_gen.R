#' generate dominance matrix with specified steepness
#'
#' @param n_ind integer, the number of individuals
#' @param n_int integer, the number of interactions
#' @param steep numeric (between 0 and 1), the desired steepness value
#' @param id_bias logical, should some individuals be more observed than others
#'        (default is \code{TRUE}). Can also be a numeric of the same length as
#'        \code{n_ind} with actual weights.
#' @param rank_bias logical, should closer ranked dyads interact more often than
#'        dyads that are more distant in rank (default is \code{TRUE})
#' @importFrom stats runif pnorm
#' @importFrom utils combn
#' @return a list with the first item being the interactions in sequence form,
#'         the second item as matrix and the last a list with input settings
#' @export
#'
#' @examples
#' res <- simple_steep_gen(n_ind = 5, n_int = 30, steep = 0.99)
#' res$sequence
#' res$matrix
#'
#' # with prespecified weights
#' res <- simple_steep_gen(n_ind = 5,
#'                         n_int = 100,
#'                         steep = 0.95,
#'                         rank_bias = TRUE,
#'                         id_bias = c(1, 0.1, 0.1, 0.1, 0.1))
#' res$matrix
#'
#' library(EloRating)
#' steeps <- runif(20, 0, 1)
#' nids <- sample(6:10, length(steeps), TRUE)
#' mats <- sapply(1:length(steeps), function(x) {
#'   simple_steep_gen(nids[x], nids[x] ^ 2.5, steeps[x], FALSE)[[2]]
#'  })
#' obs_steeps <- unlist(lapply(mats, function(x)steepness(x)[1]))
#' plot(steeps, obs_steeps, xlim = c(0, 1), ylim = c(0, 1))
#' abline(0, 1)


simple_steep_gen <- function(n_ind,
                             n_int,
                             steep,
                             id_bias = TRUE,
                             rank_bias = TRUE) {
  # weights for ids
  if (is.logical(id_bias)) {
    id_weights <- runif(n = n_ind)
    if (!id_bias) {
      id_weights <- rep(1, n_ind)
    }
  } else {
    if (length(id_bias) != n_ind) {
      stop("id_bias must be the same length as n_ind", call. = FALSE)
    }
    id_weights <- id_bias
  }

  # make dyads and give dyadic weight
  dyads <- t(combn(seq_len(n_ind), 2))
  dyads <- cbind(dyads, id_weights[dyads[, 1]] * id_weights[dyads[, 2]])
  rank_weights <- rep(1, nrow(dyads))
  if (rank_bias) {
    rank_weights <- as.numeric(pnorm(scale(-abs(dyads[, 1] - dyads[, 2])),
                                     sd = 2))
    dyads[, 3] <- rank_weights
  }
  colnames(dyads) <- c("id1", "id2", "dyad_weights")

  # the number of dyads
  d <- seq_len(nrow(dyads))
  # output matrix
  m <- matrix(ncol = n_ind, nrow = n_ind, 0)
  # output sequence
  s <- matrix(ncol = 2, nrow = n_int)

  for (k in seq_len(n_int)) {
    x <- dyads[sample(d, size = 1, prob = dyads[, 3]), 1:2]
    if (runif(1, 0, 1) > (steep + 1) / 2) {
      x <- rev(x)
    }
    m[x[1], x[2]] <- m[x[1], x[2]] + 1
    s[k, ] <- x
  }

  colnames(m) <- paste0("i_", seq_len(n_ind))
  rownames(m) <- colnames(m)

  s[, 1] <- colnames(m)[s[, 1]]
  s[, 2] <- colnames(m)[as.numeric(s[, 2])]
  list(sequence = s,
       matrix = m,
       settings = list(dyads = dyads,
                       id_weights = id_weights,
                       rank_weights = rank_weights,
                       steep = steep))
}
