#' generate dominance matrix with specified steepness
#'
#' @param n_ind integer, the number of individuals
#' @param n_int integer, the number of interactions
#' @param steeps numeric (between 0 and 1), the desired steepness values. The
#'        length of this vector determines how many periods will be generated.
#' @param id_bias logical, should some individuals be more observed than others
#'        (default is \code{TRUE}). Can also be a numeric of the same length as
#'        \code{n_ind}.
#' @return a list with the first item being the interactions in sequence form
#'         and the second item as matrix
#' @export
#'
#' @examples
#' simple_steep_gen_alternating(n_ind = 3, n_int = 10, steeps = c(0.3, 0.9))


simple_steep_gen_alternating <- function(n_ind, n_int, steeps, id_bias = TRUE) {
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
  id_weights <- id_weights / sum(id_weights)

  # make dyads and give dyadic weight
  dyads <- t(combn(seq_len(n_ind), 2))
  dyads <- cbind(dyads, id_weights[dyads[, 1]] * id_weights[dyads[, 2]])

  # generate the individual periods
  res <- lapply(steeps,
                simple_steep_gen,
                n_ind = n_ind,
                n_int = n_int,
                id_bias = id_weights)

  matrices <- lapply(res, function(x)x$matrix)
  xsequence <- do.call("rbind", lapply(res, function(x)x$sequence))
  # generate dates for the sequence (1 year = 1 period = 1 element of steeps)
  dates <- c()
  for (i in seq_len(length(steeps))) {
    x <- seq.Date(from = as.Date(paste0(1970 + i, "-01-01")),
                  to = as.Date(paste0(1970 + i, "-12-31")),
                  by = "day")
    dates <- c(dates, sort(as.character(sample(x, n_int, replace = TRUE))))
  }
  xsequence <- data.frame(winner = xsequence[, 1],
                          loser = xsequence[, 2],
                          date = as.Date(dates))

  list(id_weights = id_weights,
       steeps = steeps,
       matrices = matrices,
       sequence = xsequence)
}
