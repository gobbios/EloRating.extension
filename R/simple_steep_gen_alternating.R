#' generate dominance matrix with specified steepness
#'
#' @param n_ind integer, the number of individuals
#' @param n_int integer, the number of interactions
#' @param steeps numeric (between 0 and 1), the desired steepness values. The
#'        length of this vector determines how many periods will be generated.
#' @param id_bias numeric, between 0 and 1. If 0 all individual are equally
#'        likely to interact. If 1, some individuals have higher propensities
#'        to interact
#' @param rank_bias numeric, between 0 and 1. If 0 there is no relationship
#'        between rank distance and interaction propensity. If 1 there is a
#'        strong relationship: dyads closer in rank interact more often.
#' @return a list with the first item being the interactions in sequence form
#'         and the second item as matrix
#' @export
#'
#' @importFrom grDevices hcl.colors adjustcolor
#' @importFrom graphics par axis polygon points abline
#' @importFrom stats dnorm
#'
#' @details In essence, this is just a wrapper for
#'          \code{\link{simple_steep_gen}}.
#'
#' @examples
#' simple_steep_gen_alternating(n_ind = 3, n_int = 10, steeps = c(0.3, 0.9))


simple_steep_gen_alternating <- function(n_ind,
                                         n_int,
                                         steeps,
                                         id_bias = 0,
                                         rank_bias = 0) {

  dyad_weights <- generate_interaction_probs(n_ind = n_ind,
                                             id_bias = id_bias,
                                             rank_bias = rank_bias)

  # generate the individual periods
  res <- lapply(steeps,
                simple_steep_gen,
                n_ind = n_ind,
                n_int = n_int,
                id_bias = id_bias,
                rank_bias = rank_bias)

  # and extract matrices as list and sequence
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

  list(dyad_weights = dyad_weights,
       steeps = steeps,
       matrices = matrices,
       sequence = xsequence)
}
