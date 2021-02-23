#' steepness based on winning probabilities from Elo-rating
#'
#' @param x either a matrix or the results of \code{\link{elo.seq}}
#' @param date \code{Date} or character (YYYY-MM-DD), optional date to extract
#'        steepness on when \code{x} is the result of \code{\link{elo.seq}}.
#'        Default is \code{NULL}, i.e. the last date in the sequence. Ignored
#'        if \code{x} is a matrix.
#' @param sliding numeric, optional number of days centered around \code{date}
#'        over which steepness is averaged. Default is \code{0}, only the
#'        ratings on the specified date are used.
#' @param n_rand numeric, the number of randomized sequences generated if
#'        \code{x} is a matrix. Otherwise ignored. Default is \code{1000}.
#' @param return_mean logical, should the mean of the randomizations be
#'        returned or the individual values from randomizations or sliding
#' @param fixed_ranks logical, should the same ranking be used for the
#'        regression (as derived by David's score)
#'
#' @return an (average) steepness value ranging between 0 and 1.
#'
#' @importFrom EloRating DS fastelo winprob
#' @importFrom stats coef lm
#'
#' @export
#'
#' @examples
#' library(EloRating)
#' data(bonobos)
#' # from a matrix with ranomized sequences
#' steepness_elo(bonobos)
#' steepness_elo(bonobos, fixed_ranks = TRUE)
#'
#' # from sequence
#' data(adv)
#' e <- elo.seq(adv$winner, adv$loser, adv$Date)
#' steepness_elo(e)
#' steepness_elo(e, date = "2010-01-12")
#' steepness_elo(e, date = "2010-01-12", sliding = 3)

steepness_elo <- function(x,
                          date = NULL,
                          sliding = 0,
                          n_rand = 1000,
                          return_mean = TRUE,
                          fixed_ranks = FALSE) {
  if (inherits(x = x, what = "matrix")) {
    s <- mat2seq(x)
    nids <- nrow(x)
    nint <- nrow(s)
    steeps <- numeric(n_rand)
    xnames <- colnames(x)
    # get DS for ranking
    ds <- DS(x)
    rownames(ds) <- xnames
    dsranks <- rank(ds[xnames, "DS"])

    for (i in seq_len(n_rand)) {
      s <- s[sample(sample(seq_len(nrow(s)))), ]
      r <- fastelo(WINNER = s[, 1],
                   LOSER = s[, 2],
                   ALLIDS = rownames(x),
                   KVALS = rep(100, nint),
                   STARTVALUES = rep(0, nids),
                   NORMPROB = TRUE,
                   ROUND = TRUE)
      r <- r$ratings[xnames]
      wp <- outer(r, r, winprob)
      diag(wp) <- 0
      wp <- rowSums(wp)

      if (fixed_ranks) {
        steeps[i] <- coef(lm(wp ~ dsranks))[2]
      } else {
        ranks <- rank(wp)
        steeps[i] <- coef(lm(wp ~ ranks))[2]
      }

    }
  }


  if (inherits(x = x, what = "elo")) {
    eloobject <- x
    if (is.null(date)) {
      date <- rev(eloobject$truedates)[1]
    }
    xline <- which(eloobject$truedates == date)
    if (sliding > 0) {
      xline <- (xline - sliding) : (xline + sliding)
      xline <- xline[xline > 0 & xline < length(eloobject$truedates)]
    }
    steeps <- numeric(length(xline))
    for (i in seq_len(length(xline))) {
      r <- eloobject$cmat[xline[i], ]
      wp <- outer(r, r, winprob)
      diag(wp) <- 0
      wp <- rowSums(wp)
      ranks <- rank(wp)
      steeps[i] <- as.numeric(coef(lm(wp ~ ranks))[2])
    }
  }

  if (return_mean) {
    return(mean(steeps))
  } else {
    return(steeps)
  }
}
