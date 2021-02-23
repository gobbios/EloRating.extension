#' steepness via repeatability (cf aniDom package)
#'
#' @param xdata matrix or eloobject
#' @param reps number of randomized sequences
#' @param normprob logical, use exponential curve (default) or normal curve
#'        for calculating winning probabilities
#' @references Sanchez-Tojar et al 2018
#' @return a steepness value
#' @importFrom rptR rptGaussian
#' @importFrom EloChoice elochoice
#' @importFrom EloRating mat2seq
#' @export
#'
#' @examples
#' data(bonobos, package = "EloRating")
#' steep_anidom(bonobos, reps = 20)

steep_anidom <- function(xdata, reps = 1000, normprob = FALSE) {
  output <- NULL

  if (inherits(xdata, what = "matrix")) {
    s <- mat2seq(xdata)
    res <- elochoice(winner = s[, 1],
                     loser = s[, 2],
                     normprob = normprob,
                     runs = reps + 1)
    res <- res$ratmat[-1, ]
    xdata <- data.frame(id = rep(colnames(res), each = reps),
                        rating = as.numeric(res))
    r <- rptGaussian(rating ~ (1 | id),
                     "id",
                     data = xdata,
                     nboot = 0,
                     npermut = 0)
    output <- as.numeric(r$R)
  }

  output
}
