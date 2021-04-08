

#' initial values for Bayesian Elo steepness
#'
#' @param mat square interaction matrix
#' @param rem numeric, the proportion of interactions to be removed if any.
#'        Default is 0, i.e. all interactions are used.
#' @param chains the number of sets of initial values
#'
#' @return a list of length \code{chains}
#' @export
#' @importFrom EloRating fastelo
#'
#' @examples
#' data("bonobos", package = "EloRating")
#' initial_values(bonobos)
initial_values <- function(mat, rem = 0, chains = 4) {
  locs <- as.matrix(which(upper.tri(mat), arr.ind = TRUE))
  dimnames(mat)
  if (!is.null(unlist(dimnames(mat)))) {
    xn <- colnames(mat)
  } else {
    xn <- replicate(ncol(mat), paste(sample(c(letters, LETTERS), 10, replace = TRUE), collapse = ""))
    # colnames(mat) <- xn
    # rownames(mat) <- xn
  }


  w1 <- unlist(apply(locs, 1, function(x) rep(xn[x[1]], mat[x[1], x[2]])))
  l1 <- unlist(apply(locs, 1, function(x) rep(xn[x[2]], mat[x[1], x[2]])))
  w2 <- unlist(apply(locs, 1, function(x) rep(xn[x[2]], mat[x[2], x[1]])))
  l2 <- unlist(apply(locs, 1, function(x) rep(xn[x[1]], mat[x[2], x[1]])))
  xwinner <- c(w1, w2)
  xloser <- c(l1, l2)

  n <- length(xwinner)
  res <- vector("list", chains)
  for (i in seq_len(chains)) {
    # select and shuffle interactions
    xsel <- sample(seq_len(n), ceiling(n * (1 - rem)), replace = FALSE)
    allids <- unique(c(xwinner[xsel], xloser[xsel]))
    length(allids)

    # create ratings
    xres <- fastelo(WINNER = xwinner[xsel], LOSER = xloser[xsel], ALLIDS = allids, KVALS = rep(100, n), STARTVALUES = rep(0, length(allids)), NORMPROB = FALSE, ROUND = FALSE)$ratings

    # make sure that any individuals that were removed during the sampling receive a start value of 0
    restemplate <- rep(0, ncol(mat))
    names(restemplate) <- xn
    restemplate[names(xres)] <- xres

    xres <- round(as.numeric(scale(restemplate)), 4)
    if (!is.null(unlist(dimnames(mat)))) {
      names(xres) <- xn
    }
    res[[i]] <- list(EloStart_raw = xres)
    # names(res)[i] <- "EloStart"
  }

  res
}
