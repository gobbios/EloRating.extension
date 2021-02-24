#' H steepness
#'
#' @param mat matrix
#' @param rand number of randomized matrices (default is 1000)
#'
#' @return a steepness value
#' @export
#'
#' @examples
#' data(bonobos)
#' steepness_h(bonobos)

steepness_h <- function(mat, rand = 1000) {
  res <- csteep(mat, rand)
  res$Hobs
}
