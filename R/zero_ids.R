#' individuals without interactions
#'
#' @param mat a square matrix
#' @param return_ids logical, if \code{FALSE} (default) the number of
#'        individuals is returned, otherwise the row/column indices or
#'        column names if present
#'
#' @return numeric
#' @export
#'
#' @examples
#' zero_ids(bonobos)
#' zero_ids(bonobos, return_ids = TRUE)
#' m <- matrix(nrow = 4, ncol = 4, 0)
#' m[1, 3] <- 1
#' zero_ids(m)
#' zero_ids(m, return_ids = TRUE)


zero_ids <- function(mat, return_ids = FALSE) {
  diag(mat) <- 0
  xids <- which(colSums(mat) + rowSums(mat) == 0)

  if (return_ids) {
    if (length(xids) == 0) {
      return(NULL)
    }
    if (is.null(dimnames(mat))) {
      return(xids)
    } else {
      return(colnames(mat)[xids])
    }
  } else {
    return(length(xids))
  }
}
