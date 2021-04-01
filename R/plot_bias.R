#' illustrate rank-distance bias in interaction matrix
#'
#' @param x results of \code{\link{simple_steep_gen}}
#'
#' @return a plot
#' @export
#' @importFrom grDevices hcl.colors
#' @importFrom graphics rect text
#'
#' @examples
#' x <- simple_steep_gen(n_ind = 20, n_int = 1000, steep = 0.8, id_bias = 0, rank_bias = 1)
#' plot_bias(x)
#' interaction_bias(x$matrix)
#' x <- simple_steep_gen(n_ind = 20, n_int = 5000, steep = 0.8, id_bias = 0, rank_bias = 0)
#' plot_bias(x)
#' interaction_bias(x$matrix)

plot_bias <- function(x) {
  mo <- x$matrix
  n <- ncol(mo)

  m <- mo + t(mo)
  m[lower.tri(m)] <- 0

  xrange <- range(m[upper.tri(m)])
  splits <- seq(from = xrange[1] - 0.0001, to = xrange[2], length.out = 21)

  p <- cut(m, breaks = splits)
  pmat <- matrix(as.numeric(p), ncol = n)

  # colours
  xcols <- hcl.colors(n = 20, palette = "viridis", alpha = 0.6, rev = TRUE)

  plot(0, 0, type = "n", xlim = c(0, ncol(m) + 1), ylim = c(ncol(m) + 1, 0), axes = FALSE, ann = FALSE)
  i=1
  j=2
  for (i in 1:(nrow(m) - 1)) {
    for (j in (i + 1):nrow(m)) {
      rect(xleft = i - 0.5, ybottom = j - 0.5, xright = i + 0.5, ytop = j + 0.5, border = NA, col = xcols[pmat[i, j]])
      text(x = i, y = j, mo[j, i], col = "white", cex = 0.7)
      rect(xleft = j - 0.5, ybottom = i - 0.5, xright = j + 0.5, ytop = i + 0.5, border = NA, col = xcols[pmat[i, j]])
      text(x = j, y = i, mo[i, j], col = "white", cex = 0.7)
    }
  }
}
