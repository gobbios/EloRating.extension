#' winning probabilities with individual performance spread
#'
#' @param r1,r2 numeric, the ratings of two individuals
#' @param s1,s2 numeric, the standard deviation of the individual performance spread. Must be positive. Its default, 200, corresponds to the standard.
#' @param do_plot logical, if \code{TRUE} produces a density plot of the two performance ratings and a plot of the distributions of the rating difference
#'
#' @return a winning probability for the first individual and optionally a plot
#' @export
#'
#' @examples
#' winprob_alt(1200, 1000)
#' winprob(1200, 1000)
#'
#' winprob_alt(1200, 1400)
#' winprob(1200, 1400)
#'
#' winprob_alt(1200, 1400, 100, 250)
#' winprob_alt(1200, 1400, 50, 250)

winprob_alt <- function(r1, r2, s1 = 200, s2 = 200, do_plot = TRUE) {
  # rating difference
  rdiff <- r1 - r2
  # SD of the difference based on the individual performance spread
  sdx <- (s1 ^ 2 + s2 ^ 2)^0.5
  # the winning probability
  res <- pnorm(rdiff / sdx)
  cols <- sample(hcl.colors(5, "Zissou1"), 2)
  if (do_plot) {
    par(mfrow = c(1, 2), family = "serif", mar = c(3, 2, 1, 1))
    spread <- 4
    xl <- min(r1 - spread * s1, r2 - spread * s2)
    xr <- max(r1 + spread * s1, r2 + spread * s2)
    x1 <- seq(r1 - spread * s1, r1 + spread * s1, length.out = 201)
    y1 <- dnorm(x1, r1, s1)
    x2 <- seq(r2 - spread * s2, r2 + spread * s2, length.out = 201)
    y2 <- dnorm(x2, r2, s2)
    plot(0, 0, "n", xlim = c(xl, xr), ylim = c(0, max(y1, y2) * 1.05), yaxs = "i", axes = FALSE, xlab = "individual ratings", ylab = "", mgp = c(1.5, 0, 0))
    axis(1, tcl = -0.2, mgp = c(1.5, 0.5, 0))

    polygon(c(x1, rev(x1)), c(rep(0, length(y1)), rev(y1)), border = NA, col = adjustcolor(cols[1], 0.8))
    points(x1, y1, type = "l")
    polygon(c(x2, rev(x2)), c(rep(0, length(y2)), rev(y2)), border = NA, col = adjustcolor(cols[2], 0.8))
    points(x2, y2, type = "l")

    x3 <- seq(rdiff - spread * sdx, rdiff + spread * sdx, length.out = 201)
    y3 <- dnorm(x3, rdiff, sdx)
    plot(0, 0, "n", xlim = range(x3), ylim = c(0, max(y3) * 1.05), yaxs = "i", mgp = c(1.5, 0, 0), ylab = "", xlab = "rating difference", axes = FALSE)
    axis(1, tcl = -0.2, mgp = c(1.5, 0.5, 0))
    abline(h = 0)

    x4 <- seq(0, rdiff + spread * sdx, length.out = 201)
    y4 <- dnorm(x4, rdiff, sdx)
    polygon(c(x4, rev(x4)), c(y4, rep(0, length(y4))), border = NA, col = cols[1])
    x5 <- seq(0, rdiff - spread * sdx, length.out = 201)
    y5 <- dnorm(x5, rdiff, sdx)
    polygon(c(x5, rev(x5)), c(y5, rep(0, length(y5))), border = NA, col = cols[2])
    abline(h = 0)

    points(x3, y3, type = "l")
  }

  res
}
