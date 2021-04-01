
#' interaction bias in interaction matrix
#'
#' correlation between rank/score distance and interaction numbers and
#' coefficient of variation between interactions of individuals
#'
#' @param m input matrix (square)
#' @param ranking_by character, the ranking method to be used (so far, only
#'        \code{"DS"} is supported)
#' @param n_rand numeric, number of randomizations (if \code{n_rand = 0}, the
#'        the default, no randomizations are performed)
#' @param concise_results logical, if \code{TRUE} (default), no randomziation
#'        estimates are returned. If \code{FALSE}, the results of the individual
#'        randomizations are returned in the results
#' @importFrom EloRating DS
#' @importFrom stats cor sd
#' @return a list
#' @export
#'
#' @examples
#' data(bonobos)
#' interaction_bias(bonobos, "DS", 100)
#'
#' res <- interaction_bias(bonobos, "DS", 1000, FALSE)
#' hist(res$randos[, "ordinal_cor"])
#' abline(v = res$output$value[1], col = "red")
interaction_bias <- function(m,
                             ranking_by = "DS",
                             n_rand = 0,
                             concise_results = TRUE) {
  # do the ranking
  if (ranking_by == "DS") {
    r <- DS(m)
    m <- m[as.character(r$ID), as.character(r$ID)]
  }
  # closer ranked interact more?
  tm <- m + t(m)
  out <- t(combn(seq_len(nrow(m)), 2))
  out <- cbind(out, apply(out, 1, function(x)tm[x[1], x[2]]))
  out <- cbind(out, out[, 2] - out[, 1])
  out1 <- cor(out[, 3], out[, 4], method = "s")

  # variation in individual number of interactions
  x <- rowSums(m) + colSums(m)
  out2 <- sd(x) / mean(x)

  # if cardinal ranking, cardinal differences are possible
  out3 <- NA
  if (ranking_by %in% c("DS")) {
    foo <- function(x) r$normDS[x[1]] - r$normDS[x[2]]
    out <- cbind(out, apply(out, 1, foo))
    out3 <- cor(out[, 4], out[, 5], method = "s")
  }

  # output
  output <- data.frame(measure = c("ordinal_cor", "CV", "cardinal_cor"),
                       value = c(out1, out2, out3),
                       p_right = NA)
  if (n_rand > 0) {
    p1 <- numeric(n_rand)
    p1[1] <- out1
    p2 <- numeric(n_rand)
    p2[1] <- out2
    for (i in 2:n_rand) {
      out[, 3] <- sample(out[, 3])
      p1[i] <- cor(out[, 3], out[, 4], method = "s")
      foo <- function(x) sum(out[out[, 1] == x | out[, 2] == x, 3])
      x <- sapply(seq_len(nrow(m)), foo)
      p2[i] <- sd(x) / mean(x)
    }


    p3 <- NA
    if (ranking_by %in% c("DS")) {
      p3 <- numeric(n_rand)
      p3[1] <- out3
      for (i in 2:n_rand) {
        out[, 3] <- sample(out[, 3])
        p3[i] <- cor(out[, 3], out[, 5], method = "s")
      }
    }
  }

  if (n_rand > 0) {
    output$p_right[1] <- sum(p1 >= out1) / n_rand
    output$p_right[2] <- sum(p2 >= out2) / n_rand
    if (!is.na(out3)) {
      output$p_right[3] <- sum(p3 >= out3) / n_rand
    }
  }

  if (concise_results) {
    output <- list(output = output,
                   randos = NULL)
  } else {
    output <- list(output = output,
                   randos = cbind(p1, p2, p3))
    colnames(output$randos) <- c("ordinal_cor", "CV", "cardinal_cor")
  }

  output
}
