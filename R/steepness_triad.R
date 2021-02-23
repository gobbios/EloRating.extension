# prop="Pij"
# x=bonobos
# steepness(x)
# x <- matrix(ncol = 5, nrow = 5, 0)
# x[upper.tri(x)] <- 3
# # x[1, 2] <- 0
# steepness(x, Dij = FALSE)

steepness_triad <- function(x, prop = "Dij", mode = 1) {
  if (is.null(dimnames(x))) {
    colnames(x) <- paste0("i", seq_len(ncol(x)))
    rownames(x) <- colnames(x)
  }

  total_triads <- choose(ncol(x), 3)
  xtri <- t(combn(x = nrow(x), m = 3))
  max_tri_per_id <- sum(xtri[, 1] == 1)

  m <- matrix(ncol = 3, nrow = 3)
  tri <- combn(x = seq_len(nrow(x)),
               m = 3,
               FUN = function(y) x[y, y],
               simplify = FALSE)
  dslist <- lapply(tri, DS, prop = prop)
  tri_known <- unlist(lapply(tri, function(x) all((x + t(x))[upper.tri(m)] > 0)))
  # xtri <- cbind(xtri, tri_known)

  res <- data.frame(id = colnames(x), cumscore = NA, relscore = NA, ntriads = NA)
  for (i in seq_len(nrow(x))) {
    # which triads include i
    loc <- which(apply(xtri, 1, function(x) i %in% x ) & tri_known)
    if (length(loc) > 0) {
      # assemble scores for those triads
      s <- do.call("rbind", dslist[loc])
      # filter for i
      s <- s[s$ID == colnames(x)[i], ]
      res$cumscore[i] <- sum(s$normDS)
      res$relscore[i] <- sum(s$normDS)/nrow(s)
      res$ntriads[i] <- nrow(s)
    } else {
      res$ntriads[i] <- 0
    }

  }

  if (mode == 1) {
    ind_scores <- (res$cumscore / ( res$ntriads/ max_tri_per_id) ) / (max_tri_per_id / 2)
  } else {
    ind_scores <- (res$cumscore / res$ntriads) * ( max_tri_per_id / 3)
  }

  coef(lm(ind_scores ~ rank(ind_scores)))[2]
}

# DS(bonobos)
# steepness_triad(bonobos)
# x <- matrix(ncol = 5, nrow = 5, 0)
# x[upper.tri(x)] <- 3
# rownames(x) <- colnames(x) <- letters[1:5]
# DS(x)
# steepness(x)
# steepness_triad(x)
#
# x[1, 2] <- 0
# DS(x)
# steepness(x)
# steepness_triad(x)
# x[1, 3] <- 0
# DS(x)
# steepness(x)
# steepness_triad(x)
#
#
# x <- remove_dyads(bonobos)
# unlist(lapply(x$matrices, steepness_triad))
# for (i in length(x$matrices)) {
#   steepness_triad(x$matrices[[i]])
# }



x <- simple_steep_gen(8, 300, 0.8)
prunk(x$matrix)
steepness(x$matrix)
steepness_triad(x$matrix)
steepness(x$matrix, Dij = FALSE)
steepness_triad(x$matrix, prop = "Pij")
steepness_triad(x$matrix, prop = "Pij", mode =0)


y <- remove_dyads(x$matrix)
unlist(lapply(y$matrices, steepness_triad, prop = "Dij"))
