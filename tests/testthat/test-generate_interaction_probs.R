# test biases:
ntests <- 100
res <- data.frame(run = 1:ntests, idbias = runif(ntests), rankbias = runif(ntests), cv_id = NA, cor_w_rank = NA)
for (i in seq_len(nrow(res))) {
  x <- generate_interaction_probs(n_ind = sample(5:20, 1), id_bias = res$idbias[i], rank_bias = res$rankbias[i])
  # generate 10000 'interactions'
  s <- sample(seq_len(nrow(x)), 10000, TRUE, prob = x[, "final"])

  xtab <- as.numeric(table(as.numeric(x[s, 1:2])))
  res$cv_id[i] <- sd(xtab) / mean(xtab)

  xtab <- table((x[s, 2] - x[s, 1]))
  rds <- as.numeric(names(xtab))

  # coef(lm(scale(xtab) ~ scale(rds)))[2]
  # res$cor_w_rank[i] <- cor(xtab, as.numeric(names(xtab)))

  rankdiff <- x[, 2] - x[, 1]
  interactprob <- x[, "final"]
  res$cor_w_rank[i] <- cor(rankdiff, interactprob)


}

plot(res$idbias, res$cv_id)
cor(res$idbias, res$cv_id)

plot(res$rankbias, res$cor_w_rank)
cor(res$rankbias, res$cor_w_rank)



test_that("multiplication works", {
  expect_true(cor(res$idbias, res$cv_id) > 0.3)
  expect_true(cor(res$rankbias, res$cor_w_rank) < -0.3)
})
