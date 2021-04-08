
library(EloRating)
xdata <- randomsequence(nID = 12, avgIA = 100)
mat <- creatematrix(winners = xdata$seqdat$winner, losers = xdata$seqdat$loser)

# orders
o1 <- sample(colnames(mat))
m1 <- mat[o1, o1]
o2 <- sample(colnames(mat))
m2 <- mat[o2, o2]
o3 <- sample(colnames(mat))
m3 <- mat[o3, o3]

i0 <- initial_values(mat, chains = 10)
i1 <- initial_values(m1, chains = 10)
i2 <- initial_values(m2, chains = 10)
i3 <- initial_values(m3, chains = 10)


c0 <- cor(matrix(unlist(i0), ncol = 10))
c1 <- cor(matrix(unlist(i1), ncol = 10))
c2 <- cor(matrix(unlist(i2), ncol = 10))
c3 <- cor(matrix(unlist(i3), ncol = 10))
test_that("multiplication works", {
  expect_true(mean(c0[upper.tri(c0)]) > 0.9)
  expect_true(mean(c1[upper.tri(c1)]) > 0.9)
  expect_true(mean(c2[upper.tri(c2)]) > 0.9)
  expect_true(mean(c3[upper.tri(c3)]) > 0.9)
})


mx <- mat
colnames(mx) <- NULL
rownames(mx) <- NULL

ix <- initial_values(mx, chains = 10)
cx <- cor(cbind(matrix(unlist(i0), ncol = 10), matrix(unlist(ix), ncol = 10)))

test_that("multiplication works", {
  expect_true(mean(cx[upper.tri(cx)]) > 0.9)
})



# remove interactions of one
xid <- sample(ncol(mat), 1)
mx <- mat
mx[xid, ] <- 0
mx[, xid] <- 0

xres <- initial_values(mx)

test_that("multiplication works", {
  xres <- unlist(lapply(xres, function(x)x$EloStart[xid] == 0))
  expect_true(all(xres))
})







