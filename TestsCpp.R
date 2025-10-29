
# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)
library(testthat)

# Source your C++ funcitons
sourceCpp("LassoInC.cpp")

# Source your LASSO functions from HW4 (make sure to move the corresponding .R file in the current project folder)
source("LassoFunctions.R")

# Do at least 2 tests for soft-thresholding function below. You are checking output agreements on at least 2 separate inputs
#################################################
test_that("soft_c returns same values as R version", {
  expect_equal(soft(5,3), soft_c(5,3))
  expect_equal(soft(5,-4), soft_c(5,-4))
  expect_equal(soft(5,10), soft_c(5,10))
  expect_equal(soft(-5,10), soft_c(-5,10))
  expect_equal(soft(0,10), soft_c(0,10))
  expect_equal(soft(0,0), soft_c(0,0))
  expect_equal(soft(-1,-1), soft_c(-1,-1))
})

# Do at least 2 tests for lasso objective function below. You are checking output agreements on at least 2 separate inputs
#################################################
test_that("lasso_c returns same values as R version", {
  set.seed(12)
  X <- matrix(rnorm(10 * 2), 10, 2)
  Y <- rnorm(10)
  beta <- 1:2
  lambda <- 2
  expect_equal(lasso_c(X,Y, beta, lambda), lasso(X,Y, beta, lambda))
  
  
  X <- matrix(rnorm(20 * 5), 20, 5)
  Y <- rnorm(20)
  beta <- rnorm(5)
  lambda <- 0.5
  expect_equal(lasso_c(X,Y, beta, lambda), lasso(X,Y, beta, lambda))
  
  X <- matrix(rnorm(100 * 50), 100, 50)
  Y <- rnorm(100)
  beta <- rnorm(50)
  lambda <- 1
  expect_equal(lasso_c(X,Y, beta, lambda), lasso(X,Y, beta, lambda))
})

# Do at least 2 tests for fitLASSOstandardized function below. You are checking output agreements on at least 2 separate inputs
#################################################
test_that("fitLassostandardized returns same estimates in C and R", {
  set.seed(12)
  X <- matrix(rnorm(50), 50)
  Y <- X %*% rnorm(1)
  stand <- standardizeXY(X, Y)
  Xtilde <- stand$Xtilde
  Ytilde <- stand$Ytilde
  
  expect_equal(
    fitLASSOstandardized(Xtilde, Ytilde, lambda = 0)$beta,
    as.vector(fitLASSOstandardized_c(Xtilde, Ytilde, lambda = 0, c(1)))
  )
  
  
  set.seed(13)
  X <- matrix(rnorm(50 * 3), 50, 3)
  Y <- X %*% c(10, 3, 6)
  stand <- standardizeXY(X, Y)
  Xtilde <- stand$Xtilde
  Ytilde <- stand$Ytilde
  
  expect_equal(
    fitLASSOstandardized(Xtilde, Ytilde, lambda = 1)$beta,
    as.vector(fitLASSOstandardized_c(Xtilde, Ytilde, lambda = 1, c(0,0,0)))
  )
  
  set.seed(13)
  X <- matrix(rnorm(50 * 3), 50, 3)
  Y <- X %*% c(-1,-2,-3)
  stand <- standardizeXY(X, Y)
  Xtilde <- stand$Xtilde
  Ytilde <- stand$Ytilde
  
  expect_equal(
    fitLASSOstandardized(Xtilde, Ytilde, lambda = 1)$beta,
    as.vector(fitLASSOstandardized_c(Xtilde, Ytilde, lambda = 1, c(0,0,0)))
  )
  
  set.seed(11)
  X <- matrix(rnorm(100 * 10), 100, 10)
  Y <- X %*% rnorm(10)
  stand <- standardizeXY(X, Y)
  Xtilde <- stand$Xtilde
  Ytilde <- stand$Ytilde
  
  expect_equal(
    fitLASSOstandardized(Xtilde, Ytilde, lambda = 0.75)$beta,
    as.vector(fitLASSOstandardized_c(Xtilde, Ytilde, lambda = 0.75, rep(0, 10)))
  )
})

# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################
set.seed(14)
X <- matrix(rnorm(50 * 3), 50, 3)
Y <- X %*% c(10, 3, 6)
stand <- standardizeXY(X, Y)
Xtilde <- stand$Xtilde
Ytilde <- stand$Ytilde

microbenchmark::microbenchmark(
  fitLASSOstandardized(Xtilde, Ytilde, lambda = 1),
  fitLASSOstandardized_c(Xtilde, Ytilde, lambda = 1, c(0,0,0))
)
# Median R time 226.95 microseconds
# Median C++ time 13.10 microseconds


# Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
#################################################
test_that("fitLASSOstandardized_seq_c works properly",{
  set.seed(15)
  X <- matrix(rnorm(50 * 3), 50, 3)
  Y <- X %*% c(10, 3, 6)
  stand <- standardizeXY(X, Y)
  Xtilde <- stand$Xtilde
  Ytilde <- stand$Ytilde
  out1 <- fitLASSOstandardized_seq(Xtilde, Ytilde)
  
  expect_equal(fitLASSOstandardized_seq_c(Xtilde, Ytilde, out1$lambda_seq), out1$beta_mat)
  
  set.seed(16)
  X <- matrix(rnorm(100 * 5), 100, 5)
  Y <- X %*% c(10, 3, 6, 4, 2)
  stand <- standardizeXY(X, Y)
  Xtilde <- stand$Xtilde
  Ytilde <- stand$Ytilde
  out2 <- fitLASSOstandardized_seq(Xtilde, Ytilde)
  
  expect_equal(fitLASSOstandardized_seq_c(Xtilde, Ytilde, out2$lambda_seq), out2$beta_mat)
  
  set.seed(17)
  X <- matrix(rnorm(100 * 50), 100, 50)
  Y <- X %*% rnorm(50)
  stand <- standardizeXY(X, Y)
  Xtilde <- stand$Xtilde
  Ytilde <- stand$Ytilde
  out2 <- fitLASSOstandardized_seq(Xtilde, Ytilde)
  
  expect_equal(fitLASSOstandardized_seq_c(Xtilde, Ytilde, out2$lambda_seq), out2$beta_mat)
})
# Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
######################################################################
set.seed(16)
X <- matrix(rnorm(100 * 5), 100, 5)
Y <- X %*% c(10, 3, 6, 4, 2)
stand <- standardizeXY(X, Y)
Xtilde <- stand$Xtilde
Ytilde <- stand$Ytilde
out2 <- fitLASSOstandardized_seq(Xtilde, Ytilde)
microbenchmark::microbenchmark(
  fitLASSOstandardized_seq(Xtilde, Ytilde),
  fitLASSOstandardized_seq_c(Xtilde, Ytilde, out2$lambda_seq)
)
# R version has median time of 10800.25 microseconds
# C++ version has median time of 388.20 microseconds


# Tests on riboflavin data
##########################
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression

# Make sure riboflavin$x is treated as matrix later in the code for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]

# Standardize the data
out <- standardizeXY(riboflavin$x, riboflavin$y)

# This is just to create lambda_seq, can be done faster, but this is simpler
outl <- fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, n_lambda = 30)

# The code below should assess your speed improvement on riboflavin data
microbenchmark::microbenchmark(
  fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, outl$lambda_seq),
  fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outl$lambda_seq),
  times = 10
)
# Median time for R is 2806.91395
# Median time for C++ is 69.10695