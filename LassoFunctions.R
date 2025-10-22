#  Standardize X and Y: center both X and Y; scale centered X
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X, Y){
  n <- nrow(X)
  
  #  Center Y
  Ymean <- mean(Y)
  Ytilde <- Y - Ymean
  
  #  Center and scale X
  Xmeans <- colMeans(X)
  Xcentered <- X - matrix(Xmeans, n, ncol(X), byrow = T)
  weights <- sqrt(diag(crossprod(Xcentered) / n))
  Xtilde <- Xcentered %*% diag(1 / weights)
  
  # Return:
  # Xtilde - centered and appropriately scaled X
  # Ytilde - centered Y
  # Ymean - the mean of original Y
  # Xmeans - means of columns of X (vector)
  # weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

#  Soft-thresholding of a scalar a at level lambda 
# [OK to have vector version as long as works correctly on scalar; will only test on scalars]
soft <- function(a, lambda){
  return(sign(a) * max(abs(a) - lambda, 0))
}

#  Calculate objective function of lasso given current values of Xtilde, Ytilde, beta and lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
lasso <- function(Xtilde, Ytilde, beta, lambda){
 n <- length(Ytilde)
 return(
   sum(((Ytilde - Xtilde %*% beta) * (2 * n)^(-1/2) )^2) + lambda * sum(abs(beta))
 )
}

#  Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1 (vector)
# lamdba - tuning parameter
# beta_start - p vector, an optional starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.001){
  # Conver Xtilde to a matrix
  Xtilde <- as.matrix(Xtilde)
  
  #  Check that n is the same between Xtilde and Ytilde
  if (nrow(Xtilde) != length(Ytilde)) {
    stop("Xtilde must have number of rows equal to the length of Ytilde")
  }
  #  Check that lambda is non-negative
  if (lambda < 0) {
    stop("Lambda must be non-negative")
  }
  #  Check for starting point beta_start. 
  # If none supplied, initialize with a vector of zeros.
  # If supplied, check for compatibility with Xtilde in terms of p
  if (is.null(beta_start)){
    beta_start = numeric(ncol(Xtilde))
  } else if (ncol(Xtilde) != length(beta_start)) {
    stop("beta_start must have length equal to the number of rows in Xtilde")
  }
  
  #  Coordinate-descent implementation. 
  # Stop when the difference between objective functions is less than eps for the first time.
  # For example, if you have 3 iterations with objectives 3, 1, 0.99999,
  # your should return fmin = 0.99999, and not have another iteration
  
  # Initialize starting points
  beta <- beta_start
  fprevious <- lasso(Xtilde, Ytilde, beta, lambda)
  p <- ncol(Xtilde)
  n <- nrow(Xtilde)
  resid <- Ytilde - Xtilde %*% beta
  
  # Start loop, break at end if difference between new and old is less than eps
  while(TRUE){
    beta_old <- beta
    for (j in 1:p){
      beta_old <- beta[j]
      beta[j] <- soft((1 / n) * crossprod(Xtilde[ , j], resid + Xtilde[ , j] * beta[j]), lambda)
      resid <- resid - Xtilde[, j] * (beta[j] - beta_old)
    }
    
    # Check break condition
    fmin <- lasso(Xtilde, Ytilde, beta, lambda)
    
    if(abs(fprevious - fmin) < eps) {
      break
    }
    
    # Update previous if no break is needed
    fprevious <- fmin
  }
  
  
  # Return 
  # beta - the solution (a vector)
  # fmin - optimal function value (value of objective at beta, scalar)
  return(list(beta = beta, fmin = fmin))
}

#  Fit LASSO on standardized data for a sequence of lambda values. Sequential version of a previous function.
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence,
#             is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized_seq <- function(Xtilde, Ytilde, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  p <- ncol(Xtilde)
  n <- nrow(Xtilde)
  
  #  Check that n is the same between Xtilde and Ytilde
  if (nrow(Xtilde) != length(Ytilde)) {
    stop("Xtilde must have number of rows equal to the length of Ytilde")
  }
  #  Check for the user-supplied lambda-seq (see below)
  # If lambda_seq is supplied, only keep values that are >= 0,
  # and make sure the values are sorted from largest to smallest.
  # If none of the supplied values satisfy the requirement,
  # print the warning message and proceed as if the values were not supplied.
  if (!is.null(lambda_seq)) {
    lambda_seq <- sort(lambda_seq[lambda_seq >= 0], decreasing = T)
  }
  if (length(lambda_seq) == 0) {
  # If lambda_seq is not supplied, calculate lambda_max 
  # (the minimal value of lambda that gives zero solution),
  # and create a sequence of length n_lambda as
    lambda_max <- max(abs(crossprod(Xtilde, Ytilde)/n))
    lambda_seq <- exp(seq(log(lambda_max), log(0.01), length = n_lambda))
  }
  #  Apply fitLASSOstandardized going from largest to smallest lambda 
  # (make sure supplied eps is carried over). 
  # Use warm starts strategy discussed in class for setting the starting values.
  
  # Initialize beta_mat and fmin_vec
  beta_mat <- matrix(0, nrow = p, ncol = length(lambda_seq))
  fmin_vec <- numeric(length(lambda_seq))
  
  for (j in 1:length(lambda_seq)) {
    LASSOfit <- fitLASSOstandardized(Xtilde = Xtilde, Ytilde = Ytilde,
                                     lambda = lambda_seq[j], beta_start = beta_mat[ , max(j - 1, 1)],
                                     eps = eps)
    beta_mat[ , j] <- LASSOfit$beta
    fmin_vec[j] <- LASSOfit$fmin
  }
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, fmin_vec = fmin_vec))
}

#  Fit LASSO on original data using a sequence of lambda values
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  #  Center and standardize X,Y based on standardizeXY function
  stand <- standardizeXY(X, Y)
  Xtilde <- stand$Xtilde
  Ytilde <- stand$Ytilde
  
  #  Fit Lasso on a sequence of values using fitLASSOstandardized_seq
  # (make sure the parameters carry over)
  LASSOStand <- fitLASSOstandardized_seq(Xtilde = Xtilde, Ytilde = Ytilde, lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
  #  Perform back scaling and centering to get original intercept and coefficient vector
  # for each lambda
  beta_mat <- LASSOStand$beta_mat / stand$weights
  beta0_vec <- stand$Ymean - colSums(beta_mat * stand$Xmeans)
  lambda_seq <- LASSOStand$lambda_seq
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec))
}


#  Fit LASSO and perform cross-validation to select the best fit
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# k - number of folds for k-fold cross-validation, default is 5
# fold_ids - (optional) vector of length n specifying the folds assignment (from 1 to max(fold_ids)), if supplied the value of k is ignored 
# eps - precision level for convergence assessment, default 0.001
cvLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, k = 5, fold_ids = NULL, eps = 0.001){
  #  Fit Lasso on original data using fitLASSO
  LASSOfitWHOLE <- fitLASSO(X, Y, lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
  lambda_seq <- LASSOfitWHOLE$lambda_seq
  #  If fold_ids is NULL, split the data randomly into k folds.
  # If fold_ids is not NULL, split the data according to supplied fold_ids.
  if (is.null(fold_ids)){
  fold_ids <- sample(rep_len(1:k, nrow(X)))
  } else {
    k <- max(fold_ids)
  }
  
  # Stop if k > n
  if(k > nrow(X)) {
    stop("Number of folds, k, must be less than or equal to the number of samples, n.")
  }
  
  #  Calculate LASSO on each fold using fitLASSO,
  # and perform any additional calculations needed for CV(lambda) and SE_CV(lambda)
  CV_mat <- matrix(0, ncol = length(lambda_seq), nrow = k)
  for (fold in 1:k) {
    foldfit <- fitLASSO(X[fold_ids != fold, ], Y[fold_ids != fold], lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
    CV_mat[fold, ] <- colMeans(
       (outer(Y[fold_ids == fold], foldfit$beta0_vec, '-') - apply(foldfit$beta_mat, 2, \(beta_i) {X[fold_ids == fold, ] %*% beta_i}))^2
    )
  }
  #  Find lambda_min
  cvm <- colMeans(CV_mat)
  lambda_min <- lambda_seq[which.min(cvm)]
  #  Find lambda_1SE
  cvse <- apply(CV_mat, 2, \(cv_i) {sd(cv_i)/sqrt(k)})
  lambda_1se <- max(
    lambda_seq[cvm <= min(cvm) + cvse[which.min(cvm)]]
  )
  
  # Return output
  # Output from fitLASSO on the whole data
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  # fold_ids - used splitting into folds from 1 to k (either as supplied or as generated in the beginning)
  # lambda_min - selected lambda based on minimal rule
  # lambda_1se - selected lambda based on 1SE rule
  # cvm - values of CV(lambda) for each lambda
  # cvse - values of SE_CV(lambda) for each lambda
  beta_mat <- LASSOfitWHOLE$beta_mat
  beta0_vec <- LASSOfitWHOLE$beta0_vec
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec, fold_ids = fold_ids, lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse))
}

