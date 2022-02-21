
#' @title Looney's Corrected z-test
#' @description Performs a corrected z-test for partially matched samples with a matrix or two-vector input that includes missing data.
#'
#' @param X A 2 x k (or k x 2) matrix with partially matched pairs. If parameter y is specified, then X should instead be a vector with length equal to y.
#' @param y A vector of data paired with the X vector of data. Length should be equal to length of X. Default is NULL if X is a matrix of the correct dimensions.
#' @param alt The test's alternative hypothesis. Deaults to 'two.sided'.
#'
#' @return The test statistic (z.stat) and p-value (p.val) of the test.
#' @export
#'
#' @examples
#' #INSTALL THE LIBRARY 'mvtnorm'
#' #Generate bivariate gaussian sample data
#'
#' require(mvtnorm)
#' set.seed(123) # simulation seed
#'
#' rho <- c(0.2, 0.7) # correlation test values
#' n <- c(20, 50) # test sample sizes
#' sig2 <- c(1, 1) # homoscedastic test standard deviations
#'
#'
#' #Homoscedastic data with n=50, rho=0.7
#'
#' homo_50_0.7 <- rmvnorm(n[2], mean=c(0, 0), sigma=matrix(c(sig2[1]^2, rho[2]*sig2[1]*sig2[2], rho[2]*sig2[1]*sig2[2], sig2[2]^2), nrow=2, byrow=T))
#'
#' #Generate indices to make NA
#' NA_index50 <- sample(1:50, 20) # 20 NA indices generated
#'
#'
#' #Use generated indices to make random entries NA
#'
#' homo_50_0.7[,1][NA_index50[1:10]] <- NA
#' homo_50_0.7[,2][NA_index50[11:20]] <- NA
#'
#' #We now have a sample data set called homo_50_0.7 that is homoscedastic, n=50, and rho=0.7. A complete dataset for testing can be found in the vignette.
#'
#'
#' #Example of function use, two vectors
#'
#' looney_z_corr(homo_50_0.7[,1], homo_50_0.7[,2], alt='less')
#'
#' #Example of function use, one matrix
#'
#' looney_z_corr(homo_50_0.7, alt='less')
#'
looney_z_corr <- function(X, y=NULL, alt = c('two.sided', 'less', 'greater')) { # Where X is a matrix of partially matched pairs

  if (missing(alt)) {
    alt = 'two.sided'
  }

  if (is.null(dim(X)) & !is.null(y)) {
    if (length(X) != length(y)) {
      print("Warning: input vectors not of equal length, may result in computational errors")
    }
    X <- matrix(c(X, y), nrow=2, byrow=T)
    # This is the case where X and y arguments are vectors of partially matched data
    # redefine X as a 2 x k matrix
  }

  if (dim(X)[1] != 2) { # Take transpose of incorrectly oriented input matrix (should be 2 x k)
    X <- t(X)
  }

  if (sum(colSums(is.na(X))==2)!=0) { # condition is if there is at least one column with double NA entries
    print("Warning: pairs found with two NA values. Omitting these pairs.")
    X <- X[,-c(which(colSums(is.na(X))==2))] # remove columns with double NA entries
  }

  r1_NA <- which(is.na(X[1,])) # row 1 indices that contain NA
  r2_NA <- which(is.na(X[2,])) # row 2 indices that contain NA

  if (length(r2_NA) < 3 | length(r1_NA)  < 3) {
    print("Warning: not enough NA values in at least one row (fewer than 3), ommitting all unmatched pairs and using matched pair t-test")
    X <- X[,-c(which(colSums(is.na(X))==1))]
    t_test <- t.test(X[1,], X[2,], paired=T, alternative = alt)
    return(t_test)
  }

  n1_dat <- X[,colSums(is.na(X))==0] # Matrix of complete pairs

  if (dim(n1_dat)[2] < 3) {
    print("Warning: not enough complete pairs, perform two-sample t-test with lone values")
    n2_dat <- X[1,][is.na(X[2,])] # Lone values in first row
    n3_dat <- X[2,][is.na(X[1,])] # Lone values in second row
    t_test <- t.test(n2_dat, n3_dat, alternative = alt)
    return(t_test)
  }

  n1_n2_dat <- X[1,][!is.na(X[1,])] # First row of X without NA
  n1_n3_dat <- X[2,][!is.na(X[2,])] # Second row of X without NA


  numerator <- mean(n1_n2_dat) - mean(n1_n3_dat)
  if (dim(n1_dat)[2] < 2) { # Covariance is undefined for zero or one pair
    denominator <- sqrt(var(n1_n2_dat)/(length(n1_n2_dat)) + var(n1_n3_dat)/(length(n1_n3_dat)))
  }
  else { # Case where covariance is defined
    denominator <- sqrt( var(n1_n2_dat)/(length(n1_n2_dat)) + var(n1_n3_dat)/(length(n1_n3_dat))
                         - (2*dim(n1_dat)[2]*cov(n1_dat[1,], n1_dat[2,])/(length(n1_n2_dat)*length(n1_n3_dat))))
  }

  Z <- numerator/denominator
  if (alt == 'less') {
    p.val <- pnorm(Z)
  } else if (alt == 'greater') {
    p.val <- 1 - pnorm(Z)
  } else {
    if (Z <= 0){
      p.val <- 2*pnorm(Z)
    } else {
      p.val <- 2*(1-pnorm(Z))
    }
  }

  return(list(z.stat=Z, p.val=p.val))

}

























