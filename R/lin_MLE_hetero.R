
#' @title Lin and Stiver's MLE Test with Heteroscedasticity
#' @description Performs an MLE based test for heteroscedastic partially matched samples with a matrix or two-vector input that includes missing data.
#'
#' @param X A 2 x k (or k x 2) matrix with partially matched pairs. If parameter y is specified, then X should instead be a vector with length equal to y.
#' @param y A vector of data paired with the X vector of data. Length should be equal to length of X. Default is NULL if X is a matrix of the correct dimensions.
#' @param alt The test's alternative hypothesis. Deaults to 'two.sided'.
#'
#' @return The test statistic (z.stat) and the p-value (p.val) of the test.
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
#' sig1 <- c(1, 3) # heteroscedastic test standard deviations
#'
#'
#' # Heterosedastic data with n=50, rho=0.7
#'
#' hetero_50_0.7 <- rmvnorm(n[2], mean=c(0, 0), sigma=matrix(c(sig1[1]^2, rho[2]*sig1[1]*sig1[2], rho[2]*sig1[1]*sig1[2], sig1[2]^2), nrow=2, byrow=T))
#'
#'
#' #Generate indices to make NA
#' NA_index50 <- sample(1:50, 20) # 20 NA indices generated
#'
#'
#' #Use generated indices to make random entries NA
#'
#' hetero_50_0.7[,1][NA_index50[1:10]] <- NA
#' hetero_50_0.7[,2][NA_index50[11:20]] <- NA
#'
#'
#' #We now have a sample data set called hetero_50_0.7 that is heteroscedastic, n=50, and rho=0.7. A complete dataset for testing can be found in the vignette.
#'
#'
#' #Example of function use, two vectors
#'
#' lin_MLE_hetero(hetero_50_0.7[,1], hetero_50_0.7[,2], alt='less')
#'
#' #Example of function use, one matrix
#'
#' lin_MLE_hetero(hetero_50_0.7, alt='less')
#'
lin_MLE_hetero <- function(X, y=NULL, alt= c('two.sided', 'less', 'greater')) {

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
  n2_dat <- X[1,][is.na(X[2,])] # Lone values in first row
  n3_dat <- X[2,][is.na(X[1,])] # Lone values in second row

  if (dim(n1_dat)[2] < 3) {
    print("Warning: not enough complete pairs, perform two-sample t-test with lone values")
    t_test <- t.test(n2_dat, n3_dat, alternative = alt)
    return(t_test)
  }

  n1 <- dim(n1_dat)[2]
  n2 <- length(n2_dat) # Lone values in first row length
  n3 <- length(n3_dat) # Lone values in second row length

  r <- cov(n1_dat[1,], n1_dat[2,])/(sd(n1_dat[1,])*sd(n1_dat[2,]))
  g <- n1*((n1 + n2 + n3*cov(n1_dat[1,], n1_dat[2,])/var(n1_dat[2,])))/((n1 + n2)*(n1 + n3)-n2*n3*r^2)
  f <- n1*((n1 + n3 + n2*cov(n1_dat[1,], n1_dat[2,])/var(n1_dat[1,])))/((n1 + n2)*(n1 + n3)-n2*n3*r^2)
  V1 <- ((f^2/n1 + (1-f)^2/n2)*var(n1_dat[1,])*(n1-1)+(g^2/n1 + (1-g)^2/n3)*var(n1_dat[2,])*(n1-1)-2*f*g*cov(n1_dat[1,], n1_dat[2,])*(n1-1)/n1)/(n1-1)

  Z_LS <- (f*(mean(n1_dat[1,]) - mean(n2_dat)) - g*(mean(n1_dat[2,]) - mean(n3_dat)) + mean(n2_dat) - mean(n3_dat))/sqrt(V1) # test statistic

  if (alt == 'less') {
    p.val <- pt(Z_LS, n1)
  } else if (alt == 'greater') {
    p.val <- 1 - pt(Z_LS, n1)
  } else {
    if (Z_LS <= 0){
      p.val <- 2*pt(Z_LS, n1)
    } else {
      p.val <- 2*(1-pt(Z_LS, n1))
    }
  }

  return(list(z.stat=Z_LS, p.val=p.val))
}




























