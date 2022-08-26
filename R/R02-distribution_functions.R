# Distribution functions for data generation
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you
# have any questions or comments
# Last updated 2022-08-25

# Table of contents
# 1) Categorical outcomes
#   1.1) sim.bernoulli
#   1.2) sim.binomial
#   1.3) sim.beta_binomial
#   1.4) sim.ordered_probit
# 2) Continuous outcomes
#   2.1) sim.beta
#   2.2) sim.gaussian
#   2.3) sim.normal
#   2.4) sim.log_normal
#   2.5) sim.student_t
# 3) Multivariate outcomes
#   3.1) sim.multivariate_normal
# 4) Miscellaneous

#### 1) Categorical outcomes ####

#### 1.1) sim.bernoulli ####
#' Simulate From a Bernoulli Distribution
#'
#' Generates draws from a Bernoulli distribution.
#'
#' @param theta A vector of probabilities between 0 and 1.
#'
#' @returns A vector of binary values, either 0 or 1.
#'
#' @examples
#'
#' # For reproducibility
#' set.seed( 111 )
#' # Generate three draws
#' theta <- c( 0.1, 0.5, 0.9 )
#' sim.bernoulli( theta )
#'
#' @export

sim.bernoulli <- function( theta ) {

  # Check for valid inputs
  if ( any( min( theta ) < 0 | max( theta ) > 1 ) ) {
    stop( 'Argument "theta" must be between 0 and 1' )
  }

  # Number of observations
  n <- length( theta )

  # Generate draws
  y <- rbinom( n, 1, theta )

  return( y )
}

#### 1.2) sim.binomial ####
#' Simulate From a Binomial Distribution
#'
#' Generates draws from a Binomial distribution.
#'
#' @param theta A vector of probabilities (between 0 and 1).
#' @param nu A vector of integers (1 or higher), the maximum
#'   possible number of counts/successes.
#'
#' @returns A vector of counts between 0 and \code{nu}.
#'
#' @examples
#'
#' # For reproducibility
#' set.seed( 111 )
#' # Generate three draws
#' theta <- c( 0.1, 0.5, 0.9 )
#' nu <- c( 10, 20, 30 )
#' sim.binomial( theta, nu )
#'
#' @export

sim.binomial <- function( theta, nu ) {

  # Check for valid inputs

  if ( any( min( theta ) < 0 | max( theta ) > 1 ) ) {
    stop( 'Argument "theta" must be between 0 and 1' )
  }

  if ( min( nu ) < 1 ) {
    stop( 'Argument "nu" must be 1 or higher' )
  }

  # Ensure 'nu' is an integer
  nu <- round( nu )

  # Determine max number of observations
  n <- max( c(
    length( theta ),
    length( nu )
  ) )

  # Ensure input vectors are same length
  theta <- rep_len( theta, n )
  nu <- rep_len( nu, n )

  # Generate draws
  y <- rbinom( n, nu, theta )

  return( y )
}

#### 1.3) sim.beta_binomial ####
#' Simulate From a Beta-Binomial Distribution
#'
#' Generate draws from a Beta-Binomial distribution.
#'
#' @param mu A vector of location parameters (between 0 and 1).
#' @param phi A vector of dispersion parameters (greater than 0).
#' @param nu A vector of integers (1 or higher), the maximum
#'   possible number of counts/successes.
#'
#' @returns A vector of counts between 0 and \code{nu}.
#'
#' @examples
#'
#' # For reproducibility
#' set.seed( 111 )
#' # Generate three draws
#' mu <- c( 0.1, 0.5, 0.9 )
#' phi <- c( 10, 35, 70 )
#' nu <- c( 10, 20, 30 )
#' sim.beta_binomial( mu, phi, nu )
#'
#' @export

sim.beta_binomial <- function( mu, phi, nu ) {

  # Check for valid inputs

  if ( any( min( mu ) < 0 | max( mu ) > 1 ) ) {
    stop( 'Argument "mu" must be between 0 and 1' )
  }

  if ( min( phi ) <= 0 ) {
    stop( 'Argument "phi" must be greater than 0' )
  }

  if ( min( nu ) < 1 ) {
    stop( 'Argument "nu" must be 1 or higher' )
  }

  # Ensure 'nu' is an integer
  nu <- round( nu )

  # Determine max number of observations
  n <- max( c(
    length( mu ),
    length( phi ),
    length( nu )
  ) )

  # Ensure input vectors are same length
  mu <- rep_len( mu, n )
  phi <- rep_len( phi, n )
  nu <- rep_len( nu, n )

  # Generate probabilities from Beta distribution
  theta <- rbeta(
    n,
    mu * phi,
    (1 - mu ) * phi
  )

  # Generate draws from Binomial distribution
  y <- rbinom( n, nu, theta )

  return( y )
}

#### 1.3) sim.ordered_probit ####
#' Simulate From an Ordered-Probit Model
#'
#' Generate draws from a Ordered-Probit model (also
#' known as a thresholded cumulative normal model).
#'
#' @param eta A vector of shift parameters.
#' @param tau A N x (K-1) matrix of thresholds where
#'   \code{tau[i,j] > tau[i,j-1]}.
#' @param sigma A vector of dispersion parameters
#'   (greater than 0), the latent distribution's
#'   standard deviation.
#'
#' @details
#'
#' Given N observations, let \eqn{Y = y_i} be the
#' observed ordinal ranking for a random variable
#' \eqn{Y = \{1, ..., K\} }, with \eqn{i = \{1, ..., N\} }.
#' The ordered-probit model assumes that
#' \eqn{y_i} is produced via categorization of a
#' latent (unobservable) continuous variable,
#' \eqn{Z = z_i}. The distribution for the latent
#' variable is partitioned into \eqn{K} regions via
#' \eqn{K - 1} thresholds \eqn{\tau_1, ..., \tau_{K-1}}
#' where \eqn{\tau_j > \tau_{j-1}}. Furthermore,
#' we assume that \eqn{\tau_K = \infty} and
#' \eqn{\tau_0 = -\infty}. The region that
#' \eqn{z_i} falls within determines the value
#' of \eqn{y_i}. For example, if \eqn{z_i} is
#' between \eqn{\tau_1} and \eqn{\tau_2}, then
#' \eqn{y_i = 2}.
#'
#' The ordered-probit model assumes that \eqn{Z}
#' follows a standard normal distribution (mean of 0
#' and standard devation of 1), which means:
#'
#' \deqn{ P( y_i = k) = \Phi( tau_k - eta_i ) - \Phi( tau_{k-1} - eta_i ). }
#'
#' Here \eqn{ \Phi(x) } is the cumulative distribution
#' function for the standard normal model and \eqn{\eta_i}
#' is a shift parameter (which can be a linear combination
#' of predictors and regression coefficients).
#'
#' @returns A vector of integers between 1 and K.
#'
#' @examples
#'
#' # For reproducibility
#' set.seed( 111 )
#'
#' # Number of thresholds determines
#' # number of ordered categories
#'
#' # 2 thresholds = 3 categories
#' sim.ordered_probit( 0, c( -1, 1 ) )
#'
#' # 4 thresholds = 5 categories
#' sim.ordered_probit( 0, c( -.84, -.25, .25, .84 ) )
#'
#' # 2 draws
#' eta <- c( 0, -.25 )
#' tau <- rbind( c( -1, 1 ), c( -.5, .5) )
#' sim.ordered_probit( eta, tau )
#'
#' @export

sim.ordered_probit <- function( eta, tau, sigma = 1 ) {

  # Ensure 'tau' is a matrix
  if ( !is.matrix(tau) ) {
    tau <- rbind( tau )
  }

  # Check for valid inputs

  if ( min( sigma ) <= 0 ) {
    stop( 'Argument "sigma" must be greater than 0' )
  }

  for ( k in 2:ncol(tau) ) {
    if ( any( tau[,k] < tau[,k-1] ) ) {
      stop( paste0(
        'Argument "tau" must be ordered from lowest to highest by column'
      ) )
    }
  }

  # Determine max number of observations
  n <- max( c(
    length( eta ),
    nrow( tau ),
    length( sigma )
  ) )

  # Ensure input vectors are same length
  eta <- rep_len( eta, n )
  if ( n > 1 ) {
    tau <- tau[ rep_len( 1:nrow(tau), n ), ]
  }
  sigma <- rep_len( sigma, n )

  N <- nrow( tau )
  K <- ncol( tau ) + 1

  theta <- matrix( NA, N, K )

  if ( N == 1 ) theta <- rbind( theta )

  for ( k in 1:K ) {

    if ( k == 1 ) {
      theta[,k] <- pnorm( tau[,k] - eta, 0, sigma )
    }

    if ( k == K ) {
      theta[,k] <- 1 - pnorm( tau[,k-1] - eta, 0, sigma )
    }

    if ( k > 1 & k < K ) {
      theta[,k] <-
        pnorm( tau[,k] - eta, 0, sigma ) -
        pnorm( tau[,k-1] - eta, 0, sigma )
    }

  }

  y <- sapply( 1:N, function(n) {
    sample( 1:K, size = 1, prob = theta[n,] )
  } )

  return( y )
}

#### 2) Continuous outcomes ####

#### 2.1) sim.beta ####
#' Simulate From a Beta Distribution
#'
#' Generates draws from a Beta distribution.
#'
#' @param mu A vector of location parameters (between 0 and 1).
#' @param phi A vector of dispersion parameters (greater then 0).
#'
#' @returns A vector of values between 0 and 1.
#'
#' @examples
#'
#' # For reproducibility
#' set.seed( 111 )
#' # Generate three draws
#' mu <- c( 0.5, 0.1, 0.9 )
#' phi <- c( 2, 5, 10 )
#' sim.beta( mu, phi )
#'
#' # Convert shape parameters to location/dispersion
#' shape_param <- c( alpha = 2, beta = 5 )
#' loc_disp_param <- calc.beta_binomial( shape_param, 'mu and phi' )
#' sim.beta( loc_disp_param[1], loc_disp_param[2] )
#'
#' @export

sim.beta <- function( mu, phi ) {

  # Check for valid inputs

  if ( any( min( mu ) < 0 | max( mu ) > 1 ) ) {
    stop( 'Argument "mu" must be between 0 and 1' )
  }

  if ( min( phi ) <= 0 ) {
    stop( 'Argument "phi" must greater than 0' )
  }

  # Determine max number of observations
  n <- max( c(
    length( mu ),
    length( phi )
  ) )

  # Ensure input vectors are same length
  mu <- rep_len( mu, n )
  phi <- rep_len( phi, n )

  # Generate draws
  y <- rbeta( n, mu * phi, (1 - mu ) * phi )

  return( y )
}

#### 2.2) sim.gaussian ####
#' Simulate From a Gaussian Distribution
#'
#' Generate draws from a Gaussian (Normal) distribution.
#'
#' @param mu A vector of location parameters, the distribution
#'   mean.
#' @param sigma A vector of dispersion parameters (greater
#'   than 0), the distribution standard deviation.
#'
#' @return A vector of values.
#'
#' @examples
#'
#' # For reproducibility
#' set.seed( 111 )
#' # Generate two draws
#' mu <- c( 0, 100 )
#' sigma <- c( 1, 15 )
#' round( sim.gaussian( mu, sigma ), 2 )
#'
#' @export

sim.gaussian <- function( mu, sigma ) {

  # Check for valid inputs

  if ( min( sigma ) <= 0 ) {
    stop( 'Argument "sigma" must be greater than 0' )
  }

  # Determine max number of observations
  n <- max( c(
    length( mu ),
    length( sigma )
  ) )

  # Ensure input vectors are same length
  mu <- rep_len( mu, n )
  sigma <- rep_len( sigma, n )

  # Generate draws
  y <- rnorm( n, mu, sigma )

  return( y )
}

#### 2.3) sim.normal ####
#' @rdname sim.gaussian
#' @export

sim.normal <- function( mu, sigma ) {
  return( sim.gaussian( mu, sigma ) )
}

#### 2.4) sim.log_normal ####
#' Simulate From a Log-Normal Distribution
#'
#' Generate draws from a Log-Normal distribution.
#'
#' @param mu A vector of location parameters, the
#'   mean of the log of the outcome.
#' @param sigma A vector of dispersion parameters (greater than
#'   0).
#'
#' @return A vector of values.
#'
#' @examples
#'
#' # For reproducibility
#' set.seed( 111 )
#' # Generate two draws
#' mu <- c( 0, 0.194 )
#' sigma <- 1
#' round( sim.log_normal( mu, sigma ), 2 )
#'
#' @export

sim.log_normal <- function( mu, sigma ) {

  # Check for valid inputs

  if ( min( sigma ) <= 0 ) {
    stop( 'Argument "sigma" must be greater than 0' )
  }

  # Determine max number of observations
  n <- max( c(
    length( mu ),
    length( sigma )
  ) )

  # Ensure input vectors are same length
  mu <- rep_len( mu, n )
  sigma <- rep_len( sigma, n )

  log_y <- rnorm( n, mu, sigma )
  y <- exp( log_y )

  return( y )
}

#### 2.5) sim.student_t ####
#' Simulate From a Student-T Distribution
#'
#' Generate draws from a Student-T distribution.
#'
#' @param mu A vector of location parameters, the distribution
#'   mean.
#' @param sigma A vector of dispersion parameters (greater than
#'   0).
#' @param nu A vector of parameters for the degrees of freedom
#'   (greater than 2).
#'
#' @return A vector of values.
#'
#' @examples
#'
#' # For reproducibility
#' set.seed( 111 )
#' # Generate two draws
#' mu <- c( 0, 100 )
#' sigma <- c( 1, 15 )
#' nu <- c( 3, 99 )
#' round( sim.student_t( mu, sigma, nu ), 2 )
#'
#' @export

sim.student_t <- function( mu, sigma, nu ) {

  # Check for valid inputs

  if ( min( sigma ) <= 0 ) {
    stop( 'Argument "sigma" must be greater than 0' )
  }

  if ( min( nu ) <= 2 ) {
    stop( 'Argument "nu" must be greater than 2' )
  }

  # Determine max number of observations
  n <- max( c(
    length( mu ),
    length( sigma ),
    length( nu )
  ) )

  # Ensure input vectors are same length
  mu <- rep_len( mu, n )
  sigma <- rep_len( sigma, n )
  nu <- rep_len( nu, n )

  # Generate draws from Student-T distribution
  # with mean of 0, SD of 1, and 'nu' degrees
  # of freedom
  y2 <- rt( n, df = nu )

  # Rescale distribution
  y1 <- y2 / ( sqrt( nu / (nu - 2) ) )
  y <- y1 * sigma + mu

  return( y )
}

#### 3) Multivariate outcomes ####

#### 3.1) sim.multivariate_normal ####
#' Simulate From a Multivariate Normal Distribution
#'
#' Generate N draws from a Multivariate Normal distribution.
#'
#' @param Mu Either (a) a vector of K location parameters,
#'   or (b) an N x K matrix of location parameters.
#' @param Sigma Either (a) a K x K covariance matrix, or (b)
#'   a K x K x N array of covariance matrices.
#' @param Omega Either (a) a K x K correlation matrix, or (b)
#'   a K x K x N array of correlation matrices. The argument
#'   \code{tau} must also be provided.
#' @param tau Either (a) a vector of K standard deviations,
#'   or (b) an N x K matrix of standard deviations. The argument
#'   \code{Omega} must also be provided.
#'
#' @return A N x K matrix of values.
#'
#' @examples
#'
#' # For reproducibility
#' set.seed( 111 )
#'
#' # Generate two draws for 3 correlated variables
#' N <- 2; K <- 3
#' # Location parameters
#' Mu <- matrix( c( 0, 100 ), N, K )
#' # Correlation matrix
#' Omega <- rbind( c( 1, .2, .2 ), c( .2, 1, .4 ), c( .2, .4, 1 ) )
#' Omega <- array( Omega, dim = c( K, K, N ) )
#' # Standard deviations
#' tau <- matrix( c( 1, 15 ), N, K )
#' # Covariance matrix
#' Sigma <- Omega
#' Sigma[,,1] <- diag(tau[1,]) %*% Sigma[,,1] %*% diag(tau[1,])
#' Sigma[,,2] <- diag(tau[2,]) %*% Sigma[,,2] %*% diag(tau[2,])
#'
#' round( sim.multivariate_normal( Mu, Sigma ), 2 )
#' round( sim.multivariate_normal( Mu, Omega = Omega, tau = tau ), 2 )
#'
#' @export

sim.multivariate_normal <- function( Mu, Sigma = NULL,
                                     Omega = NULL, tau = NULL ) {

  # Check for valid inputs

  # If covariance matrix not provided
  if ( is.null( Sigma ) ) {

    # Check that both 'Omega' and 'tau' were provided
    if ( is.null( Omega ) | is.null( tau ) ) {

      stop( 'Must provide either "Sigma" or "Omega" & "tau"' )

      # Close 'Check that both "Omega" and "tau" were provided'
    }

    # Check for valid inputs for "Omega"

    err_msg <- paste0(
      '"Omega" must be a square positive-definite matrix or ',
      'an array of said matrices'
    )

    dmn <- dim( Omega )

    if ( is.null( dmn ) ) stop( err_msg )

    if ( !length( dmn ) %in% c( 2, 3 ) ) stop( err_msg )

    if ( dmn[1] != dmn[2] ) stop( err_msg )

    err_msg <- paste0(
      '"Omega" must be a correlation matrix'
    )
    if ( min( Omega ) < -1 |
         max( Omega ) > 1 ) {
      stop( err_msg )
    }

    if ( length( dmn ) == 2 ) {
      if ( any( diag( Omega ) != 1 ) ) stop( err_msg )
    } else {
      vec <- sapply( 1:dmn[3], function(i) all( diag( Omega[,,i] ) != 1 ) )
      if ( any( vec ) ) stop( err_msg )
    }

    if ( length( dmn ) == 2 ) {
      if ( det( Omega <= 0 ) ) stop( err_msg )
    } else {
      vec <- sapply( 1:dmn[3], function(i) det( Omega[,,i] ) <= 0 )
      if ( any( vec ) ) stop( err_msg )
    }

    K <- dmn[1]

    # Check for valid inputs for "tau"

    err_msg <- paste0(
      '"tau" must be vector of K (or matrix of N x K) ',
      'standard deviations'
    )

    if ( !is.matrix( tau ) ) {
      tau <- rbind( tau )
    }

    if ( ncol( tau ) != K ) {
      stop( err_msg )
    }

    if ( min(tau) <= 0 ) {
      stop( err_msg )
    }

    # Close 'If covariance matrix not provided'
  } else {

    # Check for valid inputs for "Sigma"

    err_msg <- paste0(
      '"Sigma" must be a square positive-definite matrix or ',
      'an array of said matrices'
    )

    dmn <- dim( Sigma )

    if ( is.null( dmn ) ) stop( err_msg )

    if ( !length( dmn ) %in% c( 2, 3 ) ) stop( err_msg )

    if ( dmn[1] != dmn[2] ) stop( err_msg )

    err_msg <- paste0(
      '"Sigma" must be a covariance matrix'
    )

    if ( length( dmn ) == 2 ) {
      if ( det( Sigma ) <= 0 ) stop( err_msg )
    } else {
      vec <- sapply( 1:dmn[3], function(i) det( Sigma[,,i] ) <= 0 )
      if ( any( vec ) ) stop( err_msg )
    }

    K <- dmn[1]

  }

  if ( !is.matrix( Mu ) ) {
    Mu <- rbind( Mu )
  }

  if ( ncol( Mu ) != K ) {
    stop( paste0(
      '"Mu" must be a vector of K (or a N x K matrix) location parameters'
      )
    )
  }

  # Determine max number of observations

  n1 <- nrow( Mu )

  n2 <- 1
  if ( length( dmn ) > 2 ) n2 <- dmn[3]

  if ( is.null( Sigma ) ) {

    n3 <- nrow( tau )

    n <- max( n1, n2, n3 )

  } else {

    n <- max( n1, n2 )

  }

  # Ensure input matrices and arrays are same length

  if ( nrow( Mu ) == 1 ) {
    Mu <- matrix( Mu, n, K, byrow = TRUE )
  } else {
    Mu <- matrix( Mu, n, K )
  }


  if ( is.null( Sigma ) ) {

    if ( nrow( tau ) == 1 ) {
      tau <- matrix( tau, n, K, byrow = TRUE )
    } else {
      tau <- matrix( tau, n, K )
    }

    if ( length( dmn ) == 2 ) {
      Omega_orig <- array( Omega, dim = c( K, K, 1) )
      index <- rep( 1, n )
    } else {
      Omega_orig <- Omega
      index <- rep_len( 1:dmn[3], n )
    }

    Sigma <- array( NA, dim = c( K, K, n ) )
    for ( i in 1:n ) {
      Sigma[,,i] <-
        diag( tau[i,] ) %*%
        Omega_orig[,,index[i]] %*%
        diag( tau[i,] )
    }

  } else {

    if ( length( dmn ) == 2 ) {
      Sigma_orig <- array( Sigma, dim = c( K, K, 1) )
      index <- rep( 1, n )
    } else {
      Sigma_orig <- Sigma
      index <- rep_len( 1:dmn[3], n )
    }

    Sigma <- array( NA, dim = c( K, K, n ) )
    for ( i in 1:n ) {
      Sigma[,,i] <- Sigma_orig[,,index[i]]
    }

  }

  # Generate probabilities from Multivariate Normal distribution

  out <- sapply( 1:n, function(i) {
    MASS::mvrnorm( n = 1, mu = Mu[i,], Sigma = Sigma[,,i] )
  } ) |> t()

  return( out )
}

#### 3.2) sim.lkj ####
#' Simulate From an LKJ Distribution
#'
#' Generate draws from an LKJ distribution. Implementation
#' taken from code in McElreath's (2021) R package 'rethinking'
#' and based on work by Lewandowski et al. (2009).
#'
#' @param eta A vector of N shape parameters (values less than 1
#'   produce correlation matrices with more extreme off-diagonals
#'   closer to -1 or 1, while values greater than 1 produce
#'   correlation matrices with less extreme off-diagonals closer
#'   to 0).
#' @param K An integer equal to 2 or larger, the number of
#'   rows/columns for the correlation matrices.
#'
#' @references
#' * Lewandowski, D., Kurowicka, D., & Joe, H. (2009). Generating
#'     random correlation matrices based on vines and extended
#'     onion method. Journal of Multivariate Analysis, 100 (9),
#'     1989–2001. https://doi.org/10.1016/j.jmva.2009.04.008.
#' * McElreath, R. (2021). rethinking [Computer software].
#'     Retrieved from https://github.com/rmcelreath/rethinking.
#'
#' @returns An array of N correlation matrices of size K x K.
#'
#' @examples
#'
#' # For reproducibility
#' set.seed( 111 )
#' # Generate two draws
#' eta <- c( 1, 4 ); K <- 3
#' round( sim.lkj( eta, K ), 2 )
#'
#' @export

sim.lkj <- function( eta, K ) {

  # Check for valid inputs

  if ( min( eta ) <= 0 ) {
    stop( 'Argument "eta" must be greater than 0' )
  }

  if ( length(K) > 1 ) {
    K <- K[1]
    warning(
      "Dimensions of correlation matrix set using only first value of 'K'"
    )
  }

  if ( K < 2 ) {
    stop( 'Argument "K" must be an integer greater than 1' )
  }

  # Ensure 'K' is an integer
  K <- round( K )

  # Determine max number of observations
  N <- length( eta )

  # Initialize output
  Omega <- array( 0, dim = c( K, K, N ) )

  # Loop over observations
  for ( n in 1:N ) {

    # Generate random correlation matrix
    Omega[,,n] <- extended_onion_method( eta[n], K )

    # Close 'Loop over observations'
  }

  if ( N == 1 ) {
    Omega <- Omega[,,1]
  }

  return( Omega )
}

#### 3.2.1) extended_onion_method ####
# Extended Onion Method for Generating Correlation Matrices
#
# McElreath's (2021) implementation of the extended
# onion method from Lewandowski et al. (2009) for
# generating random correlation matrices
#
# @param 'eta' A shape parameter (values less than 1
#   produce correlation matrices with more extreme off-diagonals
#   closer to -1 or 1, while values greater than 1 produce
#   correlation matrices with less extreme off-diagonals closer
#   to 0).
# @param 'K' An integer equal to 2 or larger, the number of
#   rows/columns for the correlation matrix.
#
# @references
# * Lewandowski, D., Kurowicka, D., & Joe, H. (2009). Generating
#     random correlation matrices based on vines and extended
#     onion method. Journal of Multivariate Analysis, 100 (9),
#     1989–2001. https://doi.org/10.1016/j.jmva.2009.04.008.
# * McElreath, R. (2021). rethinking [Computer software].
#     Retrieved from https://github.com/rmcelreath/rethinking.
#
# @returns A K x K correlation matrix R with density
# proportional to [ det(R) ]^(eta-1) for eta > 1.

extended_onion_method <- function( eta, K ) {

  # Initialize parameter for simulating
  # from the beta distribution
  alpha <- eta + (K - 2)/2

  # Simulate from the beta distribution and
  # rescale draw to be between [-1, 1]
  r12 <- 2 * rbeta(1, alpha, alpha) - 1

  # Initialize output
  R <- matrix(0, K, K)

  # Upper triangular Cholesky factor
  R[1,1] <- 1
  R[1,2] <- r12
  R[2,2] <- sqrt(1 - r12^2)

  # If larger than bivariate
  if (K > 2) {

    # Loop over additional correlations
    for (m in 2:(K - 1) ) {

      alpha <- alpha - 0.5
      y <- rbeta(1, m / 2, alpha)

      # Draw uniformly on a hypersphere
      z <- rnorm(m, 0, 1)
      z <- z / sqrt(crossprod(z)[1])

      R[1:m,m+1] <- sqrt(y) * z
      R[m+1,m+1] <- sqrt(1 - y)

      # Close 'Loop over additional correlations'
    }

    # Close 'If larger than bivariate'
  }

  return( crossprod(R) )
}

#### 4) Miscellaneous ####

# sim.mixture( prob, ... ) {
#
#   args <- list(...)
#   K <- length( args )
#
#   for ( k in 1:K ) {
#
#     if ( k == 1 ) {
#
#       M <- matrix( NA, length( args[[k]] ), K )
#
#     }
#
#     M[,k] <- args[[k]]
#   }
#
# }

