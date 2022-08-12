# Distribution functions for data generation
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you
# have any questions or comments
# Last updated 2022-08-11

# Table of contents
# 1) sim.bernoulli
# 2) sim.binomial
# 3) sim.beta_binomial
# 4) sim.gaussian
# 5) sim.student_t

#### 1) sim.bernoulli ####
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

#### 2) sim.binomial ####
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
  n <- max( length( theta ), length( nu ) )

  # Ensure input vectors are same length
  theta <- rep_len( theta, n )
  nu <- rep_len( nu, n )

  # Generate draws
  y <- rbinom( n, nu, theta )

  return( y )
}

#### 3) sim.beta_binomial ####
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

#### 4) sim.gaussian ####
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

#### 5) sim.student_t ####
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


