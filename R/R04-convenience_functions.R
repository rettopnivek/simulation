# Convenience/helper functions
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you
# have any questions or comments
# Last updated 2022-08-24

# Table of contents
# 1) Functions for calculations
#   1.1) calc.beta_binomial
# 2) Functions for matrix operations
#   2.1) matrix.identity
#   2.2) matrix.corr
# 3) run_in_parallel

#### 1) Functions for calculations ####

#### 1.1) calc.distribution_char ####
#' Calculations for Distribution Characteristics
#'
#' Function to carry out a variety of useful calculations
#' to describe various statistical distributions
#'
#' @param parameters A vector of numeric values, the parameter
#'   values for a distribution (e.g., the mean and standard
#'   deviation for a normal distribution).
#' @param distribution A character string, the statistical
#'   distribution for which to calculate characteristics,
#'   either...
#'   * \code{'bernoulli'}
#'   * \code{'binomial'}
#'   * \code{'beta-binomial'}
#'   * \code{'beta'}
#'   * \code{'gaussian'} or \code{'normal'}
#'   * \code{'log_normal'}
#' @param statistic A character string, the type of
#'   characteristic to calculate, either...
#'   * \code{'mean'}
#'   * \code{'variance'}
#'   * \code{'SD'}
#'   * \code{'skew'}
#'   * \code{'kurtosis'}
#'
#' @return A numeric value, results of the specified
#' calculation.
#'
#' @examples
#'
#' # FORTHCOMING
#'
#' @export

calc.distribution_char <- function(
    parameters,
    distribution,
    statistic = 'mean' ) {

  #### 1.1.1) Setup ####

  distributions <- list(
    bernoulli = c(
      'Bernoulli', 'bernoulli'
    ),
    binomial = c(
      'Binomial', 'binomial'
    ),
    beta = c(
      'Beta', 'beta'
    ),
    log_normal = c(
      'Log-Normal', 'log-normal',
      'Lognormal', 'lognormal',
      'log_normal'
    )
  )

  statistics <- list(
    mean = c( 'Mean', 'mean', 'M', 'm' ),
    median = c( 'Median', 'median', 'Md', 'md' ),
    mode = c( 'Mode', 'mode', 'Mo', 'mo' ),
    variance = c( 'Variance', 'variance', 'V', 'v' ),
    SD = c( 'Standard deviation', 'standard deviation', 'SD', 'sd' ),
    skew = c( 'Skewness', 'skewness', 'Skew', 'skew' ),
    kurtosis = c( 'Ex. kurtosis', 'ex. kurtosis', 'Kurtosis', 'kurtosis' )
  )

  out <- NA

  stat_to_compute <- sapply( seq_along( statistics ), function(i) {
    statistic %in% statistics[[i]]
  } )
  names( out ) <- names( statistics )[ stat_to_compute ]

  #### 1.1.2) Distributions ####

  #### Bernoulli distribution ####
  if ( distribution %in% distributions$bernoulli ) {

    p <- parameters; q <- 1 - p

    # Compute mean
    if ( statistic %in% statistics$mean ) {

      out[1] <- p

      # Close 'Compute mean'
    }

    # Compute variance
    if ( statistic %in% statistics$variance ) {

      out[1] <- p*q

      # Close 'variance'
    }

    # Compute standard deviation
    if ( statistic %in% statistics$SD ) {

      out[1] <- sqrt( p*q )

      # Close 'Compute standard deviation'
    }

    # Compute skewness
    if ( statistic %in% statistics$skew ) {

      out[1] <- (q - p) / sqrt( p*q )

      # Close 'Compute skewness'
    }

    # Compute kurtosis
    if ( statistic %in% statistics$kurtosis ) {

      out[1] <- (1 - 6*p*q) / (p*q)

      # Close 'Compute kurtosis'
    }

    if ( is.na( out ) ) {
      stop(
        paste0(
          "Cannot calculate ", statistic,
          " for the Bernoulli distribution"
        )
      )
    }

    # Close 'Bernoulli distribution'
  }

  #### Binomial distribution ####
  if ( distribution %in% distributions$binomial ) {

    p <- parameters[1]; q <- 1 - p; n <- parameters[2]

    # Compute mean
    if ( statistic %in% statistics$mean ) {

      out[1] <- n*p

      # Close 'Compute mean'
    }

    # Compute variance
    if ( statistic %in% statistics$variance ) {

      out[1] <- n*p*q

      # Close 'variance'
    }

    # Compute standard deviation
    if ( statistic %in% statistics$SD ) {

      out[1] <- sqrt( n*p*q )

      # Close 'Compute standard deviation'
    }

    # Compute skewness
    if ( statistic %in% statistics$skew ) {

      out[1] <- (q - p) / sqrt( n*p*q )

      # Close 'Compute skewness'
    }

    # Compute kurtosis
    if ( statistic %in% statistics$kurtosis ) {

      out[1] <- (1 - 6*p*q) / (n*p*q)

      # Close 'Compute kurtosis'
    }

    if ( is.na( out ) ) {
      stop(
        paste0(
          "Cannot calculate ", statistic,
          " for the Binomial distribution"
        )
      )
    }

    # Close 'Binomial distribution'
  }

  #### Beta distribution ####
  if ( distribution %in% distributions$beta ) {

    mu <- parameters[1]; phi <- parameters[2]

    shape_alpha <- mu * phi
    shape_beta <- (1 - mu) * phi
    a_x_b <- shape_alpha * shape_beta
    a_plus_b <- shape_alpha + shape_beta
    b_minus_a <- shape_beta - shape_alpha
    a_minus_b <- shape_alpha - shape_beta

    # Compute mean
    if ( statistic %in% statistics$mean ) {

      out[1] <- mu

      # Close 'Compute mean'
    }

    # Compute variance
    if ( statistic %in% statistics$variance ) {

      num <- a_x_b
      denom <- ( a_plus_b^2 ) * ( a_plus_b + 1 )

      out[1] <- num / denom

      # Close 'variance'
    }

    # Compute standard deviation
    if ( statistic %in% statistics$SD ) {

      num <- a_x_b
      denom <- ( a_plus_b^2 ) * ( a_plus_b + 1 )

      out[1] <- sqrt( num / denom )

      # Close 'Compute standard deviation'
    }

    # Compute skewness
    if ( statistic %in% statistics$skew ) {

      num <- 2*b_minus_a*sqrt( a_plus_b + 1 )
      denom <- ( a_plus_b + 2 ) * sqrt( a_x_b )

      out[1] <- num / denom

      # Close 'Compute skewness'
    }

    # Compute kurtosis
    if ( statistic %in% statistics$skew ) {

      part_1 <- ( a_minus_b^2 ) * ( a_plus_b + 1 )
      part_2 <- a_x_b * ( a_plus_b + 2 )
      part_3 <- a_x_b * ( a_plus_b + 2 ) * ( a_plus_b + 3 )

      out[1] <- 6*( part_1 - part_2 ) / part_3

      # Close 'Compute kurtosis'
    }

    # Close 'Beta distribution'
  }

  #### Log-Normal distribution ####
  if ( distribution %in% distributions$log_normal ) {

    mu <- parameters[1]; sigma <- parameters[2]; sigma2 <- sigma^2

    # Compute mean
    if ( statistic %in% statistics$mean ) {

      out[1] <- exp( mu + sigma2/2 )

      # Close 'Compute mean'
    }

    # Compute median
    if ( statistic %in% statistics$median ) {

      out[1] <- exp( mu )

      # Close 'Compute median'
    }

    # Compute variance
    if ( statistic %in% statistics$variance ) {

      out[1] <- ( exp( sigma2 ) - 1 ) * exp( 2 * mu + sigma2 )

      # Close 'variance'
    }

    # Compute standard deviation
    if ( statistic %in% statistics$SD ) {

      out[1] <- sqrt( ( exp( sigma2 ) - 1 ) * exp( 2 * mu + sigma2 ) )

      # Close 'Compute standard deviation'
    }

    # Compute skewness
    if ( statistic %in% statistics$skew ) {

      out[1] <- ( exp( sigma2 ) + 2 ) * sqrt( exp( sigma2) - 1 )

      # Close 'Compute skewness'
    }

    # Compute kurtosis
    if ( statistic %in% statistics$kurtosis ) {

      out[1] <-
        exp( 4 * sigma2 ) +
        2 * exp( 3 * sigma2 ) +
        3 * exp( 2 * sigma2 ) - 6

      # Close 'Compute kurtosis'
    }

    if ( is.na( out ) ) {
      stop(
        paste0(
          "Cannot calculate ", statistic,
          " for the Log-Normal distribution"
        )
      )
    }

    # Close 'Log-Normal distribution'
  }

  return( out )
}


#### 1.1) calc.beta_binomial ####
#' Calculations for Beta Distribution and Beta-Binomial Model
#'
#' Function to carry out a variety of useful calculations
#' for the beta distribution and beta-binomial model.
#'
#' @param inputs A vector of numeric values, either...
#'   * The two shape parameters (\code{'mean'},
#'     \code{'variance'}, and \code{'mu & phi'}).
#'   * The mean and dispersion parameters (\code{'alpha & beta'}).
#'   * The number of counts, max possible number of counts,
#'     prior for the alpha parameter, prior for the beta
#'     parameter, and width of the credible interval.
#' @param type A character string, the type of calculation
#'   to carry out, either...
#'   * \code{'mean'} (compute mean of beta distribution).
#'   * \code{'variance'} (compute variance of beta distribution).
#'   * \code{'mu & phi'} (convert shape parameters to mean and
#'     dispersion parameters).
#'   * \code{'alpha & beta'} (convert mean and dispersion
#'     parameters to mean and dispersion parameters).
#'   * \code{'credible interval'} (compute credible interval
#'     for beta-binomial model).
#'
#' @return A vector of values, results of the specified
#' calculation.
#'
#' @examples
#'
#' # Mean and variance of beta distribution
#' calc.beta_binomial( c( 1, 4 ), 'mean' )
#' calc.beta_binomial( c( 1, 4 ), 'variance' )
#'
#' # Convert shape parameters to mean/dispersion parameters
#' calc.beta_binomial( c( 1, 4 ), 'mu & phi' )
#'
#' # Convert mean/dispersion parameters to shape parameters
#' calc.beta_binomial( c( 0.2, 5 ), 'alpha & beta' )
#'
#' # Credible interval for beta-binomial model
#' set.seed( 111 )
#' x <- sim.beta_binomial( .2, 5, 100 )
#' # Default priors ( alpha = .5, beta = .5 )
#' calc.beta_binomial( c( x, 100 ), 'credible interval' )
#' # Custom priors
#' calc.beta_binomial( c( x, 100, 1, 1 ), 'credible interval' )
#' # Different interval width
#' calc.beta_binomial( c( x, 100, 1, 1, .68 ), 'credible interval' )
#'
#' @export

calc.beta_binomial <- function( inputs, type = 'mean' ) {

  # Possible options for argument 'type'
  types <- list(
    mean = c(
      'Mean', 'mean', 'M', 'm'
    ),
    variance = c(
      'Variance', 'variance', 'V', 'v'
    ),
    mu_phi = c(
      'mu and phi', 'mu & phi', 'MP', 'mp', 'mu_phi'
    ),
    alpha_beta = c(
      'alpha and beta', 'alpha & beta', 'AB', 'ab', 'alpha_beta'
    ),
    ci = c(
      'Credible interval', 'credible interval', 'CI', 'ci'
    )
  )

  possible_types <- sapply(
    seq_along( types), function(i) types[[i]][2]
  )
  possible_types <- paste(
    paste0( "  '", possible_types, "'" ), collapse = '\n'
  )

  # If no input is provided
  if ( is.null( inputs ) ) {

    message(
      paste0(
        "Options for argument 'type' include:\n",
        possible_types
      )
    )
    return( invisible(NULL) )

    # Close 'If no input is provided'
  }

  # Compute mean of beta distribution
  if ( type %in% types$mean ) {

    out <- inputs[1] / sum( inputs )

    return( out )

    # Close 'Compute mean of beta-binomial distribution'
  }

  # Compute variance of beta distribution
  if ( type %in% types$variance ) {

    num <- prod( inputs )

    denom <- ( sum( inputs )^2 ) * ( sum( inputs ) + 1 )

    out <- num / denom

    return( out )

    # Close 'Compute variance of beta distribution'
  }

  # Convert shape parameters to mean/dispersion parameters
  if ( type %in% types$mu_phi ) {

    mu <- inputs[1] / sum( inputs )
    phi <- inputs[1] / mu

    out <- c( mu = mu, phi = phi )

    return( out )

    # Close 'Convert shape parameters to mean/dispersion parameters'
  }

  # Convert mean/dispersion parameters to shape parameters
  if ( type %in% types$alpha_beta ) {

    mu <- inputs[1] / sum( inputs )
    phi <- inputs[1] / mu

    out <- c(
      alpha = inputs[1] * inputs[2],
      beta = (1 - inputs[1])*inputs[2]
    )

    return( out )

    # Close 'Convert mean/dispersion parameters to shape parameters'
  }

  # Compute credible interval
  if ( type %in% types$ci ) {

    di <- c( 1, 1, .5, .5, .95 )

    di[ seq_along( inputs) ] <- inputs

    alpha_post <- di[3] + di[1]
    beta_post <- di[4] + di[2] - di[1]

    out <- qbeta( .5 + c(-1,1)*di[5]/2, alpha_post, beta_post )

    return( out )

    # Close 'Compute credible interval'
  }

  stop(
    "Check argument 'type' - try 'calc.beta_binomial( NULL )' for options"
  )

}

#### 1.2) calc.sigma ####
#' Calculations for Variance Terms for Linear Regression
#'
#' Assorted calculations for variance terms in the
#' standard linear regression model (i.e., the
#' residual standard deviation, the standard error
#' for coefficients, etc.).
#'
#' @param inputs Either...
#'   * A vector with the coefficient for the predictor,
#'     the variance of the predictor, and the desired
#'     R-squared value.
#'   * A list with the vector of coefficients, the
#'     covariance matrix for the predictors, and
#'     the desired R-squared value.
#'   * A vector with the total number of observations,
#'     the variance for the outcome, the
#'     variance for the given predictor,
#'     and the R-squared value.
#' @param type A character string, either...
#'   * \code{'residual'} (compute standard deviation for
#'     the residuals and outcome).
#'   * \code{'SE'} (compute the standard error for a
#'     regression coefficient).
#'
#' @return A vector of values.
#'
#' @examples
#'
#' # Given a correlation of .5 and standardized
#' # scores compute residual SD
#' calc.sigma( c( beta = .5, var_x = 1, R2 = sqrt(.5) ), 'residual' )
#' calc.sigma( c( n = 30, var_y = 1, var_x = 1, R2 = .1 ), 'SE' )
#'
#' # Multiple correlated predictors
#'
#' # Create covariance matrix
#' Sigma <- matrix.corr( c( .2, .2, .4 ), 'initialize' )
#'
#' # Compute standard deviation for residuals and outcome
#' calc.sigma(
#'   list( beta = c( .2, -.3, .25 ), Sigma = Sigma, R2 = .35 ), 'residual'
#' )
#'
#' @export

calc.sigma <- function( inputs, type = 'residual' ) {

  # Possible options for argument 'type'
  types <- list(
    residual = c(
      'Residual', 'residual', 'resid',
      'Epsilon', 'epsilon'
    ),
    SE = c(
      'Sampling distribution for beta',
      'Beta sampling distribution',
      'Sampling distribution',
      'sampling distribution',
      'Sampling', 'sampling',
      'Standard error',
      'standard error',
      'SE'
    )
  )

  possible_types <- paste(
    paste0( "  '", names(types), "'" ), collapse = '\n'
  )

  # If no input is provided
  if ( is.null( inputs ) ) {

    message(
      paste0(
        "Options for argument 'type' include:\n",
        possible_types
      )
    )
    return( invisible(NULL) )

    # Close 'If no input is provided'
  }

  # Residual and outcome standard deviation
  if ( type %in% types$residual ) {

    # Simple linear regression
    if ( is.numeric( inputs ) ) {

      beta_sq <- inputs[1]^2
      var_x <- inputs[2]
      R2 <- inputs[3]

      sigma <- sqrt( ( beta_sq * var_x * (1 - R2) ) / R2 )
      sigma_Y <- sqrt( (sigma^2) / ( 1 - R2 ) )

      out <- c( sigma, sigma_Y )
      names( out ) <- c( 'residual', 'Y' )

      return( out )

      # Close 'Simple linear regression'
    }

    # Multiple linear regression
    if ( is.list( inputs ) ) {

      betas <- inputs[[1]]
      Sigma_X <- inputs[[2]]
      R2 <- inputs[[3]]

      W <- diag( betas ) %*% Sigma_X %*% diag( betas )

      var_e <- (
        ( sum( diag( W ) ) + 2*sum( W[lower.tri(W)] ) ) * (1 - R2 )
      ) / R2
      sigma <- sqrt( var_e )

      sigma_Y <- sqrt( var_e / ( 1 - R2 ) )

      out <- c( sigma, sigma_Y )
      names( out ) <- c( 'residual', 'Y' )

      return( out )

      # Close 'Multiple linear regression'
    }

    # Close 'Residual and outcome standard deviation'
  }

  # Standard error of the coefficient
  if ( type %in% types$SE ) {

    n <- inputs[1]
    var_y <- inputs[2]
    var_x <- inputs[3]
    R2 <- inputs[4]

    sigma_b1 <- (1/(n-2)) * (var_y/var_x) * (1 - R2^2)
    names( sigma_b1 ) <- 'SE'

    return( sigma_b1 )

    # Close 'Standard error of the coefficient'
  }

  stop(
    "Check argument 'type' - try 'calc.sigma( NULL )' for options"
  )

}

# calc.log_normal <- function() {}

#### 2) Functions for matrix operations ####

#### 2.1) matrix.identity ####
#' Create an Identify Matrix
#'
#' Function to create an identity matrix
#' (square matrix with ones on diagonals,
#' zeros on off-diagonals) of a specified
#' size.
#'
#' @param K Number of elements on the diagonal.
#'
#' @return An identity matrix.
#'
#' @examples
#'
#' matrix.identity( 4 )
#'
#' @export

matrix.identity <- function( K ) {

  return( diag( K ) )

}

#### 2.2) matrix.corr ####
#' Calculations for Correlation Matrices
#'
#' Function to carry out a variety of useful calculations
#' for correlation matrices.
#'
#' @param inputs Either...
#'   * A vector of correlations for the lower triangle of the
#'     correlation matrix (by column).
#'   * A list with (a) a K x K correlation matrix, and (b) a
#'     vector of K standard deviations.
#' @param type A character string, the type of calculation
#'   to carry out, either...
#'   * \code{'initialize'} (given a vector of correlations creates
#'     a correlation matrix).
#'   * \code{'cor2cov'} (given a list with a correlation matrix
#'     and vector of standard deviations computes the corresponding
#'     covariance matrix).
#'
#' @returns A matrix, either a correlation or covariance matrix.
#'
#' @examples
#'
#' Omega <- matrix.corr( c( .2, .2, .4 ), 'init' )
#'
#' Sigma <- matrix.corr( list( Omega, c( 15, 2, 0.5 ) ), 'cor2cov' )
#' Sigma
#' cov2cor( Sigma )
#'
#' @export

matrix.corr <- function( inputs, type ) {

  types <- list(
    initialize = c(
      'Initialize', 'initialize',
      'Initial', 'initial', 'init'
    ),
    cor2Cov = c(
      'Correlation to covariance',
      'correlation to covariance',
      'cor to cov',
      'Cor2cov', 'cor2cov'
    )
  )

  possible_types <- paste(
    paste0( "  '", names( types ), "'" ), collapse = '\n'
  )

  if ( is.null( inputs) ) {
    message( paste0(
      "Possible types are:\n",
      possible_types
    ) )

  }

  if ( type %in% types$initialize ) {

    total <- length( inputs )

    K <- 0.5 * ( sqrt( 8 * total + 1 ) + 1 )

    Omega <- diag( K )

    inc <- 1
    for ( k in 2:K ) {

      index <- inc:( (inc - 1) + (k-1) )

      Omega[ k, 1:(k-1) ] <- inputs[ index ]
      Omega[ 1:(k-1), k ] <- inputs[ index ]

      inc <- max( index ) + 1

    }

    is_pos_def <- det( Omega ) > 0
    btw_pm1 <- min( Omega ) >= -1 & max( Omega ) <= 1
    diag_all_1 <- all( diag( Omega ) == 1 )

    return( Omega )
  }

  if ( type %in% types$cor2Cov ) {

    Sigma <- diag( inputs[[2]] ) %*% inputs[[1]] %*% diag( inputs[[2]] )

    return( Sigma )
  }

  stop(
    paste0(
      'Argument "type" must be:\n',
      possible_types
    )
  )

}

#### 3) run_in_parallel ####
#' Run a Function in Parallel Across Multiple Cores
#'
#' Function to run a user-specified function in
#' parallel across multiple cores via the
#' \link[parallel]{parLapply} function. The
#' function is designed to streamline the steps
#' for parallel processing on Windows machines.
#'
#' @param func_for_iter A function that takes as
#'   a first argument an integer, the current
#'   iteration. Can take additional arguments
#'   (make sure to export relevant variables).
#' @param to_export A character vector, the
#'   variables to export to clusters.
#' @param packages_to_load A character vector,
#'   the packages to load on each cluster.
#' @param n_cores The number of cores to
#'   use for parallel processing (default is
#'   ~75% of available number of cores).
#' @param desired_iter The desired number of
#'   iterations (if not a multiple of \code{n_cores}
#'   will be rounded up internally to the nearest
#'   multiple).
#' @param status Logical; if \code{TRUE} indicates
#'   the setup, parallel processing, and conclusion
#'   of the function.
#' @param ... Additional arguments for the
#'   \code{func_for_iter} function.
#'
#' @return A list of size ~\code{desired_iter} with
#'   the results from each iteration run.
#'
#' @examples
#'
#' func_for_iter <- function( nr ) {
#'   # Pause system for 1 second
#'   Sys.sleep(1)
#'   return(nr)
#' }
#'
#' start_time <- Sys.time()
#' results <- lapply( 1:8, func_for_iter )
#' run_time <- Sys.time() - start_time
#' run_time
#'
#'
#' start_time <- Sys.time()
#' results <- run_in_parallel(
#'   func_for_iter, n_cores = 4, desired_iter = 8
#' )
#' run_time <- Sys.time() - start_time
#' run_time
#'
#' @export

run_in_parallel <- function( func_for_iter,
                             to_export = '',
                             packages_to_load = NULL,
                             n_cores = NULL,
                             desired_iter = 1000,
                             status = FALSE,
                             ... ) {

  if ( status ) message( 'Start' )

  start_time <- Sys.time()

  # Check if 'parallel' package is installed
  all_pck <- installed.packages()
  if ( 'parallel' %in% all_pck[,1] ) {
    library( parallel )
  } else {
    stop( "Requires the 'parallel' package" )
  }

  if ( status ) message( '  Specify cores' )

  if ( is.null( n_cores ) ) {
    n_cores <- detectCores()
    n_cores <- round( n_cores * .75 )
  }

  if ( n_cores <= 1 ) {
    stop( 'Not enough cores for parallel processing' )
  }

  # Initialize cores
  clust <- makeCluster( n_cores )

  if ( status ) message( '  Export variables and packages' )

  packages_to_load <<- packages_to_load
  func_for_iter <<- func_for_iter

  to_export <- c( to_export, 'func_for_iter', 'packages_to_load' )
  to_export <- to_export[ to_export != '' ]

  # Make sure cores have necessary variables
  clusterExport(
    clust,
    varlist = to_export
  )

  # Load in required packages on each core
  if ( !is.null( packages_to_load ) ) {
    check <- clusterEvalQ( clust, {
      lapply( packages_to_load, library, character.only = TRUE)
    } )
  }

  if ( status ) message( '  Run function in parallel' )

  iter <- ceiling( desired_iter / n_cores )
  NR <- iter * n_cores

  # Run simulations in parallel
  list_of_results <- parLapply(
    clust,
    1:NR,
    function(nr) {
      return( func_for_iter( nr, ... ))
    }
  )

  run_time <- Sys.time() - start_time

  stopCluster( clust )

  rm( packages_to_load, envir = .GlobalEnv )

  if ( status ) {
    message( 'Finished' )
    message( run_time )
  }

  return( list_of_results )
}


