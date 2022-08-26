# Distribution functions for data generation
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you
# have any questions or comments
# Last updated 2022-08-21

# Table of contents
# 1) transformation
# 2) linear_combo.pop_level
# 3) linear_combo.cluster

#### 1) transformation ####
#' Transform a Vector of Values
#'
#' Function to transform a vector of values.
#'
#' @param values A vector of numeric values.
#' @param type A character string, the type of
#'   transformation function to use, either...
#'   * The \code{'identity'} function;
#'   * The \code{'exponential'} function;
#'   * The \code{'logistic'} function.
#'   * The \code{'logit'} function.
#'   * The \code{'linear'} function.
#' @param inverse Logical; if \code{TRUE}
#'   returns the results of the inverse function
#'   (e.g., applies the logit function instead of
#'   the logistic function).
#'
#' @details
#' This function carries out a variety of standard
#' transformations of values, and can be used to
#' implement link functions for a variety of
#' statistical models.
#'
#' # Identity function
#'
#' Implements the function \eqn{ f(x) = x } or
#' \eqn{ f^{-1}(y) = y. }.
#'
#' # Exponential function
#'
#' Implements the function \eqn{ f(x) = e^x } or
#' \eqn{ f^{-1}(y) = log(y). }.
#'
#' # Logistic function
#'
#' Implements the function \eqn{ f(x) = \frac{1}{1 + e^{-x}} } or
#' \eqn{ f^{-1}(y) = log( y / (1-y) ). }.
#'
#' # Logit function
#'
#' Implements the function \eqn{ f(x) = log( x / (1-x) ) } or
#' \eqn{ f^{-1}(y) = \frac{1}{1 + e^{-y}}. }.
#'
#' # Linear function
#'
#' Implements the function \eqn{ f(x) = a + b(x) } or
#' \eqn{ f^{-1}(y) = \frac{y - a}{b}. }.
#'
#' @returns A vector of transformed values.
#'
#' @examples
#'
#' # Exponential function
#' round( transformation( -1.204, type = 'exponential' ), 3 )
#' round( transformation(  0.300, type = 'exponential', inverse = TRUE ), 3 )
#'
#' # Logistic function
#' transformation( 0.0, type = 'logistic' )
#' transformation( 0.5, type = 'logistic', inverse = TRUE )
#'
#' # Linear function
#' transformation( 115, type = 'linear', parameters = c(-100/15, 1/15) )
#' transformation(
#'   -1, type = 'linear', parameters = c(-100/15, 1/15), inverse = TRUE
#' )
#'
#' @export

transformation <- function( values, type = 'identity',
                            inverse = FALSE,
                            parameters = NULL ) {

  out <- NULL

  types <- list(
    identity = 'identity',
    exponential = 'exponential',
    logistic = 'logistic',
    logit = 'logit',
    linear = 'linear'
  )

  possible_types <- paste(
    paste0( "  '", names( types ), "'" ), collapse = '\n'
  )

  # Display options for type of transformation
  if ( is.null( values ) ) {

    message(
      paste0(
        "Options for argument 'type' include:\n",
        possible_types
      )
    )

    return( invisible(NULL) )

    # Close 'Display options for type of transformation'
  }

  # Transformation: identity
  if ( type %in% types$identity ) {

    out <- values

    # Close 'Transformation: identity'
  }

  # Transformation: exponential
  if ( type %in% types$exponential ) {

    # Apply transformation
    if ( !inverse ) {

      out <- exp( values )

      # Close 'Apply transformation'
    }

    # Apply inverse of transformation
    if ( inverse ) {

      out <- log( values )

      # Close 'Apply inverse of transformation'
    }

    # Close 'Transformation: exponential'
  }

  # Transformation: logistic
  if ( type %in% types$logistic ) {

    # Apply transformation
    if ( !inverse ) {

      out <- 1 / ( 1 + exp( -values ) )

      # Close 'Apply transformation'
    }

    # Apply inverse of transformation
    if ( inverse ) {

      out <- log( values / ( 1 - values ) )

      # Close 'Apply inverse of transformation'
    }

    # Close 'Transformation: logistic'
  }


  # Transformation: logit
  if ( type %in% types$logit ) {

    # Apply transformation
    if ( !inverse ) {

      out <- log( values / ( 1 - values ) )

      # Close 'Apply transformation'
    }

    # Apply inverse of transformation
    if ( inverse ) {

      out <- 1 / ( 1 + exp( -values ) )

      # Close 'Apply inverse of transformation'
    }

    # Close 'Transformation: logit'
  }


  # Transformation: linear
  if ( type %in% types$linear ) {

    if ( is.null( parameters) ) {
      stop( paste0(
        "Linear transformation requires a vector of two values ",
        "to be passed to the argument 'parameters"
      ) )
    }

    a <- parameters[1]
    b <- parameters[2]

    # Apply transformation
    if ( !inverse ) {

      out <- a + b * values

      # Close 'Apply transformation'
    }

    # Apply inverse of transformation
    if ( inverse ) {

      out <- ( values - a ) / b

      # Close 'Apply inverse of transformation'
    }

    # Close 'Transformation: linear'
  }

  # If type of transformation misspecified
  if ( is.null( out ) ) {

    stop(
      paste0(
        'Argument "type" must be:\n',
        possible_types
      )
    )

    # Close 'If type of transformation misspecified'
  }

  return( out )
}

#### 2) linear_combo.population ####
#' Linear Combination of Design Matrix and Coefficients
#'
#' Computes the linear combination of a design
#' matrix \code{X} and a vector of regression coefficients
#' \code{betas} (equivalent to \code{X %*% betas}).
#' Applicable for population-level coefficients
#' (often called fixed effects).
#'
#' @param X A N x K design matrix.
#' @param betas A vector of K regression coefficients.
#'
#' @returns A vector of N values.
#'
#' @examples
#' X <- rbind(
#'   c( 1, 0, 0 ),
#'   c( 1, 1, 0 ),
#'   c( 1, 0, 1 )
#' )
#' linear_combo.pop_level( X, c( 1, 1, 2) )
#'
#' @export

linear_combo.population <- function( X, betas ) {

  X <- as.matrix( X )
  N <- nrow( X )
  K <- ncol( X )

  if ( length( betas ) != K ) {
    stop(
      "Argument 'betas' must have same number as columns as argument 'X'"
    )
  }

  betas <- matrix( betas, K, 1, byrow = FALSE )

  alpha <- X %*% betas

  return( alpha[,1] )
}

#### 3) linear_combo.cluster ####
#' Linear Combination of Design Matrix and Coefficients for Clusters
#'
#' Given S clusters, N observations, and K variables,
#' computes the linear combination of K coefficients and
#' the relevant subset of the design matrix X for a given
#' cluster. Applicable for cluster-level coefficients
#' (often called random effects).
#'
#' @param etas A S x K matrix of coefficients.
#' @param index A vector of N integers from 1 to S used to
#'   expand the matrix \code{etas} to match in number of
#'   rows to the design matrix \code{X}.
#' @param X A N x K design matrix.
#'
#' @returns A vector of N values.
#'
#' @examples
#' index <- c( 1, 1, 2, 2, 3, 3 )
#' X <- cbind( 1, rep( 0:1, 3) )
#' etas <- MASS::mvrnorm( 3, rep( 0, 2 ), diag(2) )
#' etas |> linear_combo.cluster( index, X )
#'
#' @export

linear_combo.cluster <- function( etas, index, X ) {

  X <- as.matrix( X )
  N <- nrow( X )
  K <- ncol( X )

  if ( length( index ) != N ) {
    stop(
      "Argument 'index' must be same length as number of rows in 'X'"
    )
  }

  if ( ncol(etas) != K ) {
    stop(
      "Argument 'etas' must have same number as columns as argument 'X'"
    )
  }

  if ( K == 1 ) {
    indexed_etas <- cbind( etas[index] )
  } else {
    indexed_etas <- etas[index,]
  }

  out <- rowSums( X * indexed_etas )

  return( out )
}



