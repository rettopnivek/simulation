# Distribution functions for data generation
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you
# have any questions or comments
# Last updated 2022-08-12

# Table of contents
# 1) transformation
# 2) linear_combo.pop_level

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
#' Implements the function \deqn{ f(x) = x } or
#' \deqn{ f^{-1}(y) = y. }.
#'
#' # Exponential function
#'
#' Implements the function \deqn{ f(x) = e^x } or
#' \deqn{ f^{-1}(y) = log(y). }.
#'
#' # Logistic function
#'
#' Implements the function \deqn{ f(x) = frac{1}{1 + e^{-x}} } or
#' \deqn{ f^{-1}(y) = log( y / (1-y) ). }.
#'
#' # Logit function
#'
#' Implements the function \deqn{ f(x) = log( x / (1-x) ) } or
#' \deqn{ f^{-1}(y) = frac{1}{1 + e^{-y}}. }.
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
#' @export

transformation <- function( values, type, inverse = FALSE ) {

  out <- NULL

  all_types <- c(
    identity = 'identity',
    exponential = 'exponential',
    logistic = 'logistic',
    logit = 'logit'
  )

  if ( type == all_types['identity'] ) {

    out <- values

  }

  if ( type == all_types['exponential'] ) {

    if ( !inverse ) {

      out <- exp( values )

    } else {

      out <- log( values )

    }

  }

  if ( type == all_types['logistic'] ) {

    if ( !inverse ) {

      out <- 1 / ( 1 + exp( -values ) )

    } else {

      out <- log( values / ( 1 - values ) )

    }

  }


  if ( type == all_types['logit'] ) {

    if ( !inverse ) {

      out <- log( values / ( 1 - values ) )

    } else {

      out <- 1 / ( 1 + exp( -values ) )

    }

  }

  if ( is.null( out ) ) {
    stop(
      paste0(
        'Argument "type" must be:\n',
        paste( paste0( '   ', all_types ), collapse = '\n' )
      )
    )
  }

  return( out )
}

#### 2) linear_combo.pop_level ####
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

linear_combo.pop_level <- function( X, betas ) {

  X <- as.matrix( X )
  N <- nrow( X )
  K <- ncol( X )

  betas <- matrix( betas, K, 1, byrow = FALSE )

  alpha <- X %*% betas

  return( alpha[,1] )
}

