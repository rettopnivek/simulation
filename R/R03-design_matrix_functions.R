# Design matrix functions
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you
# have any questions or comments
# Last updated 2022-08-12

# Table of contents
# 1) each
# 2) `%j%`
# 3) pse
# 4) groups
# 5) intercept
# 6) get_X

#### 1) each ####
#' Repeat Each Element in a Vector
#'
#' Convenience function that repeats
#' each element in a vector a fixed
#' number of times.
#'
#' @param x A vector of elements.
#' @param times The number of times
#'   to repeat each element.
#'
#' @returns A vector of repeated elements.
#'
#' @examples
#'
#' each( 1:3, 2 )
#' each( LETTERS[1:3], 3 )
#'
#' @export

each <- function( x, times ) {

  return( rep( x, each = times ) )

}

#### 2) `%j%` ####
#' Operator to Join a Vector to a Data Frame
#'
#' An operator that takes a data frame and a
#' list with (1) a vector of values, and (2)
#' the column name to use, and adds the vector
#' as a new named column to the data frame.
#'
#' @param x A data frame or matrix.
#' @param y A list of (1) a vector matching in
#'   length to the number of rows in the data
#'   frame \code{x}, and (2) a character
#'   string, the column name for the new
#'   column.
#'
#' @returns A data frame with the vector in
#'   \code{y} as a new named column.
#'
#' @examples
#' x <- data.frame(
#'   V1 = 1:3
#' )
#' x %j% list( 4:6, 'V2' )
#'
#' @export

`%j%` <- function(x, y) {

  # Check for valid inputs

  # If RHS is a vector
  if ( is.vector( x ) ) {

    # Convert vector to matrix
    x <- cbind( x = x )

    # Close 'If RHS is a vector'
  } else {

    # If RHS is not a matrix or data frame
    if( !is.matrix(x) & !is.data.frame(x) ) {

      stop( 'Object on RHS must be data frame or matrix' )

      # Close 'If RHS is not a matrix or data frame'
    }

    # Close else for 'If RHS is a vector'
  }

  # Number of rows
  n <- nrow( x )

  err_msg <- paste0(
    'Object on LHS must be list with (1) a vector and ',
    '(2) a character string for a column name'
  )

  # If LHS is not a list
  if ( !is.list(y) ) {

    stop( err_msg )

    # Close 'If LHS is not a list'
  }

  # If not a vector and character string
  if ( !is.vector( y[[1]] ) | !is.character( y[[2]] ) ) {

    stop( err_msg )

    # Close 'If not a vector and character string'
  }

  # If number of elements does not match number of rows
  if ( length( y[[1]] ) != n ) {

    stop( 'Length of vector must match number of rows' )

    # Close 'If number of elements does not match number of rows'
  }

  out <- cbind( x, y[[1]] )

  clm <- colnames( out )
  clm[ length(clm) ] <- y[[2]]
  colnames( out ) <- clm

  return( out )
}


#### 3) pse ####
#' Initialize Data Set for Persons x Sessions x Events
#'
#' Generates a data frame with 3 columns for (1) persons
#' (i.e., participants), (2) sessions (i.e., study visits)
#' and (3) events (i.e., repeated measurements collected
#' during a session). Useful, for example, for quickly
#' creating a data frame for simulated data from a
#' longitudinal clinical trial. Also initializes a
#' column for an outcome variable.
#'
#' @param Np Number of individuals or participants.
#' @param Ns Number of sessions.
#' @param Ne Number of events.
#' @param labels A character vector giving up to
#'   3 alternative column names.
#'
#' @return A data frame with \code{Np} x \code{Ns} x
#' \code{Ne} rows. Column contents are integer vectors
#' starting from 1. Includes an additional column
#' \code{'Y'} for an outcome variable initialized
#' to \code{NA}.
#'
#' @examples
#'
#' pse( Np = 2, Ns = 2, Ne = 2 )
#'
#' @export

pse <- function( Np, Ns = 1, Ne = 1, labels = NULL ) {

  out <- data.frame(
    Index.Persons = (1:Np) |> each( Ns * Ne ),
    Index.Sessions = (1:Ns) |> each( Ne ) |> rep( Np ),
    Index.Events = (1:Ne) |> rep( Np * Ns )
  )

  if ( is.null( labels ) ) {
    labels <- colnames( out )
  }

  if ( length( labels ) < ncol( out) ) {

    L <- length( labels )

    labels <- c(
      labels,
      colnames( out )[1:3 > L]
    )

  }

  colnames( out ) <- labels

  out$Y <- NA

  return( out )
}

#### 4) groups ####
#' Split Index Variable into Groups
#'
#' Splits a index vector into separate groups.
#'
#' @param index A vector of integers starting at 1
#'   and increasing in regular increments of 1.
#' @param multiple The number of groups to create.
#' @param values The values to assign per each group.
#' @param label An optional character string - if
#'   specified, returns a list that can be passed
#'   to the operator \code{%j%} with the string
#'   used as the column name.
#'
#' @returns A vector of values. If \code{label} is
#' specified, returns a list with (1) the vector
#' of values, and (2) a character string - the
#' list can then be passed to the operator \code{%j%}.
#'
#' @examples
#'
#' # Create two groups over 4 observations
#' groups( 1:4, 2, 0:1 )
#'
#' # Example input to the %j% operator
#' x <- data.frame( Index.Persons = 1:6 )
#' x %j% groups( x[[1]], 3, LETTERS[1:3], 'X.Condition' )
#'
#' @export

groups <- function( index, multiple, values, label = NULL ) {

  n <- length( index )

  out <- rep( NA, n )

  ord <- c( 1:(multiple-1), 0 )

  for ( j in 1:multiple ) {

    matches <- index %% multiple
    matches <- matches - ord[j]

    out[ matches == 0 ] <- values[j]

  }

  if ( !is.null( label) ) {
    out <- list(
      out,
      label
    )
  }

  return( out )
}


#### 5) intercept ####
#' Create Intercept Column for Design Matrix
#'
#' Convenience function to create inputs
#' for the operator \code{%j%} and add
#' a column with all values equal to 1.
#' Useful to create the intercept column
#' for a design matrix.
#'
#' @param N Number of rows.
#' @param label A character string, the column
#'   name.
#'
#' @return A list that can be passed to the
#'   operator \code{%j%} to add a column
#'   for the intercept in a design matrix.
#'
#' @examples
#'
#' x <- data.frame( Index.Persons = 1:3 )
#' x %j% intercept( 3 )
#'
#' @export

intercept <- function( N, label = 'X.Intercept' ) {

  out <- list(
    rep( 1, N ), label
  )

}

#### 6) get_X ####
#' Extract Design Matrix
#'
#' Extracts the columns for the design matrix
#' from a data frame. Assumes standardized
#' column names (all columns for the design
#' matrix start with a common tag, typically
#' \code{'X.'}).
#'
#' @param dtf A data frame.
#' @param tag A character string, the
#'   sub-string that every column to include
#'   in the design matrix starts with.
#'
#' @returns A matrix.
#'
#' @examples
#'
#' dtf <- data.frame(
#'   Index.Persons = 1:4,
#'   X.Intercept = 1,
#'   X.Condition = rep( 0:1, 2 )
#' )
#' get_X( dtf )
#'
#' @export

get_X <- function( dtf, tag = 'X.' ) {

  i <- grepl( tag, colnames(dtf), fixed = TRUE )

  return( dtf[,i] )
}

