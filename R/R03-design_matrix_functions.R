# Design matrix functions
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you
# have any questions or comments
# Last updated 2022-08-21

# Table of contents
# 1) Operators
#   1.1) `%j%`
# 2) Functions to duplicate values
#   2.1) each
#   1.2) zeros
#   2.3) ones
# 3) Functions to generate data frames and columns
#   3.1) pse
#   3.2) groups
#   3.3) intercept
#   3.4) coding
# 4) Functions to extract and manipulate columns
#   4.1) get_X
#   4.2) z_score

#### 1) Operators ####

#### 1.1) `%j%` ####
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


#### 2) Functions to duplicate values ####

#### 2.1) each ####
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

#### 1.2) zeros ####
#' Create Vector of Zeros
#'
#' Convenience function that generates
#' a vector or matrix of zeros.
#'
#' @param N Number of elements or rows.
#' @param K Number of columns if output
#'   should be a matrix.
#'
#' @returns Either a vector or matrix of zeros.
#'
#' @examples
#'
#' zeros( 3 )
#' zeros( 3, 2 )
#'
#' @export

zeros <- function( N, K = NULL ) {

  if ( is.null( K ) ) {

    return( rep( 0, N ) )

  } else {

    return( matrix( 0, N, K ) )

  }

}

#### 2.3) ones ####
#' Create Vector of Ones
#'
#' Convenience function that generates
#' a vector or matrix of ones
#'
#' @param N Number of elements or rows.
#' @param K Number of columns if output
#'   should be a matrix.
#'
#' @returns Either a vector or matrix of ones
#'
#' @examples
#'
#' ones( 3 )
#' ones( 3, 2 )
#'
#' @export

ones <- function( N ) {

  if ( is.null( K ) ) {

    return( rep( 1, N ) )

  } else {

    return( matrix( 1, N, K ) )

  }

}

#### 3) Functions to generate data frames and columns ####

#### 3.1) pse ####
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

#### 3.2) groups ####
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

#### 3.3) intercept ####
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

#### 3.4) coding ####
#' Recode a Variable Based on Matches
#'
#' Function that recodes a variable (e.g.,
#' assigning 0 and 1 to create a dummy-coded
#' variable based on matches).
#'
#' @param x A vector of values.
#' @param to_match Either (a) a vector of values or
#'   (b) a list of vectors, the elements in \code{x}
#'   to match when assigning codes.
#' @param codes The new values or codes to assign based
#'   on the matches between \code{x} and \code{to_match}.
#'   If \code{codes == NULL} then if assigns codes of (a)
#'   0 and 1 for two levels of \code{x}, or (b) codes of
#'   -1, 0, and 1 for three levels of \code{x}.
#' @param label An optional character string, the label
#'   for the column to add to a data frame (allows the
#'   output from the function to be passed to the
#'   \code{%j%} operator).
#'
#' @returns A vector of new values. If \code{label} is
#' specified, returns a list with the vector of new
#' values and the label for the column to add to a data
#' frame.
#'
#' @examples
#'
#' # Define a vector of values
#' x <- rep( 1:4, each = 2 )
#'
#' # Create a dummy-coded variable (0 or 1)
#' coding( x, 3:4 )
#' # Create an effects-coded variable (-1, 0, 1)
#' coding( x, list( 1, 2:3, 4 ) )
#' # Custom codes (alternative variant of effects coding)
#' coding( x, list( 1:2, 3:4 ), c( -.5, .5 ) )
#'
#' # How to add a new column to a data frame using %j%
#' dtf <- data.frame( Index = x )
#' dtf <- dtf %j%
#'   coding( dtf$Index, 3:4, label = 'X.1' )
#' print( dtf )
#'
#' @export

coding <- function( x, to_match, codes = NULL, label = NULL ) {

  # Number of options
  N <- length( x )

  # Initialize output
  out <- rep( NA, N )

  # Determine structure of values to match over
  is_list <- is.list( to_match )

  # If a list of values is provided
  if ( is_list ) {

    # Number of sets to match over
    L <- length( to_match )

    # If no user-defined codes are provided
    if ( is.null( codes ) ) {

      if ( L == 2 ) {
        codes <- 0:1
      }
      if ( L == 3 ) {
        codes <- c( -1, 0, 1 )
      }

      # Close 'If no user-defined codes are provided'
    }

    N_codes <- length( codes )

    # If all codes and matching sets specified
    if ( N_codes == L ) {

      # Loop over matching sets
      for ( l in 1:L ) {

        matches <- x %in% to_match[[l]]

        out[ matches ] <- codes[l]

        # Close 'Loop over matching sets'
      }

      # Close 'If all codes and matching sets specified'
    }

    # If all non-zero codes and matching sets provided
    if ( (N_codes - 1) == L ) {

      # Vector of zeros
      out <- zeros( N )

      # Loop over matching sets
      for ( l in 1:L ) {

        matches <- x %in% to_match[[l]]

        out[ matches ] <- codes[l]

        # Close 'Loop over matching sets'
      }

      # Close 'If all non-zero codes and matching sets provided'
    }

    # Close 'If a list of values is provided'
  } else {

    matches <- x %in% to_match

    if ( is.null( codes ) ) {
      codes <- 0:1
    }

    out[ !matches ] <- codes[1]
    out[ matches ] <- codes[2]

    # Close else for 'If a list of values is provided'
  }

  if ( all( is.na( out ) ) ) {
    stop( 'FORTHCOMING')
  }

  if ( !is.null( label ) ) {

    return( list( out, label ) )

  }

  return( out )
}

#### 4) Functions to extract and manipulate columns ####

#### 4.1) get_X ####
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

#### 4.2) z_score ####
#' Compute Z-scores
#'
#' Convenience function that converts a
#' vector of values into a set of
#' z-scores (mean of 0 and SD of 1).
#'
#' @param values A vector of numeric values.
#' @param m An optional numeric value, the mean
#'   to use when computing the z-scores.
#' @param s An optional numeric value, the
#'   standard deviation to use when computing
#'   the z-scores.
#'
#' @returns A vector of z-scores.
#'
#' @examples
#'
#' # For reproducibility
#' set.seed( 111 )
#' # Simulate values with M = 100, SD = 15
#' x <- rnorm( 5, 100, 15 )
#' # Convert to z-scores
#' z <- z_score( x, m = 100, s = 15 )
#'
#' @export

z_score <- function( values, m = NULL, s = NULL ) {

  if ( is.null(m) ) {
    m <- mean(values, na.rm = TRUE)
  }
  if ( is.null(s) ) {
    s <- sd(values, na.rm = TRUE)
  }

  z <- (values - m)/s

  return(z)
}

