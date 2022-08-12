# Convenience/helper functions
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you
# have any questions or comments
# Last updated 2022-08-10

# Table of contents


#### 1) ... ####

usecalc <- function( what = '',
                     inputs = NULL,
                     distribution = '' ) {

  all_calculations <- rbind(
    c( 'mean', 'beta-binomial' ),
    c( 'mu & phi', 'beta-binomial' )
  )

  if ( what == '' ) {

    message( paste0(
      'Options for arguments:\n',
      paste(
        paste0( '  what = "', all_calculations[,1], '"\n' ),
        paste0( '  distribution = "', all_calculations[,2], '"\n\n' ),
        collapse = ''
      )

    ) )

  }

  out <- NULL

  if ( what == all_calculations[1,1] &
       distribution == all_calculations[1,2] ) {

    out <- inputs[1] / sum( inputs )

  }

  if ( what == all_calculations[1,1] &
       distribution == all_calculations[1,2] ) {

    # alpha <- mu * phi
    # beta <- (1 - mu) * phi

    mu <- inputs[1] / sum( inputs )
    phi <- inputs[1] / mu

    out <- c( mu = mu, phi = phi )

  }

}

#### 2) ... ####

covmat <- function( inputs, action = 'initialize' ) {

  actions <- list(
    Identity = c(
      'Identity matrix', 'identity matrix',
      'Identity', 'identity', 'I'
    ),
    Initialize = c(
      'Initialize', 'initialize',
      'Initial', 'initial', 'init'
    ),
    Cor2Cov = c(
      'Correlation to covariance',
      'correlation to covariance',
      'cor to cov',
      'Cor2cov', 'cor2cov'
    )
  )

  if ( is.null( inputs) ) {
    message( paste0(
      "Possible actions are:\n",
      paste( paste0( '   ', names( actions ) ), collapse = '\n' )
    ) )

    action <- ''
  }

  if ( action == '' ) { return( invisible(NULL) ) }

  if ( action %in% actions$Identity ) {

    return( diag( inputs ) )

  }

  if ( action %in% actions$Initialize ) {

    total <- length( input )

    K <- 0.5 * ( sqrt( 8 * total + 1 ) + 1 )

    Omega <- diag( K )

    inc <- 1
    for ( k in 2:K ) {

      index <- inc:( (inc - 1) + (k-1) )

      Omega[ k, 1:(k-1) ] <- input[ index ]
      Omega[ 1:(k-1), k ] <- input[ index ]

      inc <- max( index ) + 1

    }

    return( Omega )
  }

  if ( action %in% actions$Cor2Cov ) {

    Sigma <- diag( input[[2]] ) %*% input[[1]] %*% diag( input[[2]] )

    return( Sigma )
  }

  stop( 'Check argument "action"; "covmat( NULL )" gives possible options' )
}

