# Cuts a vector into separate groups according to input breaks
cutGMM <- function( x, breaks, right = TRUE ) {

  DIM <- function( ... ) {
    args <- list( ... )
    unlist( lapply( args, function( x ) {
      if ( is.null( dim( x ) ) ) {
        return( 1 )
      }
      dim( x )[2]
    } ) )
  }

  if ( DIM( x ) != 1 |
    is.numeric( x ) == FALSE |
    length( x ) == 0 ) {
    stop( "CutGMM: x should be a single numerical vector of length > 0." )
  }

  if ( length( breaks ) > 0 ) {
    if ( hasArg( right ) == FALSE ) {
      right <- TRUE
    }
    breaks <- sort( breaks )
    if ( right == TRUE ) {
      if ( length( breaks ) > 1 ) {
        compMat <- vapply( x, function( x ) x > breaks, logical( length( breaks ) ) )
        ClassesGMM <- colSums( compMat ) + 1
      } else {
        ClassesGMM <- as.numeric( x > breaks ) + 1
      }
    } else {
      if ( length( breaks ) > 1 ) {
        compMat <- vapply( x, function( x ) x >= breaks, logical( length( breaks ) ) )
        ClassesGMM <- colSums( compMat ) + 1
      } else {
        ClassesGMM <- as.numeric( x >= breaks ) + 1
      }
    }
  } else {
    warning( "cutGMM: No breaks provided. Assuming one single class.", call. = FALSE )
  }

  return( ClassesGMM )
}
