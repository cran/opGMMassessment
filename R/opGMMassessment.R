# Analysis of a Gaussian mixture structure in one-dimensional numerical data.
# Statistical justification with likelihood ratio and other tests.
#' @useDynLib(opGMMassessment, .registration = TRUE)
#' @importFrom ClusterR GMM
#' @importFrom mclust densityMclust
#' @importFrom mixtools normalmixEM
#' @importFrom parallel detectCores
#' @importFrom methods hasArg
#' @importFrom utils head tail sessionInfo
#' @importFrom stats ks.test rnorm sd
#' @importFrom DistributionOptimization DistributionOptimization
#' @importFrom mixAK NMixMCMC
#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @export
opGMMassessment <- function( Data, FitAlg = "MCMC", Criterion = "LR", MaxModes = 8,
                             MaxCores = getOption( "mc.cores", 2L ), PlotIt = FALSE, KS = TRUE, Seed ) {

  # Check input
  DIM <- function( ... ) {
    args <- list( ... )
    unlist( lapply( args, function( x ) {
      if ( is.null( dim( x ) ) ) {
        return( 1 )
      }
      dim( x )[2]
    } ) )
  }

  if ( !hasArg( "Data" ) ) {
    stop( "opGMMassessment: No data." )
  }
  if ( DIM( Data ) != 1 | is.numeric( Data ) == FALSE ) {
    stop( "opGMMassessment: Data must be a one-dimensional numerical vector." )
  }
  if ( length( Data ) < 2 ) {
    stop( "opGMMassessment: Too few data." )
  }
  list.of.FitAlgs <- c( "ClusterRGMM", "densityMclust", "DO", "MCMC", "normalmixEM" )
  if ( !FitAlg %in% list.of.FitAlgs ) {
    stop( "opGMMassessment: Fit algorithm not implemented. Use ClusterRGMM, DO, densityMclust or normalmixEM" )
  }
  list.of.Criteria <- c( "AIC", "BIC", "FM", "GAP", "LR", "NbClust", "SI" )
  if ( !Criterion %in% list.of.Criteria ) {
    stop( "opGMMassessment: Criterion not implemented. Use AIC, BIC, FM, GAP, LR, NbClust, or SI." )
  }
  list.of.Criteria.directModenumber <- c( "FM", "GAP", "NbClust", "SI" )
  if ( hasArg( "MaxModes" ) ) {
    if ( MaxModes < 1 ) {
      MaxModes <- 8
      warning( "opGMMassessment: MaxModes was < 1 and has been set to the default of 8.",
               call. = FALSE )
    }
  }

  # Create main variables
  DataOrigLength <- length( as.vector( Data ) )
  DataOrignoNA <- which( !is.na( Data ) & !is.infinite( Data ) )
  Data <- as.vector( Data[DataOrignoNA] )
  GMMdata <- Data
  list.of.Modes <- 1:MaxModes
  BestGMM <- 1
  Means <- as.vector( mean( GMMdata, na.rm = TRUE ) )
  SDs <- as.vector( sd( GMMdata, na.rm = TRUE ) )
  Weights <- 1
  Mixtures1 <- cbind( Means, SDs, Weights )
  MaxRetries <- 3

  # Internal control functions
  num_workers <- parallel::detectCores( )
  nProc <- min( num_workers - 1, MaxCores )

  if ( !missing( Seed ) ) {
    ActualSeed <- Seed
  } else {
    ActualSeed <- tail( get( ".Random.seed", envir = globalenv( ) ), 1 )
  }

  # Start of GMM fit code. For reasons of computing speed, three separate
  # versions are available, i.e., for Windows, Linux, and for single-core
  # processing.
  if ( var( GMMdata ) > 0 ) {
    if ( !Criterion %in% list.of.Criteria.directModenumber ) {
      # GMM fits were the number of modes is determined by comparing goodnbess
      # of fit
      if ( nProc > 1 & MaxModes > 1 ) {
        switch( Sys.info( )[["sysname"]], Windows = {
          doParallel::registerDoParallel( nProc )
          x <- integer( )
          GMMfit <- lapply( list.of.Modes, function( x ) {
            PerformGMMfit( GMMdata = GMMdata, FitAlg = FitAlg, Modes = list.of.Modes[x],
                           Mixtures1 = Mixtures1, ActualSeed = ActualSeed, MaxRetries = MaxRetries )
          } )
        }, {
          GMMfit <- parallel::mclapply( list.of.Modes, function( x ) {
            PerformGMMfit( GMMdata = GMMdata, FitAlg = FitAlg, Modes = list.of.Modes[x],
                           Mixtures1 = Mixtures1, ActualSeed = ActualSeed, MaxRetries = MaxRetries )
          } )
        } )
      } else {
        GMMfit <- lapply( list.of.Modes, function( x ) {
          PerformGMMfit( GMMdata = GMMdata, FitAlg = FitAlg, Modes = list.of.Modes[x],
                         Mixtures1 = Mixtures1, ActualSeed = ActualSeed, MaxRetries = MaxRetries )
        } )
      }

      BestGMM <- DetermineBestGMM( GMMdata = GMMdata, GMMfit = GMMfit, Criterion = Criterion,
                                   ActualSeed = ActualSeed, MaxModes = MaxModes, nProc = nProc, BestGMM = BestGMM )

      Means <- as.vector( GMMfit[[BestGMM]]$Mixtures[, 1] )
      SDs <- as.vector( GMMfit[[BestGMM]]$Mixtures[, 2] )
      Weights <- as.vector( GMMfit[[BestGMM]]$Mixtures[, 3] )
    } else {
      # GMM fits were the number of modes is determined prior to fitting
      BestGMM <- DetermineBestGMM( GMMdata = GMMdata, GMMfit = NA, Criterion = Criterion,
                                   ActualSeed = ActualSeed, MaxModes = MaxModes, nProc = nProc, BestGMM = BestGMM )
      GMMfit <- PerformGMMfit( GMMdata = Data, FitAlg = FitAlg, Modes = BestGMM,
                               Mixtures1 = Mixtures1, ActualSeed = ActualSeed, MaxRetries = MaxRetries )

      Means <- as.vector( GMMfit$Mixtures[, 1] )
      SDs <- as.vector( GMMfit$Mixtures[, 2] )
      Weights <- as.vector( GMMfit$Mixtures[, 3] )
    }
  }

  # Calculate Bayes boundaries
  is.integer0 <- function( x ) {
    is.integer( x ) && length( x ) == 0L
  }

  Boundaries <- c( )
  ClassesB <- rep( 1, length( GMMdata ) )
  if ( BestGMM > 1 ) {
    Boundaries <- AdaptGauss::BayesDecisionBoundaries( Means = Means, SDs = SDs,
                                                       Weights = Weights )
    if ( is.integer0( Boundaries ) == FALSE ) {
      Boundaries <- Boundaries #[Boundaries >= min(Means) & Boundaries <= max(Means)]
    }
    if ( length( Boundaries ) > 0 ) {
      ClassesB <- cutGMM( x = GMMdata, breaks = Boundaries )
    }
  }

  Classes <- rep( NA, length( DataOrigLength ) )
  Classes[DataOrignoNA] <- ClassesB

  # Do Kolmogorov-Smirnov test
  if ( KS == TRUE ) {
    set.seed( ActualSeed )
    Pred <- CreateGMM( Means = Means, SDs = SDs, Weights = Weights, n = 1000 )$Data
    KStest <- suppressWarnings( ks.test( x = GMMdata, y = Pred ), classes = "warning" )
  } else {
    KStest <- NA
  }

  # Prepare plot
  p1 <- GMMplotGG( Data = GMMdata, Means = Means, SDs = SDs, Weights = Weights,
                   Hist = TRUE )
  if ( PlotIt == TRUE ) {
    print( p1 )
  }

  # Return results
  return( list( Cls = Classes, Means = Means, SDs = SDs, Weights = Weights, Boundaries = Boundaries,
                Plot = p1, KS = KStest ) )
}
