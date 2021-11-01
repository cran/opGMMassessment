# Analysis of a Gaussian mixture structure in one-dimensional numerical data.
# Statistical justification with likelihood ratio and other tests.
#' @importFrom ClusterR GMM
#' @importFrom mclust densityMclust
#' @importFrom mixtools normalmixEM
#' @importFrom parallel detectCores
#' @importFrom NbClust NbClust
#' @importFrom AdaptGauss InformationCriteria4GMM LikelihoodRatio4Mixtures
#' @importFrom methods hasArg
#' @importFrom utils head tail sessionInfo
#' @importFrom stats ks.test rnorm sd
#' @importFrom cluster pam maxSE
#' @importFrom multimode modetest
#' @importFrom DistributionOptimization DistributionOptimization
#' @importFrom mixAK NMixMCMC
#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel stopImplicitCluster
opGMMassessment <- function(Data, FitAlg = "MCMC", Criterion = "LR", MaxModes = 8,
  MaxCores = 2048, PlotIt = FALSE, KS = TRUE, Seed) {

  # Check input
  DIM <- function(...) {
    args <- list(...)
    unlist(lapply(args, function(x) {
      if (is.null(dim(x))) {
        return(1)
      }
      dim(x)[2]
    }))
  }

  if (!hasArg("Data")) {
    stop("opGMMassessment: No data.")
  }
  if (DIM(Data) != 1 | is.numeric(Data) == FALSE) {
    stop("opGMMassessment: Data must be a one-dimensional numerical vector.")
  }
  if (length(Data) < 2) {
    stop("opGMMassessment: Too few data.")
  }
  list.of.FitAlgs <- c("ClusterRGMM", "densityMclust", "DO", "MCMC", "normalmixEM")
  if (!FitAlg %in% list.of.FitAlgs) {
    stop("opGMMassessment: Fit algorithm not implemented. Use ClusterRGMM, DO, densityMclust or normalmixEM")
  }
  list.of.Criteria <- c("AIC", "BIC", "FM", "GAP", "LR", "NbClust", "SI")
  if (!Criterion %in% list.of.Criteria) {
    stop("opGMMassessment: Criterion not implemented. Use AIC, BIC, FM, GAP, LR, NbClust, or SI.")
  }
  if (hasArg("MaxModes")) {
    if (MaxModes < 1) {
      MaxModes <- 8
      warning("opGMMassessment: MaxModes was < 1 and has been set to the default of 8.", call. = FALSE)
    }
  }

  # Create main variables
  DataOrigLength <- length(as.vector(Data))
  DataOrignoNA <- which(!is.na(Data) & !is.infinite(Data))
  Data <- as.vector(Data[DataOrignoNA])
  GMMdata <- Data
  list.of.Modes <- 1:MaxModes
  BestGMM <- 1
  Means <- as.vector(mean(GMMdata, na.rm = TRUE))
  SDs <- as.vector(sd(GMMdata, na.rm = TRUE))
  Weights <- 1
  Mixtures1 <- cbind(Means, SDs, Weights)
  MaxRetries <- 2

  # Internal control functions
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    num_workers <- 2L
  } else {
    num_workers <- parallel::detectCores()
  }
  nProc <- min(num_workers - 1, MaxCores)

  if (!missing(Seed)) {
    ActualSeed <- Seed
  } else {
    ActualSeed <- tail(get(".Random.seed", envir = globalenv()), 1)
  }

  is.integer0 <- function(x) {
    is.integer(x) && length(x) == 0L
  }

  # Functions for mode number determination methods
  idBestGMM_AICBIC <- function(GMMdata, GMMfit, Criterion, ActualSeed) {
    AICBIC <- lapply(1:length(GMMfit), function(x) {
      AICBICi <- AdaptGauss::InformationCriteria4GMM(Data = GMMdata, Means = as.matrix(GMMfit[[x]]$Mixtures)[,
        1], SDs = as.matrix(GMMfit[[x]]$Mixtures)[, 2], Weights = as.matrix(GMMfit[[x]]$Mixtures)[,
        3])
      return(AICBICi)
    })
    AICBIC <- unlist(lapply(1:length(GMMfit), function(x) {
      switch(Criterion, AIC = {
        AICBIC <- AICBIC[[x]]$AIC
      }, BIC = {
        AICBIC <- AICBIC[[x]]$BIC
      })
      return(AICBIC)
    }))
    BestGMM <- which.min(AICBIC)
    return(BestGMM)
  }

  idBestGMM_FM <- function(GMMdata, MaxModes, nProc) {
    list.of.Modes <- 1:(MaxModes - 1)
    switch(Sys.info()[["sysname"]], Windows = {
      ExcessMassi <- lapply(list.of.Modes, function(x) {
        pExM <- multimode::modetest(Data, method = "FM", mod0 = x, B = 60)$p.value
        return(pExM)
      })
    }, {
      ExcessMassi <- parallel::mclapply(list.of.Modes, function(x) {
        pExM <- multimode::modetest(Data, method = "FM", mod0 = x, B = 60)$p.value
        return(pExM)
      }, mc.cores = nProc)
    })
    firstBestGMM <- 1
    for (i in 2:length(ExcessMassi)) {
      if (ExcessMassi[i] < 0.05) {
        firstBestGMM <- i
      } else {
        break
      }
    }
    BestGMM <- firstBestGMM
    return(BestGMM)
  }

  idBestGMM_GAP <- function(GMMdata, MaxModes, ActualSeed, nProc) {
    pam1 <- function(x, k) {
      list(cluster = cluster::pam(x, 3, metric = "euclidean", stand = FALSE,
        cluster.only = TRUE))
    }
    gsPam1 <- clusGapP(x = cbind(GMMdata, GMMdata), FUNcluster = pam1, K.max = MaxModes,
      B = 60, spaceH0 = "original", nProc = nProc)
    BestGMM <- with(gsPam1, maxSE(Tab[, "gap"], Tab[, "SE.sim"], method = "globalSEmax"))
    return(BestGMM)
  }

  idBestGMM_LR <- function(GMMdata, GMMfit, ActualSeed) {
    LRi <- c(1, unlist(lapply(2:MaxModes, function(x) {
      AdaptGauss::LikelihoodRatio4Mixtures(Data = GMMdata, NullMixture = as.matrix(GMMfit[[x -
        1]]$Mixtures), OneMixture = as.matrix(GMMfit[[x]]$Mixtures), PlotIt = FALSE)$Pvalue
    })))
    firstBestGMM <- 1
    for (i in 2:length(LRi)) {
      if (LRi[i] < 0.05) {
        firstBestGMM <- i
      } else {
        break
      }
    }
    BestGMM <- firstBestGMM
    return(BestGMM)
  }

  idBestGMM_NbClust <- function(GMMdata, MaxModes) {
    NBres <- suppressWarnings(NbClust::NbClust(GMMdata, method = "kmeans", distance = "euclidean",
      max.nc = MaxModes))
    BestGMM <- length(unique(NBres$Best.partition))
    return(BestGMM)
  }

  idBestGMM_SI <- function(GMMdata, MaxModes, nProc) {
    list.of.Modes <- 1:(MaxModes - 1)
    switch(Sys.info()[["sysname"]], Windows = {
      ExcessMassi <- lapply(list.of.Modes, function(x) {
        pExM <- multimode::modetest(Data, method = "SI", mod0 = x, B = 60)$p.value
        return(pExM)
      })
    }, {
      ExcessMassi <- parallel::mclapply(list.of.Modes, function(x) {
        pExM <- multimode::modetest(Data, method = "SI", mod0 = x, B = 60)$p.value
        return(pExM)
      }, mc.cores = nProc)
    })
    firstBestGMM <- 1
    for (i in 2:length(ExcessMassi)) {
      if (ExcessMassi[i] < 0.05) {
        firstBestGMM <- i
      } else {
        break
      }
    }
    BestGMM <- firstBestGMM
    return(BestGMM)
  }

  # Start of GMM fit code For reasons of computing speed, three separate
  # versions are available, i.e., for Windows, Linux, and for single-core
  # processing.
  if (var(GMMdata) > 0 & MaxModes > 1) {
    nProc <- min(num_workers - 1, MaxModes, MaxCores)
    if (nProc > 1) {
      switch(Sys.info()[["sysname"]], Windows = {
        requireNamespace("foreach")
        doParallel::registerDoParallel(nProc)
        x <- integer()
        switch(FitAlg, ClusterRGMM = {
          GMMfit <- lapply(list.of.Modes, function(x) {
            set.seed(ActualSeed)
            GMMfit_Mode <- try(ClusterR::GMM(data = data.frame(GMMdata),
            gaussian_comps = list.of.Modes[x], dist_mode = "eucl_dist"),
            TRUE)
            if (class(GMMfit_Mode) != "try-error") {
              Mixtures <- cbind(GMMfit_Mode$centroids, sqrt(GMMfit_Mode$covariance_matrices),
            GMMfit_Mode$weights)
              if (x == 1) {
                Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]),
              ]))
                Mixtures[3] <- 1
              }
            } else {
              GMMfit_Mode <- vector(mode = "list", length = 9)
              Mixtures <- Mixtures1
            }
            return(list(GMMfit_Mode = GMMfit_Mode, Mixtures = Mixtures))
          })
        }, densityMclust = {
          GMMfit <- foreach::foreach(x = list.of.Modes) %dopar% {
            set.seed(ActualSeed)
            GMMfit_Mode <- try(mclust::densityMclust(data = GMMdata, G = x),
            TRUE)
            if (class(GMMfit_Mode) != "try-error") {
              res <- GMMfit_Mode$parameters
              Mixtures <- cbind(unname(res$mean), sqrt(unname(res$variance$sigmasq)),
            unname(res$pro))
              if (x == 1) {
                Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]),
              ]))
                Mixtures[3] <- 1
              }
            } else {
              GMMfit_Mode <- vector(mode = "list", length = 9)
              Mixtures <- Mixtures1
            }
            return(list(GMMfit_Mode = GMMfit_Mode, Mixtures = Mixtures))
          }
        }, DO = {
          GMMfit <- foreach::foreach(x = list.of.Modes) %dopar% {
            set.seed(ActualSeed)
            GMMfit_Mode <- try(DistributionOptimization::DistributionOptimization(Data = GMMdata,
            Modes = list.of.Modes[x], Monitor = 0, CrossoverRate = 0.9,
            ErrorMethod = "chisquare", Seed = ActualSeed), TRUE)
            if (class(GMMfit_Mode) != "try-error") {
              Mixtures <- cbind(GMMfit_Mode$Means, GMMfit_Mode$SDs, GMMfit_Mode$Weights)
              if (x == 1) {
                Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]),
              ]))
                Mixtures[3] <- 1
              }
            } else {
              GMMfit_Mode <- vector(mode = "list", length = 9)
              Mixtures <- Mixtures1
            }
            return(list(GMMfit_Mode = GMMfit_Mode, Mixtures = Mixtures))
          }
        }, MCMC = {
          GMMfit <- foreach::foreach(x = list.of.Modes) %dopar% {
            set.seed(ActualSeed)
            NewSeed <- ActualSeed + 1e+05
            set.seed(ActualSeed)
            FitError <- TRUE
            FitRetriesCount <- 0
            while (FitError == TRUE & FitRetriesCount < MaxRetries) {
              Prior <- list(priorK = "fixed", Kmax = x)
              nMCMC <- c(burn = 5000, keep = 10000, thin = 5, info = 1000)
              GMMfit_Mode <- try(mixAK::NMixMCMC(y0 = GMMdata, prior = Prior,
            nMCMC = nMCMC), TRUE)
              if (class(GMMfit_Mode) != "try-error") FitError <- FALSE
              FitRetriesCount <- FitRetriesCount + 1
              if (class(GMMfit_Mode) == "try-error") set.seed(NewSeed + FitRetriesCount)
            }
            if (class(GMMfit_Mode) != "try-error") {
              MeansMCMC <- GMMfit_Mode[[1]]$poster.mean.mu * GMMfit_Mode[[1]]$scale$scale +
            GMMfit_Mode[[1]]$scale$shift
              SDsMCMC <- sqrt(GMMfit_Mode[[1]]$scale$scale ^ 2 * as.numeric(GMMfit_Mode[[1]]$poster.mean.Sigma))
              WeightsMCMC <- GMMfit_Mode[[1]]$poster.mean.w[seq(0, (GMMfit_Mode[[1]]$nx_w -
            1) * GMMfit_Mode[[1]]$prior$Kmax, by = GMMfit_Mode[[1]]$prior$Kmax) +
            c(1:x)]
              Mixtures <- cbind(MeansMCMC, SDsMCMC, WeightsMCMC)
              if (x == 1) {
                Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]),
              ]))
                Mixtures[3] <- 1
              }
            } else {
              GMMfit_Mode <- vector(mode = "list", length = 9)
              Mixtures <- Mixtures1
            }
            set.seed(ActualSeed)
            return(list(GMMfit_Mode = GMMfit_Mode, Mixtures = Mixtures))
          }
        }, normalmixEM = {
          GMMfit <- foreach::foreach(x = list.of.Modes) %dopar% {
            set.seed(ActualSeed)
            NewSeed <- ActualSeed + 1e+05
            set.seed(ActualSeed)
            FitError <- TRUE
            FitRetriesCount <- 0
            while (FitError == TRUE & FitRetriesCount < MaxRetries) {
              GMMfit_Mode <- try(mixtools::normalmixEM(GMMdata, mu = kmeans(GMMdata,
            x)$centers, ECM = TRUE, maxrestarts = 1e+05), TRUE)
              if (class(GMMfit_Mode) != "try-error") FitError <- FALSE
              FitRetriesCount <- FitRetriesCount + 1
              if (class(GMMfit_Mode) == "try-error") set.seed(NewSeed + FitRetriesCount)
            }
            if (class(GMMfit_Mode) != "try-error") {
              Mixtures <- cbind(GMMfit_Mode$mu, GMMfit_Mode$sigma, GMMfit_Mode$lambda)
              if (x == 1) {
                Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]),
              ]))
                Mixtures[3] <- 1
              }
            } else {
              GMMfit_Mode <- vector(mode = "list", length = 9)
              Mixtures <- Mixtures1
            }
            set.seed(ActualSeed)
            return(list(GMMfit_Mode = GMMfit_Mode, Mixtures = Mixtures))
          }
        })
        doParallel::stopImplicitCluster()
      }, {
        switch(FitAlg, ClusterRGMM = {
          GMMfit <- lapply(list.of.Modes, function(x) {
            set.seed(ActualSeed)
            GMMfit_Mode <- try(ClusterR::GMM(data = data.frame(GMMdata),
            gaussian_comps = list.of.Modes[x], dist_mode = "eucl_dist"),
            TRUE)
            if (class(GMMfit_Mode) != "try-error") {
              Mixtures <- cbind(GMMfit_Mode$centroids, sqrt(GMMfit_Mode$covariance_matrices),
            GMMfit_Mode$weights)
              if (x == 1) {
                Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]),
              ]))
                Mixtures[3] <- 1
              }
            } else {
              GMMfit_Mode <- vector(mode = "list", length = 9)
              Mixtures <- Mixtures1
            }
            return(list(GMMfit_Mode = GMMfit_Mode, Mixtures = Mixtures))
          })
        }, densityMclust = {
          GMMfit <- parallel::mclapply(list.of.Modes, function(x) {
            set.seed(ActualSeed)
            GMMfit_Mode <- try(mclust::densityMclust(data = GMMdata, G = x),
            TRUE)
            if (class(GMMfit_Mode) != "try-error") {
              res <- GMMfit_Mode$parameters
              Mixtures <- cbind(unname(res$mean), sqrt(unname(res$variance$sigmasq)),
            unname(res$pro))
              if (x == 1) {
                Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]),
              ]))
                Mixtures[3] <- 1
              }
            } else {
              GMMfit_Mode <- vector(mode = "list", length = 9)
              Mixtures <- Mixtures1
            }
            return(list(GMMfit_Mode = GMMfit_Mode, Mixtures = Mixtures))
          }, mc.cores = nProc)
        }, DO = {
          GMMfit <- parallel::mclapply(list.of.Modes, function(x) {
            set.seed(ActualSeed)
            GMMfit_Mode <- try(DistributionOptimization::DistributionOptimization(Data = GMMdata,
            Modes = list.of.Modes[x], Monitor = 0, CrossoverRate = 0.9,
            ErrorMethod = "chisquare", Seed = ActualSeed), TRUE)
            if (class(GMMfit_Mode) != "try-error") {
              Mixtures <- cbind(GMMfit_Mode$Means, GMMfit_Mode$SDs, GMMfit_Mode$Weights)
              if (x == 1) {
                Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]),
              ]))
                Mixtures[3] <- 1
              }
            } else {
              GMMfit_Mode <- vector(mode = "list", length = 9)
              Mixtures <- Mixtures1
            }
            return(list(GMMfit_Mode = GMMfit_Mode, Mixtures = Mixtures))
          }, mc.cores = nProc)
        }, MCMC = {
          GMMfit <- parallel::mclapply(list.of.Modes, function(x) {
            set.seed(ActualSeed)
            NewSeed <- ActualSeed + 1e+05
            set.seed(ActualSeed)
            FitError <- TRUE
            FitRetriesCount <- 0
            while (FitError == TRUE & FitRetriesCount < MaxRetries) {
              Prior <- list(priorK = "fixed", Kmax = x)
              nMCMC <- c(burn = 5000, keep = 10000, thin = 5, info = 1000)
              GMMfit_Mode <- try(mixAK::NMixMCMC(y0 = GMMdata, prior = Prior,
            nMCMC = nMCMC), TRUE)
              if (class(GMMfit_Mode) != "try-error") FitError <- FALSE
              FitRetriesCount <- FitRetriesCount + 1
              if (class(GMMfit_Mode) == "try-error") set.seed(NewSeed + FitRetriesCount)
            }
            if (class(GMMfit_Mode) != "try-error") {
              MeansMCMC <- GMMfit_Mode[[1]]$poster.mean.mu * GMMfit_Mode[[1]]$scale$scale +
            GMMfit_Mode[[1]]$scale$shift
              SDsMCMC <- sqrt(GMMfit_Mode[[1]]$scale$scale ^ 2 * as.numeric(GMMfit_Mode[[1]]$poster.mean.Sigma))
              WeightsMCMC <- GMMfit_Mode[[1]]$poster.mean.w[seq(0, (GMMfit_Mode[[1]]$nx_w -
            1) * GMMfit_Mode[[1]]$prior$Kmax, by = GMMfit_Mode[[1]]$prior$Kmax) +
            c(1:x)]
              Mixtures <- cbind(MeansMCMC, SDsMCMC, WeightsMCMC)
              if (x == 1) {
                Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]),
              ]))
                Mixtures[3] <- 1
              }
            } else {
              GMMfit_Mode <- vector(mode = "list", length = 9)
              Mixtures <- Mixtures1
            }
            set.seed(ActualSeed)
            return(list(GMMfit_Mode = GMMfit_Mode, Mixtures = Mixtures))
          }, mc.cores = nProc)
        }, normalmixEM = {
          GMMfit <- parallel::mclapply(list.of.Modes, function(x) {
            set.seed(ActualSeed)
            NewSeed <- ActualSeed + 1e+05
            set.seed(ActualSeed)
            FitError <- TRUE
            FitRetriesCount <- 0
            while (FitError == TRUE & FitRetriesCount < MaxRetries) {
              GMMfit_Mode <- try(mixtools::normalmixEM(GMMdata, mu = kmeans(GMMdata,
            x)$centers, ECM = TRUE, maxrestarts = 1e+05), TRUE)
              if (class(GMMfit_Mode) != "try-error") FitError <- FALSE
              FitRetriesCount <- FitRetriesCount + 1
              if (class(GMMfit_Mode) == "try-error") set.seed(NewSeed + FitRetriesCount)
            }
            if (class(GMMfit_Mode) != "try-error") {
              Mixtures <- cbind(GMMfit_Mode$mu, GMMfit_Mode$sigma, GMMfit_Mode$lambda)
              if (x == 1) {
                Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]),
              ]))
                Mixtures[3] <- 1
              }
            } else {
              GMMfit_Mode <- vector(mode = "list", length = 9)
              Mixtures <- Mixtures1
            }
            set.seed(ActualSeed)
            return(list(GMMfit_Mode = GMMfit_Mode, Mixtures = Mixtures))
          }, mc.cores = nProc)
        })
      })
    } else {
      switch(FitAlg, ClusterRGMM = {
        GMMfit <- lapply(list.of.Modes, function(x) {
          set.seed(ActualSeed)
          GMMfit_Mode <- try(ClusterR::GMM(data = data.frame(GMMdata), gaussian_comps = list.of.Modes[x],
          dist_mode = "eucl_dist"), TRUE)
          if (class(GMMfit_Mode) != "try-error") {
            Mixtures <- cbind(GMMfit_Mode$centroids, sqrt(GMMfit_Mode$covariance_matrices),
            GMMfit_Mode$weights)
            if (x == 1) {
              Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]),
            ]))
              Mixtures[3] <- 1
            }
          } else {
            GMMfit_Mode <- vector(mode = "list", length = 9)
            Mixtures <- Mixtures1
          }
          return(list(GMMfit_Mode = GMMfit_Mode, Mixtures = Mixtures))
        })
      }, densityMclust = {
        GMMfit <- lapply(list.of.Modes, function(x) {
          set.seed(ActualSeed)
          GMMfit_Mode <- try(mclust::densityMclust(data = GMMdata, G = x),
          TRUE)
          if (class(GMMfit_Mode) != "try-error") {
            res <- GMMfit_Mode$parameters
            Mixtures <- cbind(unname(res$mean), sqrt(unname(res$variance$sigmasq)),
            unname(res$pro))
            if (x == 1) {
              Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]),
            ]))
              Mixtures[3] <- 1
            }
          } else {
            GMMfit_Mode <- vector(mode = "list", length = 9)
            Mixtures <- Mixtures1
          }
          return(list(GMMfit_Mode = GMMfit_Mode, Mixtures = Mixtures))
        })
      }, DO = {
        GMMfit <- lapply(list.of.Modes, function(x) {
          set.seed(ActualSeed)
          GMMfit_Mode <- try(DistributionOptimization::DistributionOptimization(Data = GMMdata,
          Modes = list.of.Modes[x], Monitor = 0, CrossoverRate = 0.9, ErrorMethod = "chisquare",
          Seed = ActualSeed), TRUE)
          if (class(GMMfit_Mode) != "try-error") {
            Mixtures <- cbind(GMMfit_Mode$Means, GMMfit_Mode$SDs, GMMfit_Mode$Weights)
            if (x == 1) {
              Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]),
            ]))
              Mixtures[3] <- 1
            }
          } else {
            GMMfit_Mode <- vector(mode = "list", length = 9)
            Mixtures <- Mixtures1
          }
          return(list(GMMfit_Mode = GMMfit_Mode, Mixtures = Mixtures))
        })
      }, MCMC = {
        GMMfit <- lapply(list.of.Modes, function(x) {
          set.seed(ActualSeed)
          NewSeed <- ActualSeed + 1e+05
          set.seed(ActualSeed)
          FitError <- TRUE
          FitRetriesCount <- 0
          while (FitError == TRUE & FitRetriesCount < MaxRetries) {
            Prior <- list(priorK = "fixed", Kmax = x)
            nMCMC <- c(burn = 5000, keep = 10000, thin = 5, info = 1000)
            GMMfit_Mode <- try(mixAK::NMixMCMC(y0 = GMMdata, prior = Prior,
            nMCMC = nMCMC), TRUE)
            if (class(GMMfit_Mode) != "try-error") FitError <- FALSE
            FitRetriesCount <- FitRetriesCount + 1
            if (class(GMMfit_Mode) == "try-error") set.seed(NewSeed + FitRetriesCount)
          }
          if (class(GMMfit_Mode) != "try-error") {
            MeansMCMC <- GMMfit_Mode[[1]]$poster.mean.mu * GMMfit_Mode[[1]]$scale$scale +
            GMMfit_Mode[[1]]$scale$shift
            SDsMCMC <- sqrt(GMMfit_Mode[[1]]$scale$scale ^ 2 * as.numeric(GMMfit_Mode[[1]]$poster.mean.Sigma))
            WeightsMCMC <- GMMfit_Mode[[1]]$poster.mean.w[seq(0, (GMMfit_Mode[[1]]$nx_w -
            1) * GMMfit_Mode[[1]]$prior$Kmax, by = GMMfit_Mode[[1]]$prior$Kmax) +
            c(1:x)]
            Mixtures <- cbind(MeansMCMC, SDsMCMC, WeightsMCMC)
            if (x == 1) {
              Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]),
            ]))
              Mixtures[3] <- 1
            }
          } else {
            GMMfit_Mode <- vector(mode = "list", length = 9)
            Mixtures <- Mixtures1
          }
          set.seed(ActualSeed)
          return(list(GMMfit_Mode = GMMfit_Mode, Mixtures = Mixtures))
        })
      }, normalmixEM = {
        GMMfit <- lapply(list.of.Modes, function(x) {
          set.seed(ActualSeed)
          NewSeed <- ActualSeed + 1e+05
          set.seed(ActualSeed)
          FitError <- TRUE
          FitRetriesCount <- 0
          while (FitError == TRUE & FitRetriesCount < MaxRetries) {
            GMMfit_Mode <- try(mixtools::normalmixEM(GMMdata, mu = kmeans(GMMdata,
            x)$centers, ECM = TRUE, maxrestarts = 1e+05), TRUE)
            if (class(GMMfit_Mode) != "try-error") FitError <- FALSE
            FitRetriesCount <- FitRetriesCount + 1
            if (class(GMMfit_Mode) == "try-error") set.seed(NewSeed + FitRetriesCount)
          }
          if (class(GMMfit_Mode) != "try-error") {
            Mixtures <- cbind(GMMfit_Mode$mu, GMMfit_Mode$sigma, GMMfit_Mode$lambda)
            if (x == 1) {
              Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]),
            ]))
              Mixtures[3] <- 1
            }
          } else {
            GMMfit_Mode <- vector(mode = "list", length = 9)
            Mixtures <- Mixtures1
          }
          set.seed(ActualSeed)
          return(list(GMMfit_Mode = GMMfit_Mode, Mixtures = Mixtures))
        })
      })
    }

    # Identify best fit based on selected criterion
    switch(Criterion, AIC = {
      BestGMM <- idBestGMM_AICBIC(GMMdata = GMMdata, GMMfit = GMMfit, Criterion = Criterion,
        ActualSeed = ActualSeed)
    }, BIC = {
      BestGMM <- idBestGMM_AICBIC(GMMdata = GMMdata, GMMfit = GMMfit, Criterion = Criterion,
        ActualSeed = ActualSeed)
    }, FM = {
      BestGMM <- idBestGMM_FM(GMMdata = GMMdata, MaxModes = MaxModes, nProc = nProc)
    }, GAP = {
      BestGMM <- idBestGMM_GAP(GMMdata = GMMdata, MaxModes = MaxModes, nProc = nProc,
        ActualSeed = ActualSeed)
    }, LR = {
      BestGMM <- idBestGMM_LR(GMMdata = GMMdata, GMMfit = GMMfit, ActualSeed = ActualSeed)
    }, NbClust = {
      BestGMM <- idBestGMM_NbClust(GMMdata = GMMdata, MaxModes = MaxModes)
    }, SI = {
      BestGMM <- idBestGMM_SI(GMMdata = GMMdata, MaxModes = MaxModes, nProc = nProc)
    }, BestGMM <- BestGMM)

    # Extract GMM parameters
    Means <- as.vector(GMMfit[[BestGMM]]$Mixtures[, 1])
    SDs <- as.vector(GMMfit[[BestGMM]]$Mixtures[, 2])
    Weights <- as.vector(GMMfit[[BestGMM]]$Mixtures[, 3])
  }
  # Calculate Bayes boundaries
  Boundaries <- c()
  ClassesB <- rep(1, length(GMMdata))
  if (BestGMM > 1) {
    Boundaries <- AdaptGauss::BayesDecisionBoundaries(Means = Means, SDs = SDs,
      Weights = Weights)
    if (is.integer0(Boundaries) == FALSE) {
      Boundaries <- Boundaries[Boundaries >= min(Means) & Boundaries <= max(Means)]
    }
    if (length(Boundaries) > 0) {
      ClassesB <- cutGMM(x = GMMdata, breaks = Boundaries)
    }
  }

  Classes <- rep(NA, length(DataOrigLength))
  Classes[DataOrignoNA] <- ClassesB

  # Do Kolmogorov-Smirnov test
  if (KS == TRUE) {
    set.seed(ActualSeed)
    Pred <- CreateGMM(Means = Means, SDs = SDs, Weights = Weights, n = 1000)$Data
    KStest <- suppressWarnings(ks.test(x = GMMdata, y = Pred))
  } else {
    KStest <- NA
  }

  # Prepare plot
  p1 <- GMMplotGG(Data = GMMdata, Means = Means, SDs = SDs, Weights = Weights,
    Hist = TRUE)
  if (PlotIt == TRUE) {
    print(p1)
  }

  # Return results
  return(list(Cls = Classes, Means = Means, SDs = SDs, Weights = Weights, Boundaries = Boundaries,
    Plot = p1, KS = KStest))
}
