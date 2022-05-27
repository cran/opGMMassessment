# Perfo the GMM fit using the selected algorithm.
#' @importFrom ClusterR GMM
#' @importFrom mclust densityMclust
#' @importFrom mixtools normalmixEM
#' @importFrom utils head tail sessionInfo
#' @importFrom DistributionOptimization DistributionOptimization
#' @importFrom mixAK NMixMCMC
#' @importFrom methods is
PerformGMMfit <- function(FitAlg, GMMdata, Modes, Mixtures1, ActualSeed, MaxRetries) {
  set.seed(ActualSeed)
  switch(FitAlg, ClusterRGMM = {
    GMMfit_Mode <- try(ClusterR::GMM(data = data.frame(GMMdata), gaussian_comps = Modes,
      dist_mode = "eucl_dist"), TRUE)
    if (!is(GMMfit_Mode, "try-error")) {
      Mixtures <- cbind(GMMfit_Mode$centroids, sqrt(GMMfit_Mode$covariance_matrices),
        GMMfit_Mode$weights)
      if (Modes == 1) {
        Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]), ]))
        Mixtures[3] <- 1
      }
    } else {
      GMMfit_Mode <- vector(mode = "list", length = 9)
      Mixtures <- Mixtures1
    }
  }, densityMclust = {
    GMMfit_Mode <- try(mclust::densityMclust(data = GMMdata, G = Modes), TRUE)
    if (!is(GMMfit_Mode, "try-error")) {
      res <- GMMfit_Mode$parameters
      Mixtures <- cbind(unname(res$mean), sqrt(unname(res$variance$sigmasq)),
        unname(res$pro))
      if (Modes == 1) {
        Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]), ]))
        Mixtures[3] <- 1
      }
    } else {
      GMMfit_Mode <- vector(mode = "list", length = 9)
      Mixtures <- Mixtures1
    }
  }, DO = {
    GMMfit_Mode <- try(DistributionOptimization::DistributionOptimization(Data = GMMdata,
      Modes = Modes, Monitor = 0, CrossoverRate = 0.9, ErrorMethod = "chisquare",
      Seed = ActualSeed), TRUE)
    if (!is(GMMfit_Mode, "try-error")) {
      Mixtures <- cbind(GMMfit_Mode$Means, GMMfit_Mode$SDs, GMMfit_Mode$Weights)
      if (Modes == 1) {
        Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]), ]))
        Mixtures[3] <- 1
      }
    } else {
      GMMfit_Mode <- vector(mode = "list", length = 9)
      Mixtures <- Mixtures1
    }
  }, MCMC = {
    NewSeed <- ActualSeed + 1e+05
    FitError <- TRUE
    FitRetriesCount <- 0
    while (FitError == TRUE & FitRetriesCount < MaxRetries) {
      Prior <- list(priorK = "fixed", Kmax = Modes)
      nMCMC <- c(burn = 5000, keep = 10000, thin = 5, info = 1000)
      GMMfit_Mode <- try(mixAK::NMixMCMC(y0 = GMMdata, prior = Prior, nMCMC = nMCMC),
        TRUE)
      if (!is(GMMfit_Mode, "try-error")) FitError <- FALSE
      FitRetriesCount <- FitRetriesCount + 1
      if (is(GMMfit_Mode, "try-error")) set.seed(NewSeed + FitRetriesCount)
    }
    if (!is(GMMfit_Mode, "try-error")) {
      MeansMCMC <- GMMfit_Mode[[1]]$poster.mean.mu * GMMfit_Mode[[1]]$scale$scale +
        GMMfit_Mode[[1]]$scale$shift
      SDsMCMC <- sqrt(GMMfit_Mode[[1]]$scale$scale^2 * as.numeric(GMMfit_Mode[[1]]$poster.mean.Sigma))
      WeightsMCMC <- GMMfit_Mode[[1]]$poster.mean.w[seq(0, (GMMfit_Mode[[1]]$nx_w -
        1) * GMMfit_Mode[[1]]$prior$Kmax, by = GMMfit_Mode[[1]]$prior$Kmax) +
        c(1:Modes)]
      Mixtures <- cbind(MeansMCMC, SDsMCMC, WeightsMCMC)
      if (Modes == 1) {
        Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]), ]))
        Mixtures[3] <- 1
      }
    } else {
      GMMfit_Mode <- vector(mode = "list", length = 9)
      Mixtures <- Mixtures1
    }
    set.seed(ActualSeed)
  }, normalmixEM = {
    NewSeed <- ActualSeed + 1e+05
    FitError <- TRUE
    FitRetriesCount <- 0
    while (FitError == TRUE & FitRetriesCount < MaxRetries) {
      GMMfit_Mode <- try(mixtools::normalmixEM(GMMdata, mu = kmeans(GMMdata,
        Modes)$centers, ECM = TRUE, maxrestarts = 1e+05), TRUE)
      if (!is(GMMfit_Mode, "try-error")) FitError <- FALSE
      FitRetriesCount <- FitRetriesCount + 1
      if (is(GMMfit_Mode, "try-error")) set.seed(NewSeed + FitRetriesCount)
    }
    if (!is(GMMfit_Mode, "try-error")) {
      Mixtures <- cbind(GMMfit_Mode$mu, GMMfit_Mode$sigma, GMMfit_Mode$lambda)
      if (Modes == 1) {
        Mixtures <- t(as.matrix(Mixtures[which.max(Mixtures[, 3]), ]))
        Mixtures[3] <- 1
      }
    } else {
      GMMfit_Mode <- vector(mode = "list", length = 9)
      Mixtures <- Mixtures1
    }
    set.seed(ActualSeed)
  })
  return(list(GMMfit_Mode = GMMfit_Mode, Mixtures = Mixtures))
}
