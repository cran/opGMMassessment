# Functions for mode number determination methods Identify best fit based on
# selected criterion.
#' @importFrom mixtools normalmixEM
#' @importFrom NbClust NbClust
#' @importFrom AdaptGauss InformationCriteria4GMM LikelihoodRatio4Mixtures
#' @importFrom utils head tail sessionInfo
#' @importFrom stats ks.test rnorm sd
#' @importFrom cluster pam maxSE
#' @importFrom multimode modetest
DetermineBestGMM <- function(GMMdata, GMMfit, Criterion, ActualSeed, MaxModes, BestGMM,
                             nProc) {

  # Function to identify the best GMM based on AIC/BIC
  idBestGMM_AICBIC <- function(GMMdata, GMMfit, Criterion, ActualSeed) {
    AICBIC <- lapply(seq_along(GMMfit), function(x) {
      AICBICi <- AdaptGauss::InformationCriteria4GMM(Data = GMMdata, Means = as.matrix(GMMfit[[x]]$Mixtures)[, 1],
                                                     SDs = as.matrix(GMMfit[[x]]$Mixtures)[, 2], Weights = as.matrix(GMMfit[[x]]$Mixtures)[, 3])
      return(AICBICi)
    })
    AICBIC <- unlist(lapply(seq_along(GMMfit), function(x) {
      switch(Criterion, AIC = AICBIC[[x]]$AIC, BIC = AICBIC[[x]]$BIC)
    }))
    BestGMM <- which.min(AICBIC)
    return(BestGMM)
  }

  # Function to identify the best GMM based on Excess Mass (FM)
  idBestGMM_FM <- function(GMMdata, MaxModes, nProc) {
    list.of.Modes <- 1:(MaxModes - 1)
    ExcessMassi <- switch(Sys.info()[["sysname"]],
                          Windows = lapply(list.of.Modes, function(x) multimode::modetest(GMMdata, method = "FM", mod0 = x, B = 60)$p.value),
                          parallel::mclapply(list.of.Modes, function(x) multimode::modetest(GMMdata, method = "FM", mod0 = x, B = 60)$p.value, mc.cores = nProc))
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

  # Function to identify the best GMM based on Gap Statistic
  idBestGMM_GAP <- function(GMMdata, MaxModes, ActualSeed, nProc) {
    pam1 <- function(x, k) {
      list(cluster = cluster::pam(x, k = 3, metric = "euclidean", stand = FALSE, cluster.only = TRUE))
    }
    gsPam1 <- clusGapP(x = cbind(GMMdata, GMMdata), FUNcluster = pam1, K.max = MaxModes,
                       B = 60, spaceH0 = "original", nProc = nProc)
    BestGMM <- with(gsPam1, maxSE(Tab[, "gap"], Tab[, "SE.sim"], method = "globalSEmax"))
    return(BestGMM)
  }

  # Function to identify the best GMM based on Likelihood Ratio (LR)
  idBestGMM_LR <- function(GMMdata, GMMfit, ActualSeed) {
    LRi <- c(1, unlist(lapply(2:MaxModes, function(x) {
      AdaptGauss::LikelihoodRatio4Mixtures(Data = GMMdata, NullMixture = as.matrix(GMMfit[[x - 1]]$Mixtures),
                                           OneMixture = as.matrix(GMMfit[[x]]$Mixtures), PlotIt = FALSE)$Pvalue
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

  # Function to identify the best GMM based on NbClust
  idBestGMM_NbClust <- function(GMMdata, MaxModes) {
    NBres <- suppressWarnings(NbClust::NbClust(GMMdata, method = "kmeans", distance = "euclidean",
                                               max.nc = MaxModes), classes = "warning")
    BestGMM <- length(unique(NBres$Best.partition))
    return(BestGMM)
  }

  # Function to identify the best GMM based on Silverman's Test (SI)
  idBestGMM_SI <- function(GMMdata, MaxModes, nProc) {
    list.of.Modes <- 1:(MaxModes - 1)
    ExcessMassi <- switch(Sys.info()[["sysname"]],
                          Windows = lapply(list.of.Modes, function(x) multimode::modetest(GMMdata, method = "SI", mod0 = x, B = 60)$p.value),
                          parallel::mclapply(list.of.Modes, function(x) multimode::modetest(GMMdata, method = "SI", mod0 = x, B = 60)$p.value, mc.cores = nProc))
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

  # Determine the best GMM based on the selected criterion
  switch(Criterion,
         AIC = BestGMM <- idBestGMM_AICBIC(GMMdata = GMMdata, GMMfit = GMMfit, Criterion = Criterion,
                                           ActualSeed = ActualSeed),
         BIC = BestGMM <- idBestGMM_AICBIC(GMMdata = GMMdata, GMMfit = GMMfit, Criterion = Criterion,
                                           ActualSeed = ActualSeed),
         FM = BestGMM <- idBestGMM_FM(GMMdata = GMMdata, MaxModes = MaxModes, nProc = nProc),
         GAP = BestGMM <- idBestGMM_GAP(GMMdata = GMMdata, MaxModes = MaxModes, nProc = nProc,
                                        ActualSeed = ActualSeed),
         LR = BestGMM <- idBestGMM_LR(GMMdata = GMMdata, GMMfit = GMMfit, ActualSeed = ActualSeed),
         NbClust = BestGMM <- idBestGMM_NbClust(GMMdata = GMMdata, MaxModes = MaxModes),
         SI = BestGMM <- idBestGMM_SI(GMMdata = GMMdata, MaxModes = MaxModes, nProc = nProc),
         BestGMM = BestGMM)

  return(BestGMM)
}
