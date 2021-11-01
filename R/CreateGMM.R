# Function to create a probability matrix based on multimodal Gaussian
# distribution
#' @importFrom methods hasArg
#' @importFrom stats dnorm
#' @importFrom caTools sample.split
#' @importFrom dplyr sample_n
CreateGMM <- function(Means, SDs, Weights, n = 1000, ExactN = FALSE) {
  if (!hasArg("Means") | !hasArg("SDs") | !hasArg("Weights")) {
    stop("CreateGMM: Incomplete parameters.")
  }
  if (length(c(Means, SDs, Weights)) %% 3 != 0) {
    stop("CreateGMM: Unequal number of modes in parameters.")
  }
  sumWeights <- sum(Weights)
  if (sumWeights != 1) {
    Weights <- Weights / sumWeights
  }

  if (ExactN == TRUE) {
    GMMparam <- rbind(round(Weights * 2 * n), Means, SDs)
  } else {
    GMMparam <- rbind(round(Weights * n), Means, SDs)
  }

  GMMparamRO <- split(GMMparam, rep(1:ncol(GMMparam), each = nrow(GMMparam)))
  DataDF <- cbind.data.frame(Data = as.vector(unlist(lapply(GMMparamRO, function(x) {
    do.call(rnorm, as.list(x))
  }))), Cls = rep(1:length(Weights), unlist(lapply(GMMparamRO, "[[", 1))))

  if (ExactN == TRUE) {
    sample <- caTools::sample.split(DataDF$Cls, SplitRatio = 0.51)
    DataDF <- dplyr::sample_n(subset(DataDF, sample == TRUE), 1000)
  }

  return(DataDF)
}