# Function to create a probability matrix based on multimodal Gaussian
# distribution
#' @importFrom methods hasArg
#' @importFrom stats dnorm
#' @importFrom caTools sample.split
#' @importFrom dplyr sample_n
CreateGMM <- function(Means, SDs, Weights, n = 1000, ExactN = FALSE) {
  # Check if all required parameters are provided
  if (!hasArg("Means") || !hasArg("SDs") || !hasArg("Weights")) {
    stop("CreateGMM: Incomplete parameters.")
  }

  # Check if the number of parameters is divisible by 3
  if (length(c(Means, SDs, Weights)) %% 3 != 0) {
    stop("CreateGMM: Unequal number of modes in parameters.")
  }

  # Normalize the weights if necessary
  sumWeights <- sum(Weights)
  if (sumWeights != 1) {
    Weights <- Weights / sumWeights
  }

  # Determine the number of samples per mode
  if (ExactN == TRUE) {
    GMMparam <- rbind(round(Weights * 2 * n), Means, SDs)
  } else {
    GMMparam <- rbind(round(Weights * n), Means, SDs)
  }

  # Generate the data
  GMMparamRO <- split(GMMparam, rep(seq_len(ncol(GMMparam)), each = nrow(GMMparam)))
  DataDF <- cbind.data.frame(
    Data = as.vector(unlist(lapply(GMMparamRO, function(x) {
      do.call(stats::rnorm, as.list(x))
    }))),
    Cls = rep(seq_along(Weights), unlist(lapply(GMMparamRO, "[[", 1)))
  )

  # Ensure the exact number of samples if requested
  if (ExactN == TRUE) {
    sample <- caTools::sample.split(DataDF$Cls, SplitRatio = 0.51)
    DataDF <- dplyr::sample_n(subset(DataDF, sample == TRUE), 1000)
  }

  return(DataDF)
}
