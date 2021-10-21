# Function to plot a Gaussian mixture based on the model's parameter values
#' @importFrom AdaptGauss BayesDecisionBoundaries
#' @importFrom grDevices nclass.FD
#' @importFrom methods hasArg
#' @importFrom stats density na.omit
#' @importFrom DataVisualizations ParetoDensityEstimation
#' @importFrom rlang .data
#' @import ggplot2
#' @export
GMMplotGG <- function(Data, Means, SDs, Weights, BayesBoundaries, SingleGausses = TRUE,
  Hist = FALSE, Bounds = TRUE, SumModes = TRUE, PDE = TRUE) {
  if (!hasArg("Means") | !hasArg("SDs") | !hasArg("Weights"))
    stop("GMMplotGG: Incomplete parameters.")
  if (SingleGausses == FALSE & SumModes == FALSE)
    stop("GMMplotGG: No model results requested to plot.\nUse other function to plot just raw data.")
  if (length(c(Means, SDs, Weights))%%3 != 0)
    stop("GMMplotGG: Unequal number of modes in parameters.")

  n <- 1000
  flag <- 0
  if (hasArg(Data) == TRUE) {
    Data <- na.omit(Data)
    if (length(Data) > 1) {
      if (PDE == TRUE) {
        PDF <- DataVisualizations::ParetoDensityEstimation(Data)
        x <- PDF$kernels
        PDF <- cbind.data.frame(Curve = "Data", x = PDF$kernels, y = PDF$paretoDensity)
      } else {
        PDF <- stats::density(Data)
        x <- PDF$x
        PDF <- cbind.data.frame(Curve = "Data", PDF$x, PDF$y)
      }
    } else flag <- 1
  } else flag <- 1
  if (flag == 1) {
    rangeX <- range(c(Means - 2 * SDs, Means + 2 * SDs))
    x <- seq(from = rangeX[1], to = rangeX[2], length.out = n)
    Hist = FALSE
    PDE = FALSE
  }

  Means0 <- Means[order(Means)]
  SDs0 <- SDs[order(Means)]
  Weights0 <- Weights[order(Means)]

  DataDF <- data.frame(mapply(function(w, mean, sd) w * dnorm(x, mean, sd), mean = Means0,
    sd = SDs0, w = Weights0))
  ModeMaxima <- apply(DataDF, 2, max)
  ModeMaxima <- ModeMaxima + 0.015 * max(ModeMaxima)
  SumPDF <- cbind.data.frame(Curve = "Sum", x = x, y = rowSums(DataDF))
  DataDF_longSingle <- cbind.data.frame(Curve = rep(paste0("M", c(1:length(Means0))),
    each = length(x)), x = rep(x, length(Means0)), y = as.vector(as.matrix(DataDF)))

  if (PDE == TRUE & SumModes == TRUE)
    DataDF_long <- rbind.data.frame(PDF, SumPDF)
  if (PDE == TRUE & SumModes == FALSE)
    DataDF_long <- rbind.data.frame(PDF)
  if (PDE == FALSE & SumModes == TRUE)
    DataDF_long <- rbind.data.frame(SumPDF)
  if (SingleGausses == TRUE)
    DataDF_long <- rbind.data.frame(DataDF_long, DataDF_longSingle)

  # Plot data
  requireNamespace("ggplot2")
  requireNamespace("grDevices")
  p1 <- ggplot()
  if (Hist == TRUE) {
    if (hasArg("Data") == TRUE & length(Data) > 1)
      HistData <- Data else HistData <- SumPDF
    breaks <- pretty(range(HistData), n = nclass.FD(HistData), min.n = 1)
    bwidth <- breaks[2] - breaks[1]
    p1 <- p1 + geom_histogram(data = data.frame(HistData), aes(x = HistData,
      .data$..density..), binwidth = bwidth, fill = "grey85", color = "grey95") +
      guides(fill = "none")
  }
  p1 <- p1 + geom_line(data = DataDF_long, aes(x = x, y = .data$y, colour = .data$Curve)) +
    labs(x = "Data", y = "Probability density")
  if (Bounds == TRUE & length(Means0) > 1) {
    flag <- 0
    if (hasArg(BayesBoundaries) == TRUE) {
      if (length(BayesBoundaries) == (length(Means) - 1)) {
        BayesBoundaries0 <- BayesBoundaries
      } else flag <- 1
    } else flag <- 1
    if (flag == 1)
      BayesBoundaries0 <- AdaptGauss::BayesDecisionBoundaries(Means = Means0,
        SDs = SDs0, Weights = Weights0)
    BayesBoundariesDF <- data.frame(BayesBoundaries0)
    p1 <- p1 + geom_vline(data = BayesBoundariesDF, aes(xintercept = BayesBoundaries0),
      color = "magenta", linetype = 2)
  }
  if (length(Means0) > 1)
    p1 <- p1 + annotate("text", x = Means0, y = ModeMaxima, size = rel(4), label = paste0("M#",
      1:length(Means0)))
  return(p1)
}
