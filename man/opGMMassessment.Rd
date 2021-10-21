\name{opGMMassessment}
\alias{opGMMassessment}
\title{Gaussian mixture assessment}
\description{
  The package provides the necessary functions for optimized automated evaluation of the number and parameters of Gaussian mixtures in one-dimensional data.
  It provides various methods for parameter estimation and for determining the number of modes in the mixture.
}
\usage{
opGMMassessment(Data, FitAlg = "MCMC", Criterion = "LR",
MaxModes = 8, MaxCores = 2048, PlotIt = FALSE, KS = TRUE, Seed)
}
\arguments{
 \item{Data}{the data as a vector.}
  \item{FitAlg}{which fit algorithm to use: "ClusterRGMM" = GMM from ClusterR, 
  "densityMclust" from mclust,
  "DO" from DistributionOptimization (slow),
  "MCMC" = NMixMCMC from mixAK, or
  "normalmixEM" from mixtools.}
 \item{Criterion}{which criterion should be used to etsablish the number of modes from the best GMM fit:
 "AIC", "BIC", "FM", "GAP", "LR" (likelyhhod ratio test), "NbClust" (from NbClust), "SI" (Silverman).}
 \item{MaxModes}{the maximum number of modes to be tried.}
 \item{MaxCores}{the maximum number of processor cores used under Unix.}
 \item{PlotIt}{whether to plot the fit directly (plot will be stored nevertheless).}
 \item{KS}{perform a Kolmogorov-Smirnov test of the fit versus orginal distribution.}
 \item{Seed}{optional seed parameter set internally.}
}
\value{
  Returns a list of Gaussian modes.
 \item{Cls}{the classes to which the cases are assigned according to the Gaussian mode membership.}
 \item{Means}{means of the Gaussian modes.}
 \item{SDs}{standard deviations of the Gaussian modes.}
 \item{Weights}{weights of the Gaussian modes.}
 \item{Boundaries}{Bayesian boundaries between the Gaussian modes.}
 \item{Plot}{Plot of the obtained mixture.}
 \item{KS}{Results of the Kolmogorov-Smirnov test.}
 }
\references{
 Lotsch, J. (2021):
 opGMMassessment -- an R Package for automated Guassian mixture modeling.
}
\author{
  Jorn Lotsch and Alfred Ultsch
}
\examples{
## example 1
data(iris)
opGMMassessment(Data = iris$Petal.Length,
  FitAlg = "normalmixEM", 
  Criterion = "BIC",
  PlotIt = TRUE,
  MaxModes = 5,
  MaxCores = 1,
  Seed = 42)
}
\keyword{opGMMassessment}
\keyword{Clustering}
