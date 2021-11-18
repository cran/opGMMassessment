\name{GMMplotGG}
\alias{GMMplotGG}
\title{Plot of Gaussian mixtures}
\description{
  The function plots the components of a Gaussian mixture and superimposes them on a histogram of the data.
}
\usage{
GMMplotGG(Data, Means, SDs, Weights, BayesBoundaries, 
	SingleGausses = TRUE, Hist = FALSE, Bounds = TRUE, SumModes = TRUE, PDE = TRUE)
}
\arguments{
 \item{Data}{the data as a vector.}
 \item{Means}{a list of mean values for a Gaussian mixture.}
 \item{SDs}{a list of standard deviations for a Gaussian mixture.}
 \item{Weights}{a list of weights for a Gaussian mixture.}
 \item{BayesBoundaries}{a list of Bayesian boundaries for a Gaussian mixture.}
 \item{SingleGausses}{whether to plot the single Gaussian components as separate lines.}
 \item{Hist}{whether to plot a histgram of the original data.}
 \item{Bounds}{whether to plot the Bayesian boundaries for a Gaussian mixture as vertical lines.}
 \item{SumModes}{whether to plot the summed-up mixes.}
 \item{PDE}{whether to use the Pareto density estimation instead of the standard R density function.}
}
\value{
  Returns a ggplot2 object.
 \item{p1}{the plot of Gaussian mixtures.}
 }
\references{
 Lotsch, J., Malkusch S. (2021):
 opGMMassessment -- an R Package for automated Guassian mixture modeling.
}
\author{
  Jorn Lotsch and Sebastian Malkusch
}
\examples{
## example 1
data(iris)
Means0 <- tapply(X = as.vector(iris[,3]), INDEX =  as.integer(iris$Species), FUN = mean)
SDs0 <- tapply(X = as.vector(iris[,3]), INDEX =  as.integer(iris$Species), FUN = sd)
Weights0 <- c(1/3, 1/3, 1/3)
GMM.Sepal.Length <- GMMplotGG(Data = as.vector(iris[3]), 
	Means = Means0, 
	SDs = SDs0, 
	Weights = Weights0, 
	Hist = TRUE) 
}
\keyword{GMMplotGG}
\keyword{data visualization}
