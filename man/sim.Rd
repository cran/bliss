% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simulate_Data.R
\name{sim}
\alias{sim}
\title{sim}
\usage{
sim(param, verbose = FALSE)
}
\arguments{
\item{param}{a list containing:
\describe{
 \item{beta_shapes}{a character vector. The qth item indicates the shape of
       the coefficient function associated to the qth functional covariate.}
 \item{beta_functions}{a list containing numerical vectors to define the beta
       functions}
 \item{n}{an integer, the sample size.}
 \item{p}{a vector of integers, the qth component is the number of
       times for the qth covariate.}
 \item{Q}{an integer, the number of functional covariates.}
 \item{autocorr_diag}{a list of numerical vectors (optional), the qth vector is the
       diagonal of the autocorrelation matrix of the qth functional
       covariate.}
 \item{autocorr_spread}{a vector of numerical values (optional) which are related to the
       autocorrelation of the functional covariates.}
 \item{grids}{a list of numerical vectors (optional), the qth vector is the grid
       of time points for the qth functional covariate.}
 \item{grids_lim}{a list of numerical vectors  (optional), the qth item is the lower
       and upper boundaries of the domain for the qth functional covariate.}
 \item{link}{a function (optional) to simulate data from the Generalized Functional
       Linear Regression model.}
 \item{mu}{a numerical value (optional), the 'true' intercept of the model.}
 \item{r}{a nonnegative value (optional), the signal to noise ratio.}
 \item{x_shapes}{a character vector (optional). The qth item indicates the shape of the
       functional covariate observations.}
}}

\item{verbose}{write stuff if TRUE.}
}
\value{
a list containing:
\describe{
 \item{Q}{an integer, the number of functional covariates.}
 \item{y}{a numerical vector, the outcome observations.}
\item{x}{a list of matrices, the qth matrix contains the observations of the
      qth functional covariate at time points given by \code{grids}.}
 \item{grids}{a list of numerical vectors, the qth vector is the grid of
       time points for the qth functional covariate.}
 \item{betas}{a list of numerical vectors, the qth vector is the 'true' coefficient
       function associated to the qth covariate on a grid of time points
       given with \code{grids}.}
}
}
\description{
Simulate a dataset for the Function Linear Regression model.
}
\examples{
library(RColorBrewer)
param <- list(Q=2,n=25,p=c(50,50),grids_lim=list(c(0,1),c(-1,2)))
data <- sim(param)
data$y
cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(10)
q=2
matplot(data$grids[[q]],t(data$x[[q]]),type="l",lty=1,col=cols)
plot(data$grids[[q]],data$betas[[q]],type="l")
abline(h=0,lty=2,col="gray")
}
