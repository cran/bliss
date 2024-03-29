% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic_functions.R
\name{BIC_model_choice}
\alias{BIC_model_choice}
\title{BIC_model_choice}
\usage{
BIC_model_choice(Ks, iter, data, verbose = T)
}
\arguments{
\item{Ks}{a numerical vector containing the K values.}

\item{iter}{an integer, the number of iteration for each run of \code{fit_Bliss}.}

\item{data}{a list containing:
\describe{
\item{Q}{an integer, the number of functional covariates.}
\item{y}{a numerical vector, the outcomes.}
\item{x}{a list of matrices, the qth matrix contains the observations of the
      qth functional covariate at time points given by \code{grids}.}
\item{grids}{a list of numerical vectors, the qth vector is the grid of
       time points for the qth functional covariate.}
}}

\item{verbose}{write stuff if TRUE (optional).}
}
\value{
A numerical vector, the BIC values for the Bliss model for different
        K value.
}
\description{
Model selection with BIC criterion.
}
\examples{
\donttest{
param_sim <- list(Q=1,n=100,p=c(50),grids_lim=list(c(0,1)))
data      <- sim(param_sim,verbose=TRUE)
iter = 1e2
Ks <- 1:5

res_BIC <- BIC_model_choice(Ks,iter,data)
plot(res_BIC,xlab="K",ylab="BIC")
}
}
