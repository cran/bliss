% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic_functions.R
\name{post_treatment_bliss}
\alias{post_treatment_bliss}
\title{post_treatment_bliss}
\usage{
post_treatment_bliss(posterior_sample, param, data)
}
\arguments{
\item{posterior_sample}{a list provided by the function \code{Bliss_Gibbs_Sampler}.}

\item{param}{a list containing:
\describe{
\item{K}{a vector of integers, corresponding to the numbers of intervals for
      each covariate.}
}}

\item{data}{a list containing required options to run the function
\code{dposterior}.}
}
\value{
A list of important post treatment value: BIC, the maximum of the log
likelihood and the numbre of parameters.
}
\description{
Compute the post treatment values.
}
\examples{
\donttest{
data(data1)
data(param1)
data(res_bliss1)

post_treatment_bliss(res_bliss1$posterior_sample,param1,data1)
}
}
