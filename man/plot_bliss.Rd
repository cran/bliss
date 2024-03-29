% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Auxiliary_and_graphics_functions.R
\name{plot_bliss}
\alias{plot_bliss}
\title{plot_bliss}
\usage{
plot_bliss(x, y, connect = FALSE, xlab = "", ylab = "", ylim = NULL, ...)
}
\arguments{
\item{x}{the coordinates of points in the plot.}

\item{y}{the y coordinates of points in the plot.}

\item{connect}{a logical value (optional), to handle discontinuous function.
If \code{connect} is TRUE, the plot is one line. Otherwise, several
lines are used.}

\item{xlab}{a title for the x axis.}

\item{ylab}{a title for the y axis.}

\item{ylim}{a numeric vectors of length 2, giving the y coordinate range.}

\item{...}{Arguments to be passed to methods, such as graphical parameters
(see \code{par}).}
}
\description{
A suitable representation of the Bliss estimate.
}
\examples{
\donttest{
data(data1)
data(param1)
# res_bliss1 <- fit_Bliss(data=data1,param=param1,verbose=TRUE)
}
data(res_bliss1)
### Plot the BLiss estimate on a suitable grid
plot_bliss(res_bliss1$data$grids[[1]],
           res_bliss1$Bliss_estimate[[1]],lwd=2,bound=FALSE)
}
