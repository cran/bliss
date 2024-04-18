################################# ----
#' fit_Bliss
################################# ----
#' @description Fit the Bayesian Functional
#' Linear Regression model (with Q functional covariates).
#' @return return a list containing:
#' \describe{
#'  \item{alpha}{a list of Q numerical vector. Each vector is the function
#'        alpha(t) associated to a functional covariate. For each t, alpha(t)
#'        is the posterior probabilities of the event "the support covers t".}
#'  \item{beta_posterior_density}{a list of Q items. Each item contains a list
#'        containing information to plot the posterior density of the
#'        coefficient function with the \code{image} function.
#'        \describe{
#'        \item{\code{grid_t}}{a numerical vector: the x-axis.}
#'        \item{\code{grid_beta_t}}{a numerical vector: the y-axis.}
#'        \item{\code{density}}{a matrix: the z values.}
#'        \item{\code{new_beta_sample}}{a matrix: beta sample used to compute
#'              the posterior densities.}
#'        }
#'        }
#'  \item{beta_sample}{a list of Q matrices. The qth matrix is a posterior
#'        sample of the qth functional covariates.}
#'  \item{Bliss_estimate}{a list of numerical vectors corresponding to the
#'        Bliss estimates of each functional covariates.}
#'  \item{chains}{a list of \code{posterior_sample}. \code{chains} is \code{NULL} if
#'        \code{n_chains}=1.}
#'  \item{chains_info}{a list for each chain providing: a mu estimate, a sigma_sq estimate,
#'  the Smooth estimate of the coefficient function and the autocorrelation of the
#'  Markov Chain.}
#'  \item{data}{a list containing the data.}
#'  \item{posterior_sample}{a list of information about the posterior sample:
#'        the trace matrix of the Gibbs sampler, a list of Gibbs sampler parameters
#'        and the posterior densities.}
#'  \item{support_estimate}{a list of support estimates of each functional covariate.}
#'  \item{support_estimate_fct}{another version of the support estimates.}
#'  \item{trace_sann}{a list of Q matrices which are the trace of the
#'        Simulated Annealing algorithm.}
#' }
#' @param data a list containing:
#' \describe{
#' \item{Q}{an integer, the number of functional covariates.}
#' \item{y}{a numerical vector, the outcomes.}
#' \item{x}{a list of matrices, the qth matrix contains the observations of the
#'       qth functional covariate at time points given by \code{grids}.}
#' \item{grids}{a list of numerical vectors, the qth vector is the grid of
#'        time points for the qth functional covariate.}
#' }
#' @param param a list containing:
#' \describe{
#' \item{iter}{an integer, the number of iterations of the Gibbs sampler algorithm.}
#' \item{K}{a vector of integers, corresponding to the numbers of intervals for
#'       each covariate.}
#' \item{basis}{a character vector (optional). The possible values are "uniform" (default),
#'       "epanechnikov", "gauss" and "triangular" which correspond to
#'       different basis functions to expand the coefficient function and the
#'       functional covariates}
#' \item{burnin}{an integer (optional), the number of iteration to drop from the
#'       posterior sample.}
#' \item{iter_sann}{an integer (optional), the number of iteration of the Simulated
#'       Annealing algorithm.}
#' \item{k_max}{an integer (optional), the maximal number of intervals for the
#'       Simulated Annealing algorithm.}
#' \item{l_max}{an integer (optional), the maximal interval length for the
#'       Simulated Annealing algorithm.}
#' \item{lims_kde}{an integer (optional), correspond to the \code{lims} option
#'       of the \code{kde2d} funtion.}
#' \item{n_chains}{an integer (optional) which corresponds to the number of
#'       Gibbs sampler runs.}
#' \item{new_grids}{a list of Q vectors (optional) to compute beta samples on
#'       different grids.}
#' \item{Temp_init}{a nonnegative value (optional), the initial temperature for
#'      the cooling function of the Simulated Annealing algorithm.}
#' \item{thin}{an integer (optional) to thin the posterior sample.}
#' \item{times_sann}{an integer (optional), the number of times the algorithm
#'       will be executed}
#' }
#' @param support_estimate a logical value. If TRUE, the estimate of the
#' coefficient function support is computed. (optional)
#' @param verbose write stuff if TRUE (optional).
#' @importFrom utils tail
#' @export
#' @examples
#' # see the vignette BlissIntro.
prior_Bliss <- function(data,param,support_estimate=TRUE,
                      verbose=FALSE){
  # Define Q
  Q <- data[["Q"]]
  if(is.null(Q))
    stop("Please specify Q: the number of functional covariates (in the 'data' object).")
  # Define p
  param$p  <- sapply(data$grids,length)
  # Centering the data
  data$x_save <- data$x
  for(q in 1:Q){
    data$x[[q]] <- scale(data$x[[q]],scale=F)
  }

  # How many chains i have to do ?
  if(is.null(param[["n_chains"]])){
    param[["n_chains"]] <- 1
  }
  n_chains <- param[["n_chains"]]
  # Initialize the list "chains"
  chains <- list()
  chains_info <- list()

  if(verbose) cat("Sample from the prior distribution.\n")
  # For each chain :
  for(j in 1:n_chains){
    if(verbose & n_chains > 1) cat("Chain ",j,": \n",sep="")
    chains[[j]] <- list()

    # Execute the Gibbs Sampler algorithm to sample the posterior distribution
    param_Gibbs_Sampler <- list(iter  = param[["iter"]],
                                K     = param[["K"]],
                                basis = param[["basis"]],
                                p     = param[["p"]],
                                phi_l = param[["phi_l"]],
                                grids = data[["grids"]])
    chains[[j]] <- Bliss_Gibbs_Sampler(data,param_Gibbs_Sampler,verbose)
    # ,to_sample="prior") # PMG 2024-04-15

    chains_info[[j]] <- compute_chains_info(chains[[j]],param_Gibbs_Sampler)
  }

  # Choose a chain for inference
  j <- sample(n_chains,1)
  prior_sample <- chains[[j]]

  # Compute a posterior sample of coefficient function
  if(verbose) cat("Coefficient function sample.\n")
  beta_sample <- compute_beta_sample(prior_sample,param_Gibbs_Sampler,Q)

  # Compute the support estimate
  if(support_estimate){
    if(verbose) cat("Support estimation.\n")
    support_estimate <- list()
    support_estimate_fct <- list()
    alpha <- list()
    for(q in 1:Q){
      res_support <- support_estimation(beta_sample[[q]])
      support_estimate[[q]]     <- res_support$estimate
      support_estimate_fct[[q]] <- res_support$estimate_fct
      alpha[[q]]                <- res_support$alpha
    }
    rm(res_support)
  }else{
    support_estimate <- list()
    support_estimate_fct <- list()
    alpha <- list()
  }

  # Do not return the list "chains" if n_chains is 1.
  if(n_chains == 1) chains <- NULL

  # The object to return
  res <- list(alpha                  = alpha,
              beta_sample            = beta_sample,
              chains                 = chains,
              chains_info            = chains_info,
              prior_sample           = prior_sample,
              support_estimate       = support_estimate,
              support_estimate_fct   = support_estimate_fct
  )
  class(res) = c("bliss")
  return(invisible(res))
}
