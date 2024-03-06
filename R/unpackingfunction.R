#' Run the transmission model
#'
#' @param theta named vector of model parameters; see description for details of elements
#' @param tx_pars named vector of model parameters relevant to treatment
#' @param runtime model run time, in years
#' @param stepsize model stepsize (dt) in years. Relatively small step sizes (1/12 or less) are recommended to avoid integration errors
#' @param user_tx_times a numeric vector of model treatment times (only used if input_tx_times==1 in tx_pars; else, can be safely set to NA); this parameterisation then ignores other treatment-related parameters specified in tx_pars
#' @param user_cov_weight a numeric vector of model coverage weights (only used if input_tx_times==1 in tx_pars; else, can be safely set to NA); this parameterisation then ignores other treatment-related parameters specified in tx_pars
#' @param time_extract_states a numeric model time at which to extract the state matrices. Can be safely set to NA if not needed
#' @description
#' This function allows the user to simulate an intestinal schistosomiasis transmission
#' model & return the (unlisted) output. See the package vignette for further details.
#' @return A list containing the output from the simulated transmission model
#' @export
RunTransmissionModel <- function(theta, tx_pars, runtime, stepsize, 
                                 user_tx_times, user_cov_weight, 
                                 time_extract_states) {
  
  # take log of parameter values (model takes the exponent)
  log_theta <- log(theta)
  
  # run model 
  out <-  RunModel(
    theta = log_theta, 
    tx_pars = tx_pars, 
    runtime = runtime, 
    stepsize = stepsize, 
    user_tx_times = user_tx_times,
    user_cov_weight = user_cov_weight, 
    # specify time to extract state matrices
    time_extract_states = time_extract_states
  )
  
  # unpack listed output from Rcpp model function
  unpacked_out <- list()
  
  for (i in 1:length(out)) {
    sublist <- out[[i]]
    unpacked_out <- c(unpacked_out, sublist)
  }
  
  return(unpacked_out)
}
