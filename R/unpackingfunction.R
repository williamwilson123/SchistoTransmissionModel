#' Run the transmission model
#'
#' @param theta named vector of model parameters; see description for details of elements
#' @param tx_pars named vector of model parameters relevant to treatment
#' @param runtime model run time, in years
#' @param stepsize model stepsize (dt) in years. Relatively small step sizes (1/12 or less) are recommended to avoid integration errors
#' @param tx_times a numeric vector of model treatment times (only used if input_tx_times==1 in tx_pars)
#' @description
#' This function allows the user to simulate an intestinal schistosomiasis transmission
#' model. See the package vignette for further details.
#' @return A list containing the output from the simulated transmission model
#' @export
RunTransmissionModel <- function(theta, tx_pars, runtime, stepsize, tx_times) {
  
  # run model 
  out <-  RunModel(
    # take log of parameter values (model takes the exponent)
    theta = log(theta), 
    tx_pars = tx_pars, 
    runtime = runtime, 
    stepsize = stepsize, 
    tx_times = tx_times
  )
  
  # unpack listed output from Rcpp model function
  unpacked_out <- list()
  
  for (i in 1:length(out)) {
    sublist <- out[[i]]
    unpacked_out <- c(unpacked_out, sublist)
  }
  
  return(unpacked_out)
}
