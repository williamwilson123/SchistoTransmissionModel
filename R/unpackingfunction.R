#' Run the transmission model
#'
#' @param theta named vector of model parameters; see description for details of elements
#' @param tx_pars named vector of model parameters relevant to treatment
#' @param runtime model run time, in years
#' @param stepsize model stepsize (dt) in years. Relatively small step sizes (1/12 or less) are recommended to avoid integration errors
#' @param tx_times a numeric vector of model treatment times (only used if input_tx_times==1 in tx_pars)
#' @description
#' This function allows the user to simulate a Schistosoma mansoni transmission
#' model. Parameters can be specified in the 'theta' vector, comprised of the
#' following named elements (with sensible ranges of vales in parentheses)
#' - R0 : the basic reproduction number
#' - R0_weight : a weighting parameter that controls relative transmission from snails to humans vs. humans to snails
#' - NhNs: human-to-snail population density (0 to 1)
#' - kW : Negative binomial distribution shape parameter; controls aggregation of worms among hosts (0 to 1)
#' - decay_immunity : rate of decay of IgE antibodies
#' - protection_immunity : extent of protection provided by IgE antibodies (0 to 1)
#' - epg_constant : no. eggs shed per fertile adult female schistosome
#' - ige_constant : constant to map acquired immunity to IgE optical density (0 to 1)
#' @return A list containing the output from the simulated transmission model
#' @export
RunTransmissionModel <- function(theta, tx_pars, runtime, stepsize, tx_times) {
  
  # run model 
  out <-  RunModel(
    theta = theta, tx_pars = tx_pars, runtime = runtime, stepsize = stepsize, 
    tx_times = tx_times
  )
  
  # unpack listed output
  unpacked_out <- list()
  
  for (i in 1:length(out)) {
    sublist <- out[[i]]
    unpacked_out <- c(unpacked_out, sublist)
  }
  
  return(unpacked_out)
}
