#' Run PCLake management pathway with modified model parameters
#'
#' @param val_pars values of the pathway parameters
#' @param name_pars names of the pathway parameters
#' @param initial_conditions dataframe of initial values (from prepInitials() or use default)
#'
#' @returns a dataframe of PCLake output as defined in the DATM file
#' @export
#'
#' @examples

run_pathway <- function(val_pars, name_pars, initial_conditions = NULL) {
  
  # For debugging ----------------- #
  # val_pars <- lower_bound
  # name_pars <- names(lower_bound)
  # future_states <- desired_states #
  #--------------------------------#

  # Setting parameter values -----------#
  lDATM_SETTINGS_obj <- lDATM_SETTINGS
  lag_pars_lag <- grep('_lag', name_pars) # these are model parameters with lags - take these out
  lag_pars <- which(name_pars %in% gsub('_lag', '', name_pars[lag_pars_lag]))
  
  # it wasn't working if there were no lagged pars
  if (length(lag_pars) == 0) {
    none_lag_pars <- name_pars
  } else {
    none_lag_pars <- name_pars[-c(lag_pars, lag_pars_lag)]
  }
  
  # check if these non-lagged parameter values are listed in the forcings and change there if it is, if not just change the parameter value
  for (p in none_lag_pars) {
    if (p %in% names(lDATM_SETTINGS_obj$forcings$sDefault0)) {
      lDATM_SETTINGS_obj$forcings$sDefault0[[p]]$value <- val_pars[which(name_pars ==p)]
    } else {
      lDATM_SETTINGS_obj$params[name_pars[-c(lag_pars, lag_pars_lag)], "sDefault0"] <- val_pars[-c(lag_pars, lag_pars_lag)]
    }
  }
  
  
  # set the lagged variables via the forcings in the DATM file
  # extract lagged variable 
  # NOTE: this  only works if there is already something in the forcings and you are just modifying it
  # otherwise you have to recompile the model
  
  if (length(lag_pars_lag) > 0) {
    for (i in 1:length(lag_pars)) {
      
      name1 <- name_pars[lag_pars[i]]
      name2 <- name_pars[lag_pars_lag[i]]
      
      forcing_df <- data.frame(time = 0:(lDATM_SETTINGS_obj$run_settings["dReady", "Set0"]*365)) |> 
        mutate(value = current_val[name1]) |> # set with a unchanged/current loading)
        mutate(value = ifelse(time %in% 0:(val_pars[which(name_pars == name2)]*365),
                              value,
                              val_pars[which(name_pars == name1)])) # the values after the first lag are reduced to the parameter value
      
      lDATM_SETTINGS_obj$forcings$sDefault0[[name1]] <- forcing_df
      
    }
  }
  
  
  # Initialise and run model with these parameters
  
  if (!is_null(initial_conditions)) {  
    lDATM_SETTINGS_obj$states[initial_conditions$variable, "sDefaultSetTurbid0"] <- initial_conditions$value 
    # these could be extracted from a equilibrium model run
    message("Using provided initial conditions not DATM file values")
  }
  # or just reinitialise based on the DATM file

  InitStates_01 <- PCModelInitializeModel(lDATM = lDATM_SETTINGS_obj,
                                          dirSHELL = dirShell,
                                          nameWORKCASE = nameWorkCase)
  
  PCModel_run <- PCmodelSingleRun(lDATM = lDATM_SETTINGS_obj,
                                  nRUN_SET = 0,
                                  dfSTATES = InitStates_01,
                                  integrator_method = "rk45ck",
                                  dirHOME = dirHome,
                                  nameWORKCASE = nameWorkCase)
  return(PCModel_run)
  
}


#' Evaluate a parameter pathway
#'
#' @param PCLake_output output from run_pathway
#' @param future_states dataframe of variable and target value
#' @param eval_target how should the target be evaluated, list of functions to call, matched to the state names
#' @param eval_days which days in the last year should be evaluated
#'
#' @returns numeric value output from eval_target function call
#' @export
#'
#' @examples
evaluate_pathway <- function(PCLake_output, 
                             future_states,
                             eval_target = list(function(out,target){abs(out-target)/target}),
                             eval_days = 121:244) {
  
  # For debugging ----------------- #
  # PCLake_output = PCModel_run
  # future_states = desired_states
  # eval_target = list(function(x,y){(x-y)/y})
  # eval_days = 121:244
  #--------------------------------#
  
  if (sum(!future_states$variable %in% colnames(PCLake_output)) > 0) {
    stop("PCLake output doesn't contain the evaluation variables from the future states.")
  }
  
  if ( sum(!names(eval_target) %in% colnames(PCLake_output)) > 0) {
    stop("PCLake output doesn't contain the evaluation variables from the eval_target. Either remove from the obj_function or amend future_statees")
  }
  
  # Extract model output and compare with desired state
  model_output <- PCLake_output |>
    mutate(year = floor((time-1)/365) + 1,
           doy = yday(as_date(time - (year * 365) + 364, origin = '2025-01-01'))) |> 
    filter(year == max(year), # filters to summer in the last year of the simulation
           doy %in% eval_days) |> 
    select(future_states$variable) |> 
    summarise(across(any_of(future_states$variable), mean)) |> 
    pivot_longer(cols = any_of(future_states$variable),
                 names_to = 'variable',
                 values_to = 'output')
  
  # need a more complex process if the evaluation happens seperately (different function by objective)
  if (length(eval_target) == 1) {
    pathway_error <- model_output |> 
      full_join(future_states, by = 'variable') |> 
      mutate(diff = eval_target[[1]](output, target)) |> 
      summarise(total_error = sum(diff)) 
  } else {
    message('matching evaluation function by state')
    
    pathway_error <- model_output |> 
      full_join(future_states, by = 'variable') |> 
      mutate(diff= NA)
    
    for (i in 1:nrow(pathway_error)) {
      use_fun <- eval_target[[which(names(eval_target) == pathway_error$variable[i])]]
      pathway_error$diff[i] <- use_fun(pathway_error$output[i], pathway_error$target[i])
    }
    
    pathway_error <- pathway_error  |> 
      summarise(total_error = sum(diff)) ### IS THIS HOW YOU WOULD SUM THEM????
  
    
  }
  
  
  
  # print(model_output)
  
  out_val <- pathway_error |> pull(total_error)
  return(out_val)
  
}
