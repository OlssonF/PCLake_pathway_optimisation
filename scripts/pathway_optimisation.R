#--------------------------------------#
## Project: Pathway Optimisation Framework
## Script purpose: Simple example of an optimsation framework using PCLake for a single objective (chlorophyll a < 20 ug/L) and two measures (P load reduction, marsh area) 
#                  including one with a temporal component
## Update: script now has options to optimise multiple targets simultaneously using > 2 measures, runs the final iteration population to extract the state values (not just the
#          objective function values) and plot as a time series and summarised values
## Date: 2025-09-05; updated 2025-10-01
## Author: Freya Olsson
#--------------------------------------#

library(tidyverse)
library(here)
library(DEoptim)
library(doSNOW)
library(parallelly)

## 0. Settings --------------------------
##---------------------------------------#

## Global settings
options(scipen = 999) ## no scientific notation

## 1. Directory settings ---------------------------------------------------------
## using relative paths in which the project and script is saved in the work_cases
## "scripts" contains only the PCLake functions

project_location <- here()
dirHome <- str_split(project_location,  "(?=PCModel1350)", simplify = T)[1,1]	# location of the PCModel1350 folder
dirShell <- str_split(project_location,  "(?<=PCShell)", simplify = T)[1,1]	#  PCShell folder path
dirCpp_root <- list.dirs(dirHome)[which(str_detect(list.dirs(dirHome),"3.01/PCLake_plus"))] # location of C++ code
nameWorkCase <- tail(str_split_1(project_location, "/"), n = 1) # workcase name
fileDATM <- list.files(list.dirs(dirHome)[which(str_detect(list.dirs(dirHome), "PCLake\\+/6.13.16"))], "PL613162PLUS_pathway_optim.xls", full.names = T)
folderTXT <- file.path(project_location, 'input', 'drivers_txt')
dirSave <- dirShell
# ----------------------------------------------------------------------------- #

## load external functions from the scripts folder
source(file.path(dirShell, "scripts", "R_system", "functions.R"))
source(file.path(dirShell, "scripts", "R_system", "functions_PCLake.R")) 


## Order of actions to run PCLake in R
##   1. Making folder structure for running the model
##   2. Load DATM file 
##   < Make adjustments to the model > 
##   4. Make cpp files
##   5. Compile model
##   6. Run model

## For the optimisation of pathways it is more like:
##   1. Making folder structure for running the model - this is likely already done as I am using a project structure
##   2. Load DATM file 
##   < Make adjustments to the model > 
##   4. Make cpp files
##   5. Compile model
##   6. Run optimisation - runs in parallel kind of in a forloop
#   a) dataframe of parameter ranges (model parameters) and lags
#   b) define objective function (comparison of desired future and "current" state)
#   c) identify top pathways

## 2. Load DATM file  -------------------             
lDATM_SETTINGS <- PCModelReadDATMFile_PCLakePlus(fileXLS = fileDATM,
                                                 folderTXT = folderTXT,
                                                 locDATM = "excel",
                                                 locFORCING = "txt",
                                                 readAllForcings = F)
##----------------------------------------#

## Modifications can be made to the DATM file here! (those that need to be compiled)
## Might be a good idea to make sure the default lake parameters are loaded

# Report restart variables
restart_states <- read_table(file.path(project_location, 'restart_states.txt'), col_names = 'state', show_col_types = F)
lDATM_SETTINGS$auxils$iReport[which(rownames(lDATM_SETTINGS$auxils) %in% restart_states$state)] <- 1 # report these in the output


## 3.4.A Make and adjust cpp files ----------       
#  - nRUN_SET determines which forcings are switched on
PCModelAdjustCPPfiles(dirSHELL = dirShell,
                      nameWORKCASE = nameWorkCase,
                      lDATM = lDATM_SETTINGS,
                      nRUN_SET = 0)
##----------------------------------------#

## 5.A Compile model -----------------------
PCModelCompileModelWorkCase(dirSHELL = dirShell,
                            nameWORKCASE = nameWorkCase)
##----------------------------------------#

## Find equilibrium values ----------------#
# Run the model to an equilibrium before starting the optimisation
InitStates_baseline <- PCModelInitializeModel(lDATM = lDATM_SETTINGS,
                                              dirSHELL = dirShell,
                                              nameWORKCASE = nameWorkCase)

PCModel_run_baseline <- PCmodelSingleRun(lDATM = lDATM_SETTINGS,
                                         nRUN_SET = 0,
                                         dfSTATES = InitStates_baseline,
                                         integrator_method = "rk45ck",
                                         dirHOME = dirHome,
                                         nameWORKCASE = nameWorkCase)


# extract the restart variables from the end of the baseline run
equilibrium_states <- prepInitials(listPCModelRun = PCModel_run_baseline, 
                                   day =  lDATM_SETTINGS$run_settings['dReady','Set0'] * 365)


## 6.A Set optimisation values ---------------------------

### a. Define parameter values -------------------
# dataframe of parameter ranges (model parameters) and lags 

# Define the parameter values to be optimised
lower_bound <- c('mPLoadEpi' = 0.001, #minimum Ploading
                 #'mPLoadEpi_lag' = 1,
                 'fMarsh_lag' = 5, # need at least a year?
                 'fMarsh' = 0) # fraction of marsh area relative to lake
upper_bound <- c('mPLoadEpi' = 0.01,
                 #'mPLoadEpi_lag' = 10,
                 'fMarsh_lag' = 20*1, # could be at least 20 year lag
                 'fMarsh' = 0.2) 

# error_pars <- c('mPLoadEpi' = 0.0005531256, #minimum Ploading
#                 'mPLoadEpi_lag' = 4.3935696417, # need at least a year?
#                 'fMarsh' = 0.8961858496)

# these are needed to give a value for any variable that is lagged
current_val <- c('mPLoadEpi' = 0.05,
                 'fMarsh' = 0) 
# the "unchanged" value (before the measure is in place) - could also be a timeseries I guess?

### b. Define the desired future ------------------
# What is the objective
# Define the desired future state(s)
desired_states <- data.frame(variable = c('oChlaEpi', 'aDFish'),#, 'aSecchiT'),
                             target = c(20, 6))#, 0.5))


## Update the DATM file and recompile the model ---------------#
# Report variables
lDATM_SETTINGS$auxils$iReport[which(rownames(lDATM_SETTINGS$auxils) %in% restart_states$state)] <- 0 # these can be turned off
lDATM_SETTINGS$auxils$iReport[which(rownames(lDATM_SETTINGS$auxils) %in% desired_states$variable)] <- 1 # report optim vars

# forcing variables --------------#
# anything that is being lagged needs to be in the forcings before compilation
for (i in names(lower_bound)) {
  lag_var <- gsub('_lag', '', i)
  if (str_detect(i, '_lag') & is.null(lDATM_SETTINGS$forcings$sDefault0[[lag_var]])) {
    
    
    lDATM_SETTINGS$forcings$sDefault0[[lag_var]] <- data.frame(time = 0:(lDATM_SETTINGS$run_settings["dReady", "Set0"]*365),
                                                               value = -99999) 
    message('adding ', lag_var, ' forcing')
  }
}


## 3.4. B Remake and adjust cpp files ----------       
#  - nRUN_SET determines which forcings are switched on
PCModelAdjustCPPfiles(dirSHELL = dirShell,
                      nameWORKCASE = nameWorkCase,
                      lDATM = lDATM_SETTINGS,
                      nRUN_SET = 0)
##----------------------------------------#

## 5.B Compile model -----------------------
PCModelCompileModelWorkCase(dirSHELL = dirShell,
                            nameWORKCASE = nameWorkCase)
##----------------------------------------#

# The copiled model needs to have all of the forcing variables in it (e.g. marsh area, P-loading)
# otherwise the .dll and DATM file won't match



### c. Define the objective function -----------------
# 1) run model
# 2) extract output
# 3) return a minimised value by comparing with desired future (DEOptim will make the value as negative as possible)

source(file.path(project_location, "scripts/optim_functions.R")) # functions for running and evaluating the pathways

obj_function <- function(val_pars, name_pars, future_states) {
  
  # For debugging ----------------- #
  # val_pars <- upper_bound
  # name_pars <- names(lower_bound)
  # future_states <- desired_states #
  #--------------------------------#
  
  model_output <- run_pathway(val_pars, name_pars)
  
  eval_output <- evaluate_pathway(PCLake_output = model_output, 
                                  future_states = future_states,
                                  eval_target = list(oChlaEpi = function(out,target){(out-target)/target},   # below better
                                                     aDFish = function(out,target){abs(out-target)/target},  # target exact value
                                                            # function(out,target){(target-out)/target}      # above better
                                                     ))
  return(eval_output)
}


## 6.B Run optimisation ---------------

## Parallelisation with FOREACH package
{
  nC <- parallelly::availableCores() 
  cl <- makeSOCKcluster(nC-2)
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(deSolve))
  clusterExport(cl, list("lDATM_SETTINGS", 'current_val', 'equilibrium_states',
                         "PCModelInitializeModel", 
                         "dirShell", "nameWorkCase", 'dirHome',
                         "PCmodelSingleRun", "RunModel", 'run_pathway', 'evaluate_pathway'))
  
  doSNOW::registerDoSNOW(cl)
  
  deoptim_control <- list(NP = 15 * length(lower_bound),
                          CR = 0.9,
                          F = 0.80, 
                          trace = TRUE,
                          itermax = 20, 
                          reltol = 0.001, 
                          steptol = 5, 
                          #c = 0.05, strategy = 6, p = 0.10,
                          c = 0.25, strategy = 6, p = 0.40, 
                          storepopfrom = 0,
                          parallelType = 'foreach')
  
  
  set.seed(1234)
  
  opt_pathway <- DEoptim::DEoptim(lower = lower_bound,
                                  upper = upper_bound, 
                                  fn = obj_function,
                                  control = deoptim_control, 
                                  name_pars = names(lower_bound),
                                  future_states = desired_states)
  
  parallel::stopCluster(cl)
}

# The output of DEoptim is based on members, iterations, and populations
# iteration is a generation of a population
# one population is a collection of members
# one member is a single run of the objective function with specific parameter values

## Visualise the best member output from each iteration
n_iter <- opt_pathway$optim$iter
iteration_summary <- as.data.frame(opt_pathway$member$bestmemit) |>
  mutate(fn_out = opt_pathway$member$bestvalit[1:n_iter])


ggplot(iteration_summary, aes(x=fMarsh, y= fMarsh_lag, size = mPLoadEpi, colour = fn_out)) +
  geom_point() +
  scale_color_viridis_c() 

## 7. Re-run last iteration -----------------
### 7.a Exract objective values -----------
# The DEoptim does not output the results of the function (valit) only the best member per iteration and then the parameter values
# To get the output values we can rerun the last iteration of the optimisation
last_iteration <- as_tibble(opt_pathway$member$pop)
colnames(last_iteration) <- colnames(opt_pathway$member$bestmemit)

last_iteration$runID <- rownames(last_iteration)

# set up parallelisation
n_threads <- parallel::detectCores() - 2 ## leave one free for other tasks
snowCLUSTER <- makeCluster(n_threads)
# clusterExport(snowCLUSTER, c()) ## here you can push certain objects to the clusters
registerDoSNOW(snowCLUSTER)

# stuff to show progress bar
pb <- txtProgressBar(0, nrow(last_iteration), style = 3)
progress <-function(n){setTxtProgressBar(pb, n)}
opts <- list(progress = progress)

# Use foreach to run each member of the last iteration in the function
last_iteration$fn_out <- foreach(i = 1:nrow(last_iteration),
                                 .combine = 'c', 
                                 .multicombine = TRUE,
                                 .packages=c("tidyverse", "deSolve"),
                                 .options.snow = opts) %dopar% {
                                   
                                   val_pars <- last_iteration[i,] |> 
                                     select(-runID) |> unlist()
                                   
                                   name_pars <- last_iteration[i,] |> 
                                     select(-runID) |> names()
                                   
                                   save <- obj_function(val_pars = val_pars, name_pars = name_pars, future_states = desired_states)
                                 }

last_iteration |> 
  # slice_min(fn_out, prop = 0.25) |> # best (lowest) 25% of values
  ggplot(aes(x=fMarsh_lag, y= fMarsh, size = mPLoadEpi, colour = fn_out)) + 
  geom_point() + 
  scale_colour_viridis_c() 


### 7.b Extract the state values -----------
# what are the values of the optimised variables?
state_opt <- foreach(i = 1:nrow(last_iteration),
                     .combine = bind_rows, 
                     .multicombine = TRUE,
                     .packages=c("tidyverse", "deSolve"),
                     .options.snow = opts) %dopar% {
                       
                       val_pars <- last_iteration[i,names(lower_bound)] |> unlist()
                       
                       name_pars <- last_iteration[i,names(lower_bound)] |> names()
                       
                       df_pars <- data.frame(variable = name_pars, output = val_pars)
                       
                       run_pathway(val_pars, name_pars) |>
                         mutate(year = floor((time-1)/365) + 1,
                                doy = yday(as_date(time - (year * 365) + 364, origin = '2025-01-01'))) |> 
                         filter(year == max(year), # filters to summer in the last year of the simulation
                                doy %in% 121:244) |> 
                         select(desired_states$variable) |> 
                         summarise(across(any_of(desired_states$variable), mean)) |> 
                         pivot_longer(cols = any_of(desired_states$variable),
                                      names_to = 'variable',
                                      values_to = 'output') |> 
                         bind_rows(df_pars) |> 
                         mutate(ID = i)
                       
                     }


state_opt |> 
  pivot_wider(id_cols = ID, names_from = variable, values_from = output) |> 
  pivot_longer(cols = desired_states$variable, names_to = 'opt_var', values_to = 'out') |> 
  ggplot(aes(x=mPLoadEpi, y = out, size = fMarsh_lag, colour = fMarsh)) + geom_point() +
  facet_wrap(~opt_var, scales = 'free') +
  scale_color_viridis_c(option  = 'A', begin = 1, end = 0)

### 7.c Extract the time series ---------
# above we extracted the last year (the year being optimised) but in this case we want to 
# look at the full pathway to acheive the final target

# what is the "pathway" to the optimised value?
state_pathways <- foreach(i = 1:nrow(last_iteration),
                     .combine = bind_rows, 
                     .multicombine = TRUE,
                     .packages=c("tidyverse", "deSolve"),
                     .options.snow = opts) %dopar% {
                       
                       val_pars <- last_iteration[i,names(lower_bound)] |> unlist()
                       
                       name_pars <- last_iteration[i,names(lower_bound)] |> names()
                       
                       df_pars <- data.frame(variable = name_pars, output = val_pars)
                       
                       run_pathway(val_pars, name_pars) |>
                         mutate(year = floor((time-1)/365) + 1,
                                doy = yday(as_date(time - (year * 365) + 364, origin = '2025-01-01'))) |> 
                         filter(# filters to summer all years of the simulation
                                doy %in% 121:244) |> 
                         select(c('year', desired_states$variable)) |> 
                         group_by(year) |> 
                         summarise(across(any_of(desired_states$variable), mean)) |> 
                         bind_cols(pivot_wider(df_pars, names_from = variable, values_from = output)) |>  
                         mutate(ID = i)
                       
                     }
state_pathways |> 
  pivot_longer(cols = desired_states$variable, names_to = 'variable', values_to = 'state_val') |> 
  ggplot(aes(y=state_val, x=year, group = ID, colour = fMarsh)) +
  facet_wrap(~variable, scales = 'free', nrow=2)+
  # geom_vline(aes(xintercept = fMarsh_lag, colour =fMarsh), alpha = 0.7) +
  geom_line(aes(colour = fMarsh)) +
  scale_colour_viridis_c(option = 'A')

