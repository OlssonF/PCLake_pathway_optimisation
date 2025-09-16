#--------------------------------------#
## Project: Pathway Optimisation Framework
## Script purpose: Simple example of an optimsation framework using PCLake for a single objective (chlorophyll a < 20 ug/L) and two measures (P load reduction, marsh area) 
#                  including one with a temporal component
## Date: 2025-09-05
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

## Directory settings ---------------------------------------------------------
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

## For the optimisation of pathways it is:
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


## 3.4. Make and adjust cpp files ----------       
#  - nRUN_SET determines which forcings are switched on
PCModelAdjustCPPfiles(dirSHELL = dirShell,
                      nameWORKCASE = nameWorkCase,
                      lDATM = lDATM_SETTINGS,
                      nRUN_SET = 0)
##----------------------------------------#

## 5. Compile model -----------------------
PCModelCompileModelWorkCase(dirSHELL = dirShell,
                            nameWORKCASE = nameWorkCase)
##----------------------------------------#


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


## 6. Run optimisation ---------------------------

### a. Define parameter values -------------------
# dataframe of parameter ranges (model parameters) and lags 

# Define the parameter values to be optimised
lower_bound <- c('mPLoadEpi' = 0.001, #minimum Ploading
                 'mPLoadEpi_lag' = 1, # need at least a year?
                 'fMarsh' = 0) # fraction of marsh area relative to lake
upper_bound <- c('mPLoadEpi' = 0.01, 
                 'mPLoadEpi_lag' = 5*1, # could be at least 5 year lag
                 'fMarsh' = 0.2) 

# error_pars <- c('mPLoadEpi' = 0.0005531256, #minimum Ploading
#                 'mPLoadEpi_lag' = 4.3935696417, # need at least a year?
#                 'fMarsh' = 0.8961858496)

### b. Define the desired future ------------------
# What is the objective
current_val <- c('mPLoadEpi' = 0.005) # the "unchanged" P-loading (before the measure is in place) - could also be a timeseries I guess?

# Define the desired future state(s)
desired_states <- data.frame(variable = c('oChlaEpi'),
                             target = c(20))


## Update the DATM file and recompile the model ---------------
# Report variables
lDATM_SETTINGS$auxils$iReport[which(rownames(lDATM_SETTINGS$auxils) %in% restart_states$state)] <- 0 # these can be turned off

## 3.4. Make and adjust cpp files ----------       
#  - nRUN_SET determines which forcings are switched on
PCModelAdjustCPPfiles(dirSHELL = dirShell,
                      nameWORKCASE = nameWorkCase,
                      lDATM = lDATM_SETTINGS,
                      nRUN_SET = 0)
##----------------------------------------#

## 5. Compile model -----------------------
PCModelCompileModelWorkCase(dirSHELL = dirShell,
                            nameWORKCASE = nameWorkCase)
##----------------------------------------#

# The copiled model needs to have all of the forcing variables in it (e.g. marsh area, P-loading)
# otherwise the .dll and DATM file won't match






### c. Define the objective function -----------------
# 1) run model
# 2) extract output
# 3) return a minimised value by comparing with desired future (DEOptim will make the value as negative as possible)
source(file.path(project_location, "scripts/optim_functions.R"))

obj_function <- function(val_pars, name_pars, future_states) {
  
  # For debugging ----------------- #
  # val_pars <- lower_bound
  # name_pars <- names(lower_bound)
  # future_states <- desired_states #
  #--------------------------------#
  
  model_output <- run_pathway(val_pars,
                              name_pars)
  
  eval_output <- evaluate_pathway(PCLake_output = model_output, 
                                  future_states = future_states,
                                  eval_target = function(x,y){(x-y)/y})#,
                                  # eval_days = 121:244)
  return(eval_output)
}



## Parallelization with FOREACH package
{
  nC <- parallelly::availableCores() 
  cl <- makeSOCKcluster(nC-2)
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(deSolve))
  clusterExport(cl, list("lDATM_SETTINGS", 'current_val', 'equilibrium_states',
                         "PCModelInitializeModel", "dirShell", "nameWorkCase", 'dirHome',
                         "PCmodelSingleRun", "RunModel", "run_pathway", "evaluate_pathway"))
  
  doSNOW::registerDoSNOW(cl)
  
  deoptim_control <- list(NP = 15 * length(lower_bound),
                          CR = 0.9,
                          F = 0.80, 
                          trace = TRUE,
                          itermax = 50, 
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

# iteration is a generation of a population
# one population is a collection of members
# one member is a single run of the objective function with specific parameter values

## Visualise the best member output
n_iter <- opt_pathway$optim$iter
iteration_summary <- as.data.frame(opt_pathway$member$bestmemit) |>
  mutate(fn_out = opt_pathway$member$bestvalit[1:n_iter])


ggplot(iteration_summary, aes(x=fMarsh, y= mPLoadEpi_lag, size = mPLoadEpi, colour = fn_out)) +
  geom_point() +
  scale_color_viridis_c() +
  coord_cartesian(xlim = c(0,0.2), 
                  ylim = c(1,5))

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
  slice_min(fn_out, prop = 0.25) |> # best (lowest) 25% of values
  ggplot(aes(x=fMarsh, y= mPLoadEpi_lag, size = mPLoadEpi, colour = fn_out)) + 
  geom_point(alpha = 0.5) + 
  scale_colour_gradient2() +
  coord_cartesian(xlim = c(0,0.2), 
                  ylim = c(1,5))
