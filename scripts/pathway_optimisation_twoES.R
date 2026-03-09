#--------------------------------------#
## Project: Pathway Optimisation Framework - multi-objective
## Script purpose: Simple example of an optimsation framework using PCLake for multiple objectives measures 
## Date: 2025-11-25 (mostly a copy from pathway_optimisation.R); updated 30/01/2026
## Author: Freya Olsson
#--------------------------------------#

library(tidyverse)
library(here)
library(DEoptim)
library(doSNOW)
library(parallelly)
library(ggh4x)

## 0. Settings --------------------------
##---------------------------------------#

## Global settings
options(scipen = 999) ## no scientific notation
save_output <- TRUE
make_plots <- TRUE
example_name <- 'multiES'

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

# Define the parameter values to be optimised (upper and lower) as well as
# the "unchanged" value (before the measure is in place) - could also be a timeseries I guess?
possible_measures <- read_csv(file.path(project_location, 'possible_measures.csv'), show_col_types = F) |> 
  filter(parameter %in% c('mPLoadEpi',
                          'mPLoadEpi_lag',
                          'fManVeg',
                          'fManVeg_lag'))

# see if there are other things required to run the measure optimisation
for (i in 1:length(possible_measures$parameter)) {
  if (!is.na(possible_measures$note[i])) {
    message(possible_measures$parameter[i], ' has requirements: ', possible_measures$note[i])
  } 
}

# check for associated parameters required for some optimisations - see possible_measures notes
if ('fManVeg' %in% possible_measures$parameter) {
  # fManVeg requires a cDayManVeg1
  lDATM_SETTINGS$params['cDayManVeg1', "sDefault0"] <- 259 # (259 = 19th Sept)
}

if(sum(str_detect(possible_measures$parameter, '_lag')) > 0){
  
  lag_var <- gsub('_lag', '', possible_measures$parameter[str_detect(possible_measures$parameter, '_lag')])
  main_var <- possible_measures$parameter[!str_detect(possible_measures$parameter, '_lag')]
  
  if (sum(!lag_var %in% main_var) > 0) {
    stop('You are missing the main var for a lagged variable!!')
  }
}

### b. Define the desired future ------------------
# What is the objective
# Define the desired future state(s)
desired_states_df <- data.frame(opt_var = c('oChlaEpi', 'aDSubVeg', 'aDFish'),
                                lower_range = c(0, 30, 6),
                                upper_range = c(55, 50, 8))

desired_states <- list(oChlaEpi = list(target = c(0,55),
                                       weights = 0.5),
                       aDSubVeg = list(target = c(30, 50),
                                       weights = 0.25),
                       aDFish = list(target = c(6, 8),
                                       weights = 0.25))



## Update the DATM file and recompile the model ---------------#
# Report variables
lDATM_SETTINGS$auxils$iReport[which(rownames(lDATM_SETTINGS$auxils) %in% restart_states$state)] <- 0 # these can be turned off
lDATM_SETTINGS$auxils$iReport[which(rownames(lDATM_SETTINGS$auxils) %in% names(desired_states))] <- 1 # report optim vars
lDATM_SETTINGS$params$iReport[which(rownames(lDATM_SETTINGS$params) %in% possible_measures$parameter)]  # report measure params
lDATM_SETTINGS$auxils[which(rownames(lDATM_SETTINGS$auxils) %in% 'uPLoadEpi'), ] <- 1 # also report the auxillary variable for the PLoadEpi

# forcing variables --------------#
# anything that is being lagged needs to be in the forcings before compilation
for (i in possible_measures$parameter) {
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

#' Define the objective function with the output to be optimised
#'
#' @param val_pars a vector of parameter values to be evaluated
#' @param name_pars a vector of the parameter names, same order as val_pars
#' @param future_states the desired future values to be compared
#' @description the function is split into two parts - run the pathway returns PCLake output, this is then compared with the desired states
#' @returns an single evaluation value to be optimised
#' @export
#'
#' @examples obj_function(val_pars = val_pars, name_pars = name_pars, future_states = desired_states)

obj_function <- function(val_pars, name_pars, future_states) {
  
  # For debugging ----------------- #
  # val_pars <- possible_measures$upper_bound
  # name_pars <- possible_measures$parameter
  # future_states <- desired_states #
  #--------------------------------#
  
  model_output <- run_pathway(val_pars, name_pars, current_val = possible_measures$current_val,
                              initial_conditions = equilibrium_states)
  
  eval_output <- evaluate_pathway(PCLake_output = model_output, 
                                  future_states = future_states,
                                  eval_days = 50:300, 
                                  eval_funs = mean,
                                  eval_target = list("oChlaEpi" = range_obj,
                                                     "aDSubVeg" = range_obj,
                                                     "aDFish" = range_obj) # see optim_functions.R
                                  
  )
  return(eval_output)
}


## 6.B Run optimisation ---------------

## Parallelisation with FOREACH package
{
  nC <- parallelly::availableCores() 
  cl <- makeSOCKcluster(nC-2)
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(deSolve))
  clusterExport(cl, list("lDATM_SETTINGS", 'possible_measures', 'equilibrium_states',
                         "PCModelInitializeModel", 
                         "above_obj", "below_obj", "exact_obj", "range_obj",
                         "dirShell", "nameWorkCase", 'dirHome',
                         "PCmodelSingleRun", "RunModel", 'run_pathway', 'evaluate_pathway'))
  
  doSNOW::registerDoSNOW(cl)
  
  deoptim_control <- list(NP = 10 * nrow(possible_measures), #number of population members, should be at least 10 times the length of the parameter
                          CR = 0.9, # crossover probability between 0-1, default in 0.5
                          F = 0.80, # differential weighting factor between 0-2. Default to 0.8
                          itermax = 100, # the maximum iteration (population generation) allowed
                          reltol = 0.01, # relative convergence tolerance. The algorithm stops if it is unable to reduce the value by a factor of reltol * (abs(val) + reltol) after steptol steps.
                          steptol = 20, # number of minimum ssteps
                          # c = 0.05, strategy = 6, p = 0.10
                          c = 0.05,  # the speed of crossover adaptation. Between 0-1. Higher c give more weight to the current successful mutations
                          strategy = 6, # DE / current-to-p-best / 1
                          p = 0.1, # the top (100 * p)% best solutions are used in the mutation
                          storepopfrom = 0, # 0 = store all intermediate pops
                          trace = TRUE,
                          parallelType = 'foreach')
  
  # Several conditions can cause the optimization process to stop:
  # if the maximum number of iterations is reached (itermax), or
  # if a number (steptol) of consecutive iterations are unable to reduce the best function value by a certain amount (reltol * (abs(val) + reltol)). 
  # 100*reltol is approximately the percent change of the objective value required to consider the parameter set an improvement over the current best member.
  
  set.seed(1234)
  
  opt_pathway <- DEoptim::DEoptim(lower = possible_measures$lower_bound,
                                  upper = possible_measures$upper_bound, 
                                  fn = obj_function,
                                  control = deoptim_control, 
                                  name_pars = possible_measures$parameter,
                                  future_states = desired_states)
  
  parallel::stopCluster(cl)
}

# The output of DEoptim is based on members, iterations, and populations
# iteration is a generation of a population
# one population is a collection of members
# one member is a single run of the objective function with specific parameter values
# rename all the missing names
names(opt_pathway$member$lower) <- possible_measures$parameter
names(opt_pathway$member$upper) <- possible_measures$parameter
colnames(opt_pathway$member$bestmemit) <- possible_measures$parameter
names(opt_pathway$optim$bestmem) <- possible_measures$parameter

## Visualise the best member output from each iteration
n_iter <- opt_pathway$optim$iter
iteration_summary <- as.data.frame(opt_pathway$member$bestmemit) |>
  mutate(fn_out = opt_pathway$member$bestvalit[1:n_iter])

if (make_plots) {
  p1 <- ggplot(iteration_summary, aes(x=fManVeg_lag, y= fManVeg, size = mPLoadEpi, colour = fn_out)) +
    geom_point() +
    scale_color_viridis_c()  +
    theme_bw()
  if (save_output) {
    ggsave(plot = p1, filename = file.path(project_location, 'output', 'plots', paste0('bestmemit_', example_name, '.png')), 
           width = 10, height = 8, unit = 'cm')
  }
}


# write output for later?------------------------------
if (save_output) {
  write_rds(c(list(possible_measures = possible_measures, desired_states = desired_states), opt_pathway$optim),
            file = file.path(project_location, 'output', paste0('summary_',example_name, '.RData'))) # summary of the optimisation example
  write_csv(iteration_summary, file = file.path(project_location, 'output',  paste0('bestmemit_',example_name, '.csv'))) # best pathway for each population 
}

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


if (make_plots) {
  p2 <- last_iteration |> 
    # filter(fn_out == 0) |> 
    ggplot(aes(x=fManVeg_lag, y= fManVeg, size = mPLoadEpi_lag, colour = fn_out)) + 
    geom_point() + 
    scale_colour_viridis_c() +
    theme_bw()
  if (save_output) {
    ggsave(plot = p2, filename = file.path(project_location, 'output', 'plots', paste0('lastpop_', example_name, '.png')), 
           width = 12, height = 9, unit = 'cm')
  }
}

### 7.b Extract the state values -----------
# what are the values of the optimised variables?
state_opt <- foreach(i = 1:nrow(last_iteration),
                     .combine = bind_rows, 
                     .multicombine = TRUE,
                     .packages=c("tidyverse", "deSolve"),
                     .options.snow = opts) %dopar% {
                       
                       val_pars <- last_iteration[i,possible_measures$parameter] |> unlist()
                       cur_val <- possible_measures$current_val
                       name_pars <- last_iteration[i,possible_measures$parameter] |> names()
                       
                       df_pars <- data.frame(variable = name_pars, output = val_pars)
                       
                       run_pathway(val_pars, name_pars, cur_val) |>
                         mutate(year = floor((time-1)/365) + 1,
                                doy = yday(as_date(time - (year * 365) + 364, origin = '2025-01-01'))) |> 
                         filter(year == max(year), # filters to summer in the last year of the simulation
                                doy %in% 50:300) |> 
                         select(names(desired_states)) |> 
                         summarise(across(any_of(names(desired_states)), mean)) |> 
                         pivot_longer(cols = any_of(names(desired_states)),
                                      names_to = 'variable',
                                      values_to = 'output') |> 
                         bind_rows(df_pars) |> 
                         mutate(ID = i)
                       
                     }


if (make_plots) {
    p3 <- state_opt |> 
    pivot_wider(id_cols = ID, names_from = variable, values_from = output) |>
    pivot_longer(cols = names(desired_states), names_to = 'opt_var', values_to = 'out') |> 
    # full_join(desired_states_df, by = join_by(opt_var)) |> 
    # filter(between(out, lower_range, upper_range)) |> # check within range
    # filter(n() == nrow(desired_states_df), .by = ID) |> # only where both states pass
    ggplot(aes(x=mPLoadEpi, y = fManVeg_lag, size = fManVeg, colour = mPLoadEpi_lag)) + geom_point() +
    facet_wrap(~opt_var, scales = 'free') +
    scale_colour_viridis_c(option = 'A', begin = 0.3, end = 0.9) +
    theme_bw()
  
  if (save_output) {
    ggsave(plot = p3, filename = file.path(project_location, 'output', 'plots', paste0('lastpopstate_', example_name, '.png')), 
           width = 12, height = 8, unit = 'cm')
  }
}



if (save_output) {
  write_csv(last_iteration, file = file.path(project_location, 'output',  paste0('lastpop_',example_name,'.csv')))  # each member values of the final population
  
  state_opt |> 
    pivot_wider(id_cols = ID, names_from = variable, values_from = output) |> 
    pivot_longer(cols = names(desired_states), names_to = 'opt_var', values_to = 'out') |> 
    write_delim(file = file.path(project_location, 'output',  paste0('lastpopstate_',example_name, '.csv')))
}

### 7.c Extract the time series ---------
# above we extracted the last year (the year being optimised) but in this case we want to 
# look at the full pathway to acheive the final target

# what is the "pathway" to the optimised value?
state_pathways <- foreach(i = 1:nrow(last_iteration),
                          .combine = bind_rows, 
                          .multicombine = TRUE,
                          .packages=c("tidyverse", "deSolve"),
                          .options.snow = opts) %dopar% {
                            
                            val_pars <- last_iteration[i,possible_measures$parameter] |> unlist()
                            cur_val <- possible_measures$current_val
                            name_pars <- last_iteration[i,possible_measures$parameter] |> names()
                            
                            df_pars <- data.frame(variable = name_pars, output = val_pars)
                            
                            run_pathway(val_pars, name_pars, cur_val) |>
                              mutate(year = floor((time-1)/365) + 1,
                                     doy = yday(as_date(time - (year * 365) + 364, origin = '2025-01-01'))) |> 
                              filter(# filters to summer all years of the simulation
                                doy %in% 50:300) |> 
                              select(c('year', names(desired_states))) |> 
                              group_by(year) |> 
                              summarise(across(any_of(names(desired_states)), mean)) |> 
                              bind_cols(pivot_wider(df_pars, names_from = variable, values_from = output)) |>  
                              mutate(ID = i)
                            
                          }



if (make_plots) {
  p4 <- state_pathways |> 
    filter(year == max(state_pathways$year),
           oChlaEpi <= 20,
           between(aDSubVeg, 50, 100),
           between(aDFish, 6, 8)) |> 
    select(ID) |> 
    left_join(state_pathways, by = join_by(ID)) |> 
    mutate(.by = year, ID = row_number()) |> # renumber the pathways 
    mutate(fManVeg_use = ifelse(fManVeg_lag < year, fManVeg, 0),
           mPLoadEpi_use = ifelse(mPLoadEpi_lag < year, mPLoadEpi, possible_measures$current_val[which(possible_measures$parameter == 'mPLoadEpi')]))  |> 
    select(all_of(c('ID', 'year', 'aDSubVeg', 'oChlaEpi', 'mPLoadEpi_use', 'fManVeg_use'))) |> 
    pivot_longer(cols = !any_of(c('ID','year'))) |> 
    # filter(ID %in% 7:9) |>
    ggplot(aes(x=year, y = value, 
               # size = value, 
               colour = name)) +
    geom_line(lineend = 'round', linejoin = 'round', linemitre = 1, linewidth = 1) +
    ggh4x::facet_nested_wrap(vars(ID, name), scales = 'free_y',  nrow = 20, 
                             strip.position = 'right', dir = 'v', remove_labels = 'y',
                             nest_line = element_line(colour = 'black'), 
                             strip = strip_nested(text_y = list(element_text(), 
                                                                element_text(colour = 'white', 
                                                                             size = 1)),
                                                  background_y = list(element_rect(),
                                                                      element_blank()), 
                                                  by_layer_y = TRUE)) +
    theme_bw(base_size = 12) +
    theme(panel.border = element_rect(colour = 'black'),
          legend.position = 'top',
          panel.spacing.y = unit(c( rep( c( rep(0.2,4), 0.6), 3), rep(0.2,4)),
                                 "lines")
    ) +
    facetted_pos_scales(y = list(name == "fManVeg_use" ~ scale_y_continuous(limits = c(0,1),
                                                                           n.breaks = 2),
                                 name == "mPLoadEpi_use" ~ scale_y_continuous(limits = c(0,0.01),
                                                                              n.breaks = 2),
                                 name == "oChlaEpi" ~ scale_y_continuous(limits = c(0,200),
                                                                         n.breaks = 2),
                                 name == "aDSubVeg" ~ scale_y_continuous(limits = c(0,100),
                                                                         n.breaks = 2)
                                 
    )
    ) +
    scale_colour_manual(values = c('black', 'orange','seagreen', 'orchid', 'grey'))
  if (save_output) {
    ggsave(plot = p4, filename = file.path(project_location, 'output', 'plots', paste0('lastpoppathways_', example_name, '.png')), 
           height = 20, width = 30, units = 'cm')
  }
}



if (save_output) {
  write_csv(state_pathways, file = file.path(project_location, 'output',  paste0('lastpoppathways_',example_name,'.csv')))  # the annual output of the final population
}
