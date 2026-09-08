#--------------------------------------#
## Project: Pathway optimisation framework - get initialisation 
## Script purpose: runs the initialisation step for the the pathway optimisation and plots the 
##                 30 year simulation for SI
## Date: 2026-03-24
## Author: Freya Olssn
# Created with R version 4.5.2 (2025-10-31 ucrt)
#--------------------------------------#

library(tidyverse)
library(here)

## 0. Settings --------------------------
##---------------------------------------#

## Global settings
options(scipen = 999) ## no scientific notation
save_output <- TRUE
make_plots <- TRUE
example_name <- 'simple'
## 1. Directory settings ---------------------------------------------------------
## using relative paths in which the project and script is saved in the work_cases
## "scripts" contains only the PCLake functions

project_location <- here()
DATM_file <- "PL613162PLUS_pathway_optim.xls"

dirHome <- str_split(project_location,  "(?=PCModel1350)", simplify = T)[1,1]	# location of the PCModel1350 folder
dirShell <- str_split(project_location,  "(?<=PCShell)", simplify = T)[1,1]	#  PCShell folder path
dirCpp_root <- list.dirs(dirHome)[which(str_detect(list.dirs(dirHome),"3.01/PCLake_plus"))] # location of C++ code
nameWorkCase <- tail(str_split_1(project_location, "/"), n = 1) # workcase name
fileDATM <- file.path(str_remove(dirShell, "/[^/]+$"), DATM_file)
dirSave <- dirShell
# ----------------------------------------------------------------------------- #

## load external functions from the scripts folder
source(file.path(dirShell, "scripts", "R_system", "functions.r"))
source(file.path(dirShell, "scripts", "R_system", "functions_PCLake.r")) 


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
                                                 locDATM = "excel",
                                                 readAllForcings = F)
##----------------------------------------#

## Modifications can be made to the DATM file here! (those that need to be compiled)
## Might be a good idea to make sure the default lake parameters are loaded

# Report restart variables
restart_states <- read_table(file.path(project_location,'input', 'restart_states.txt'), col_names = 'state', show_col_types = F)
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


vars_plot <-c('Epilimnion chla' = 'oChlaEpi',
              'Epilimnion total phosphorus' = 'oPTotWEpi',
              'Secchi depth' = 'aSecchiT')
names(vars_plot) <- NULL

library(ggplot2)

df_labels <- data.frame(variable = vars_plot ) |> 
  mutate(facet_label = c("Epilimnion~chlorophyll~a~concentration ~ (mu*g ~ L^{-1})",
                         'Epilimnion~total~phosphorus~concentration~(g~m^{-3})',
                         'Secchi~depth~(m)'))

figureS1 <- ggpubr::ggarrange(PCModel_run_baseline |> 
                    select(all_of(c('time', vars_plot))) |> 
                    pivot_longer(cols = -time, 
                                 names_to = 'variable',
                                 values_to = 'value') |>
                    full_join(df_labels) |> 
                    ggplot(aes(x=time/365, y = value)) + 
                    geom_line() +
                    facet_wrap(~facet_label, nrow = 3, scales = 'free_y',
                               labeller = label_parsed) +
                    theme_bw(base_size = 14) +
                    scale_x_continuous(name = 'time (year)'),
                  
                  
                  PCModel_run_baseline |> 
                    slice_tail(n = 365) |> 
                    select(all_of(c('time', vars_plot))) |> 
                    pivot_longer(cols = -time, 
                                 names_to = 'variable',
                                 values_to = 'value') |>
                    full_join(df_labels) |> 
                    ggplot(aes(x=yday(as_date(time)), y = value)) + 
                    geom_line() +
                    facet_wrap(~facet_label, nrow = 3, scales = 'free_y',
                               labeller = label_parsed) +
                    theme_bw(base_size = 14) +
                    scale_x_continuous(name = 'day of year'),
                  
                  ncol = 2, widths = c(1, 0.7))
ggsave(figureS1, filename = file.path(here::here(), 'output/plots/ms/FigureS1.png'),
         height = 16, width = 27, units = 'cm')

