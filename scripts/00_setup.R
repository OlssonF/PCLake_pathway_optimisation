#--------------------------------------#
## Project: Pathway Optimisation Framework
## Script purpose: set up the directories, folders, and files 
## Date: 2026-07-30
## Author: Freya Olsson
# Created with R version 4.5.2 (2025-10-31 ucrt)
#--------------------------------------#

#------------------------------------------#
# ----- 1. Install required packages -------
#------------------------------------------#

required_packages <- c(
  'tidyverse', # 2.0.0
  'here', # 1.0.2
  'DEoptim', # 2.2.8
  'doSNOW', # 1.0.20
  'parallelly', # 1.46.1
  'ggh4x', # 0.3.1
  'ggpubr', # 0.6.2
  'GGally' # 2.4.0
)

install.packages(required_packages)

#

#------------------------------------------#
# ----- 2. Set up the workcase -------------
#------------------------------------------#
library(tidyverse)
# as per the R scripts provided in the PCLake repository, set up the workcase with
# initial copies of the source_cpp, model_code, source_cpp_adjusted folders are made
project_location <- here::here()
DATM_file <- "PL613162PLUS_pathway_optim.xls"

dirHome <- str_split(project_location,  "(?=PCModel1350)", simplify = T)[1,1]	# location of the PCModel1350 folder
dirShell <- str_split(project_location,  "(?<=PCShell)", simplify = T)[1,1]	#  PCShell folder path
dirCpp_root <- list.dirs(dirHome)[which(str_detect(list.dirs(dirHome),"3.01/PCLake_plus"))] # location of C++ code
nameWorkCase <- tail(str_split_1(project_location, "/"), n = 1) # workcase name

# Load PCLake functions from the main repository
source(file.path(dirShell, "scripts", "R_system", "functions.R"))
source(file.path(dirShell, "scripts", "R_system", "functions_PCLake.R")) 

# This will copy the relevant source cpp files across from the root dir
PCModelWorkCaseSetup(dirSHELL = dirShell, 
                     dirCPP_ROOT = dirCpp_root,
                     nameWORKCASE = nameWorkCase)

#

# -----------------------------------------#
# ----- 3. Copy the DATM file across -------
# -----------------------------------------#
copy_datm <- file.copy(file.path(project_location, DATM_file),
                       file.path(str_remove(dirShell, "/[^/]+$"), DATM_file),
                       overwrite = F) # does not overwrite an existing DATM with the same name
if (copy_datm) {
  message("DATM copied to right folder")
} else {
  message("A DATM file with that name already exists")
}



