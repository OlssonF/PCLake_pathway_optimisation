# PCLake_pathway_optimisation

This project aims to develop a framework for pathway optimisation to generate target-seeking scenarios for management. The framework combines the process-model PCLake and a differential evolution algorithm to identify combinations of measures that could be implemented for management and restoration that improve the state of the system based on sinlge or multi-objective targets.

## PCLake+ set-up:

This implementation of PCLake uses the DATM file saved in the main PCLake directory (see the set up in the main [PCLake repository](https://github.com/pcmodel/PCModel/) to show where this file sits). The file name should be ./PCModel/Licence_agreement/I_accept/PCModel1350/PCModel/3.00/Models/PCLake+/6.13.16/PL613162PLUS_pathway_optim.xls. To run PCLake in R uses the scripts from ./PCModel/Licence_agreement/I_accept/PCModel1350/PCModel/3.00/Models/PCLake+/6.13.16/PCShell/scripts.

## Notes on doing the pathway optimisation:

-   The measure are constrained by their *magnitude* (e.g. P load reduction or proportion of Marsh area) and their *timing* (using the \_lag suffix). Once a measure is in place (t \> lag) then it remains in place for the duration of the simulation.

-   At the moment this framework is limited to the optimisation of model states (e.g. `aChlaEpi` or `aDFish`) and the model parameters or drivers.

## Setting up the project:

-   This directory only includes the optimisation project files. To actually run the workflow you need a lot more files and set up

-   Start with cloning the [PCLake repository](https://github.com/pcmodel/PCModel/) or open a docker container that contains a fixed version of PCLake

-   Then clone this project repository into the *workcases* subdirectory for PCShell (R implementation) (./PCModel/Licence_agreement/I_accept/PCModel1350/PCModel/3.00/Models/PCLake+/6.13.16/PCShell/work_cases)

-   Then open the R project (PCLake_pathway_optimisation.Rproj)

## Running the optimisation examples

1.  scripts/00_setup.R organises the workcase folder - you should only do this if you have cloned the PCLake repo first and then cloned this directory into the. This script will:

-   copy the DATM file to the right location (./PCModel/Licence_agreement/I_accept/PCModel1350/PCModel/3.00/Models/PCLake+/6.13.16/) and
-   set up all the cpp, model_code folders etc.
-   You only need to do this once!

3.  scripts/01a_pathway_optimisation_simple.R and scripts/01b_pathway_optimisation_multi.R are the two scripts for generating the example cases from Olsson et al., 2026 (submitted to JEM).
4.  scripts/02_analysis.R is for generating the plots and values found in the manuscript.
5.  scripts/example_model_initialisation.R is a script that runs the initialisation step (long PCLake run) to generate "initial conditions" for optimisation run. Figure is for the SI.
6.  R/optim_functions.R contains some custom functions that help to run the model and evaluate output.
7.  scripts/archive/ contains scripts that were not used in the paper but other implementations of the framework (e.g. robust, prioritisation, two horizons etc). DO NOT USE AS IS.
