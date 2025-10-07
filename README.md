# PCLake_pathway_optimisation

This project aims to develop a framework for pathway optimisation based on the NFF. The framework uses the PCLake process-model and a differential evolution algorithm to identify combinations of measures that could be implemented for management and restoration that improve the state of the system based on sinlge or multi-objective targets.

## PCLake set-up:

This implementation of PCLake uses the DATM file saved in the main PCLake directory (see the set up in the main [PCLake repository](https://github.com/pcmodel/PCModel/) to show where this file sits). The file name should be ./PCModel/Licence_agreement/I_accept/PCModel1350/PCModel/3.00/Models/PCLake+/6.13.16/PL613162PLUS_pathway_optim.xls. To run PCLake in R uses the scripts from ./PCModel/Licence_agreement/I_accept/PCModel1350/PCModel/3.00/Models/PCLake+/6.13.16/PCShell/scripts

## Notes on doing the pathway optimisation:

-   The measure are constrained by their *magnitude* (e.g. max P load or proportion of Marsh area) and their *timing* (using the \_lag suffix). Once a measure is in place (t \> \_lag) then it remains in place for the duration of the simulation.

-   At the moment this framework is limited to the optimisation of model states (e.g. aChlaEpi or aDFish) and the model parameters or drivers.

## Setting up the project:

-   This directory only includes the optimisation project files. To actually run the workflow you need a lot more files and set up

-   Start with cloning the PCLake repository or open a docker container that contains a fixed version of PCLake

-   Then clone this project repository into the *workcases* subdirectory for PCShell (R implementation) (./PCModel/Licence_agreement/I_accept/PCModel1350/PCModel/3.00/Models/PCLake+/6.13.16/PCShell/work_cases)

-   Then open the R project (PCLake_pathway_optimisation.Rproj)
