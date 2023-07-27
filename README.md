## Pitfalls-in-Time-to-Event-Analysis-for-Registry-Data
Appendix of the study "Pitfalls in Time to Event Analysis for Registry Data: A Simulation &amp; Real Case Approaches based Tutorial."

Contains 3 main files:
  - Survival_function.R
  - CEREDIH_randomised.csv
  - Simulation_and_Real_life_data_Survival_methods.R

# Survival_function.R
Contains the following functions:
  - plot_km_trunc: from a survfit, or a Surv formula, draw the Kaplan-Meier with confidence intervals, including a display of the risk table. Note that unlike functions from other packages, this function accurately take left troncation (and any kind of late entry in the risk set) into account for the risk table.
  - plot_km_trunc_CR: Same than the previous function, but draws the cumulative incidence function (CIF) instead to fit competing risks models. Left truncation, or any kind of late entry, are accurately taken account in the risk table.
  - rec_event: Calculate the mean cumulative number of events in recurrent events situation, alongside confidence intervals. Note that to our knowledge, this is the only available tool to display the confidence interval of the mean cumulative number of events using only the survival package.
  - plot.rec_event: plot any object created from the rec_event function, including confidence interval and risk table.

# CEREDIH_randomised.csv
A randomised version of the CEREDIH registry, including only variables used in the study. Please, understant we cannot provide the real full dataset due to the sensitivy of the informations. this base should be used only as a tool to run the codes dedicated to real life data in Simulation_and_Real_life_data_Survival_methods.R.

# Simulation_and_Real_life_data_Survival_methods.R
Script showing how to use the different functions from Survival_function.R, and how to display survival curves. All simulations are also available allowing readers to reproduice the graphics for the simulation part. Please note that using the CEREDIH_randomised.csv dataset won't reproduice the survival curves displayed in the study in the CEREDIH part.

# How to use

First run "Survival_function.R".

If you want to use codes for real life data, download and import in R the CEREDIH_randomised.csv database.

You can now use "Simulation_and_Real_life_data_Survival_methods.R"
