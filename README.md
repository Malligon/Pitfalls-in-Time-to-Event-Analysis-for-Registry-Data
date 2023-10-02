## Pitfalls-in-Time-to-Event-Analysis-for-Registry-Data

<!-- badges: start -->

Appendix of the study "Pitfalls in Time to Event Analysis for Registry Data: A Simulation &amp; Real Case Approaches based Tutorial."

Contains 3 main files:
  - `Survival_function.R`
  - `CEREDIH_randomised.csv`
  - `Simulation_and_Real_life_data_Survival_methods.R`

# Survival_function.R
Contains the following functions:
  - `plot_km_trunc`: from a survfit, or a Surv formula, draw the Kaplan-Meier with confidence intervals, including a display of the risk table. Note that unlike functions from other packages, this function accurately take left troncation (and any kind of late entry in the risk set) into account for the risk table.
  - `plot_km_trunc_CR`: Same than the previous function, but draws the cumulative incidence function (CIF) instead to fit competing risks models. Left truncation, or any kind of late entry, are accurately taken account in the risk table.
  - `rec_even`t: Calculate the mean cumulative number of events in recurrent events situation, alongside confidence intervals. Note that to our knowledge, this is the only available tool to display the confidence interval of the mean cumulative number of events using only the survival package.
  - `plot.rec_event`: plot any object created from the rec_event function, including confidence interval and risk table.

# CEREDIH_randomised.csv
A randomised version of the CEREDIH registry, including only variables used in the study. Please, understant we cannot provide the real full dataset due to the sensitivy of the informations. this base should be used only as a tool to run the codes dedicated to real life data in `Simulation_and_Real_life_data_Survival_methods.R`.

# Simulation_and_Real_life_data_Survival_methods.R
Script showing how to use the different functions from `Survival_function.R`, and how to display survival curves. All simulations are also available allowing readers to reproduice the graphics for the simulation part. Please note that using the CEREDIH_randomised.csv dataset won't reproduice the survival curves displayed in the study in the CEREDIH part.

# How to use

First run `Survival_function.R`. Then, import these libraries.

```r
library(survival)
library(stringr)
library(ggsci)
```

Each section, for the simulation part, is divided in two parts. The simulation of patients, and then the survival analysis. Let's jump to the section `Effect of truncature`. First, we need to simulate the data.

```r
set.seed(37)

N = 1000

shape = 2
scale = 30

T <- rweibull(N, shape, scale) # true times
C <- runif(N, 0, 85) # censored times

Ttrunc <- runif(N, 0, 50) # Time of truncature

status <- T <= C
Tobs <- pmin(T, C)

Tobs2 <- Tobs[Tobs >= Ttrunc] # only times of event after truncature could be observed
status2 <- status[Tobs >= Ttrunc] # idem for the status of event
Ttrunc2 <- Ttrunc[Tobs >= Ttrunc] # idem for times of truncature

DB <- data.frame(Tobs = c(Tobs2, Tobs2), 
                 status = c(status2, status2),
                 start = c(Ttrunc2, rep(0, length(status2))),
                 trunc = c(rep(1, length(status2)), rep(2, length(status2))))
```

Here `DB` contains each simulated patient twice. The variable `trunc` will allow us to compare the naive analysis without truncature (`trunc = 2`) and the right one with truncature taken into account (`trunc = 1`).

Now we just have to run the function `plot_km_trunc` with the right parameters to draw the curves.

``` r
plot_km_trunc(survfit(Surv(start, Tobs, status) ~ trunc, data = DB),
              DB,
              Trunc = TRUE,
              col = c("blue", "red"),
              strata.names = c("KM", "Naive"),
              xlim = c(-1, 50),
              ylim = c(-0.05, 1.05),
              legend = FALSE,
              lwd = 0.5,
              legend.lwd = 0.5)
lines(xaxis, exp(-(xaxis/scale)^shape), type = "l",
      lwd = 0.5 ) # true survival
legend("topright", legend = c("True survival", "KM taking into account truncation", "KM ignoring truncation"), col = c("black", "blue", "red"), lty = 1)
```

![Comparison truncature](https://github.com/Malligon/Pitfalls-in-Time-to-Event-Analysis-for-Registry-Data/assets/43923608/08ec59a6-f4da-4254-af8a-9c58c0f99f7b)

In a real study, the database should contain only the patients with the group `trunc = 1` and the codes would look like that.

``` r
DB <- data.frame(Tobs = Tobs2, 
                 status = status2,
                 start = Ttrunc2)

plot_km_trunc(survfit(Surv(start, Tobs, status) ~ 1, data = DB),
              DB,
              Trunc = TRUE,
              xlim = c(-1, 50),
              ylim = c(-0.05, 1.05),
              legend = FALSE,
              lwd = 0.5,
              legend.lwd = 0.5)
```

If you want to use codes for real life data, download and import in R the `CEREDIH_randomised.csv database`. Just modify the `path` object in R with a string corresponding of the location of the database. It has to end with a `\`.

``` r
path <- "path of your file" # example "D:/Desktop/"

DB <- read.csv2(paste(path, "CEREDIH_randomised.csv", sep = ""), dec=",", na.strings=c("NA",""," ",".","--", "N/A"))
```

You can now use `Simulation_and_Real_life_data_Survival_methods.R`.
