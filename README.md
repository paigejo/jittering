# jittering
Accounting for jittered DHS observation locations

## Repository Structure
+ The ``code'' directory stores all R and C++ code, as well as compiled C++ code for Template Model Builder (TMB)
    + The `setup.R' script should be the first script run
    + The `exampleSPDEJitter.R' script shows how to run models that account for jittering as well as ones that do not
    + The `makeIntegrationPoints.R' script contains functions for generating integration points and weights for TMB
    + `modSPDEJitter.cpp' contains the TMB code for accounting for jittering given a number of integration points
+ The `Figures' directory stores all figures produced
