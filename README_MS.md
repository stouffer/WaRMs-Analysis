# MS_VariableResponsesAlpinePlants
Files pertaining to manuscript: Variable responses of alpine-plant communities to warming and loss of dominant species

Code uses R version 3.6.2, brms package version 2.11.1, and tidyverse package version 1.3.0 

Files in this repository:

example.Rdata
Contains data from Canada low elevation as an example of data structure and distribution. This Rdata file contains data prepped for poisson and beta distribution stan models, meaning species names and other factors have been indexed.

brmsworkflow.R 
Contains an example workflow for a site. Runs the models and comparisons for both poisson and beta error distributions. Poisson error distribution models are what is presented in the manuscript. In general, the beta distribution did not converge well for our data so were considered invalid, however, they may be useful for future present cover or proportional cover data. 

brmsModelFunctions.R
Contains poisson error distribution brms model code for those models described in the manuscript.

beta_distribution/brmsBetaFunctions.R
Contains beta error distribution brms model code. These models were not presented in the manuscript as they did not converge well with our current data. 

beta_distribution/brmsBetaworkflow.R
Contains beta error distribution workflow, including how to run models and compare them. These models were not presented in the manuscript as they did not converge well with our current data. 
