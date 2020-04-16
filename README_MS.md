# MS_VariableResponsesAlpinePlants
Files pertaining to manuscript: Variable responses of alpine-plant communities to warming and loss of dominant species

brmsworkflow.R 
Contains an example workflow for a site. Runs the models and comparisons for both poisson and beta error distributions. Poisson error distribution models are what is presented in the manuscript. In general, the beta distribution did not converge well for our data so were considered invalid, however, they may be useful for future present cover or proportional cover data. 

brmsModelFunctions.R
Contains poisson error distribution brms model code for those models described in the manuscript.

brmsBetaFunctions.R
Contains beta error distribution brms model code. These models were not presented in the manuscript as they did not converge well with our current data. 

