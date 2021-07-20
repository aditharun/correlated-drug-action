#metascript

setwd("scripts")

source("doses.hsa.R")

source("doses.loewe.R")
source("sp.doses.hsa.R")
source("sp.doses.loewe.R")

print("DONE with model estimates")

setwd("../figures")

#source("render_figures_doses.R")

#source("generate_dose_supp_data.R")

#print("DONE with figure and table rendering")