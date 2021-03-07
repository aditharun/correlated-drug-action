setwd("scripts")

source("clinicaltrial.R")

print("DONE WITH CLINICAL TRIALS")

source("placebo.R")

print("DONE WITH PLACEBOS")

source("summary_clintrial.R")

print("DONE WITH SUPP TABLE and REF TABLE")

setwd("../figures")

source("render_figures_clintrial.R")

print("DONE WITH FIGURES for CLINICAL TRIALS")

#merge supp file 1 and 3 into an excel sheet with two different tabs MANUALLY 