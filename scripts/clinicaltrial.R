#CODE that I used to get results for my data
#Can be rewritten to be extended such that the functions here compute the reuslts for data that a person wants to run

library(tidyverse)
library(readxl)

source("randomization.R")
source("estimation_methods.R")
source("functions.R")

set.seed(1498)

# Note: store data for each trial in separate folder, folder has 3 csv's with time and pfs as two columns. All the folders for each trial must be stored in one bigger super folder

#---CONFIGURABLE PARAMETERS--------#
indir <- getwd()

#results are stored in this directory
outdir <- file.path(indir, "..", "results.clinical")

#path to raw data
path.to.data <- file.path(indir, "..", "raw-data", "clinical_trials")

#parameters

#how many samples for generating the null distribution
nullsamples <- 100000

#for linearly interpolating PFS curve data, what interval (in months) should we use
interval <- 0.01

#censor the last .03 (right tail) PFS data because its usually very noisy
censor <- 0.97

#number of samples in bootstrap esitmate of p-value (limit of resolution is 1/n) where n is the value for n.boot
n.boot <- 5000


#Set to FALSE if we do not want to run simulation estimates (default), set to true otherwise
run_sim <- FALSE

main.AssembleData <- function(path.to.data, interval, nullsamples, censor, dir.out){

  print("starting read in of files")

  files <- list.files(path = path.to.data, full.names = TRUE, pattern="*.csv",recursive = TRUE)
  files.w.placebo <- grep("PL", files)
  files <- files[-files.w.placebo]

  data <- lapply(files, get.data)
  names.combinations <- list.files(path=path.to.data, pattern="*.csv", recursive=TRUE)[-files.w.placebo]
  names.combinations <- gsub(".csv", "", names.combinations)
  names(data) <- names.combinations

  path.to.patient.numbers <- file.path(path.to.data, "patient_numbers.xlsx")

  num.trials <- length(data)/3

  for (x in 1:num.trials){
    subset <- data[c(3*x-2, 3*x-1, 3*x)]
    assemble.data(subset, path.to.patient.numbers, interval, nullsamples, censor, outdir)
  }
}

main.compute <- function(interval, censor, dir.out, run_sim){

  trials <- list.dirs(dir.out, recursive=FALSE, full.names=TRUE)
  num.trials <- length(trials)

  track.results <- file.path(dir.out, "results.csv")

  print("starting estimation")

  num.trials <- length(track.results)
  possible_rho <- seq(-1, 1, 0.01)

  counter <- 1
  for (combo in trials){

    print(combo)

    df <- compute(name=combo, interval=interval, censor=censor, possible_rho=possible_rho, run_sim=run_sim)

    plotter(name=combo, run_sim=run_sim)

    if (!file.exists(track.results)){
        write.table(df, file=track.results, append=TRUE, sep=",", row.names=FALSE, col.names=TRUE)
    } else{
        write.table(df, file=track.results, append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
    }

    print((counter/num.trials)*100)
    counter <- counter + 1

  }
}

#make sure to call this function in the directory where the results are stored aka where results.clinical is!
main.ci <- function(interval, ci, n.boot, censor, path.to.results){
  print("starting confidence intervals")
  dir <- path.to.results
  trials <- list.dirs(dir, recursive=FALSE)
  num.trials <- length(trials)

  track.results <- file.path(dir, "ci.csv")

  for (idx in 1:num.trials){

    print("starting new ci")
    combo <- trials[idx]

    conf.int <- confidence.interval(combo, n.boot=n.boot, interval=interval, ci=ci, censor=censor, path.to.results=path.to.results)
    saveRDS(conf.int, file.path(dir, "confint.rds"))

    if (!file.exists(track.results)){
        write.table(conf.int, file=track.results, append=TRUE, sep=",", row.names=FALSE, col.names=TRUE)
    } else{
        write.table(conf.int, file=track.results, append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
    }

    print((idx/num.trials)*100)
    print("finishing a ci")
  }
}

main.AssembleData(path.to.data, interval, nullsamples, censor, outdir)

print("onto fitting of curves")

main.compute(interval, censor, outdir, run_sim)
print("onto confidence intervals")

main.ci(interval=interval, ci=95, n.boot=n.boot, censor=censor, path.to.results=outdir)




























#end
