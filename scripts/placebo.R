library(tidyverse)
library(readxl)

set.seed(1498)

#Trial IDs from collected clinical trial data that have the necessary data to impute placebo curve
placebo.trial.numbers <- c(14, 2,4, 8)


indir <- getwd()

#Directory to store results
outdir <- file.path(indir, "..", "placebo.results")

#Directory to get clinical trial data
path.to.data <- file.path(indir, "..", "raw-data", "clinical_trials")

source("randomization.R")
source("estimation_methods.R")
source("functions.R")

#file with how many patients are in a certain trial
path.to.pt <- file.path(path.to.data, "patient_numbers.xlsx")

#parameters

#for linearly interpolating PFS curve data, what interval (in months) should we use
interval <- 0.02

#how many samples for generating the null distribution
nullsamples <- 50000

#censor the last .03 (right tail) PFS data because its usually very noisy
censor <- 0.97

#to find the optimal hill curve parameters, to what number should we check?
limit <- 6

#in what intervals should we compute estimates for the hill curve using the parameters?
spacing.of.param <- 0.1


placebo.data.assembly <- function(name, interval, null.samples, path.to.pt, censor, outdir){

  files <- list.files(path=name, pattern="*.csv", full.names=TRUE, recursive=TRUE)
  files.short <- list.files(path=name, pattern="*.csv", full.names=FALSE, recursive=TRUE)
  file.w.placebo <- files.short[grep("PL", files.short)]
  non.placebo.file <- unlist(str_split(file.w.placebo, "\\+"))
  other.drug <- non.placebo.file[!grepl(".csv", non.placebo.file)]

  name.of.drug <- other.drug

  trial.name <- unlist(str_split(name, "/"))
  trial.name <- trial.name[length(trial.name)]
  print(trial.name)

  sheet <- sub(".*-", "", trial.name)

  other.drug.file <- files[grepl(paste0(other.drug, ".csv"), files.short)]
  combo.file <- files[grep("PL", files.short)]

  drug <- get.data(other.drug.file)
  comb <- get.data(combo.file)

  drug <- linear.approx(drug, interval)
  comb <- linear.approx(comb, interval)


  val <- min(tail(drug$x, n=1), tail(comb$x, n=1))
  drug$y <- drug$y[drug$x < val]
  comb$y <- comb$y[comb$x < val]
  comb$x <- comb$x[comb$x < val]
  tib <- as_tibble(data.frame(comb$x, comb$y, drug$y))
  tib <- tib %>% mutate_at(vars(-comb.x), ~./100)

  colnames(tib) <- c("time", "drugAB", "drugA")

  dir <- file.path(outdir, trial.name)
  null.file <- file.path(dir, "null.rds")
  data.file <- file.path(dir, "data.rds")

  patient.numbers <- read_excel(path.to.pt, sheet=sheet)
  num.pts.pl <- patient.numbers %>% pull(PL)
  num.pts.mono <- patient.numbers %>% pull(mono)

  attr(tib, "num.pts.mono") <- num.pts.mono

  if (!dir.exists(dir)){
    dir.create(dir, recursive=TRUE)
  }

  if (!file.exists(null.file)){
    create.null(null.file, num.pts.pl, null.samples, censor, tib, num.pts.mono)
  }

  saveRDS(tib, data.file)
}

hill <- function(t, k, n){
  y <- 1 / (1 + (t/k)^n)
  y
}

wrapper.model <- function(x, tib, censor, mono.pts){
  cutpt <- max(tib$time)*censor
  y <- tib$drugA
  z.hat <- x + y - (x*y)
  y.subsample <- sample2.from.cdf(tib$time, z.hat, n=mono.pts, interval=(1/mono.pts))
  colnames(tib)[which(colnames(tib)=="drugAB")] <- "pfs"
  ks.val <- ks.statistic(tib, y.subsample, cutpt)
  return(ks.val)
}

compute.placebo <- function(name, limit, spacing.of.param, censor, outdir){

  possible_rho <- seq(-1,1,0.01)
  data.path <- outdir
  dirs <- list.dirs(data.path, recursive=FALSE)
  case <- dirs[grepl(name, dirs)]

  data.path <- file.path(case, "data.rds")
  null.path <- file.path(case, "null.rds")

  tib <- readRDS(data.path)
  null <- readRDS(null.path)

  num.pts.mono <- attributes(tib)$num.pts.mono

  k.range <- seq(0.001, limit, spacing.of.param)
  n.range <- seq(0.001, limit, spacing.of.param)
  dimensions <- length(k.range)

  print("prep work done")

  ks.matrix <- as.data.frame(matrix(1, dimensions, dimensions))

  for (idx in 1:dimensions){

    k <- k.range[idx]
    vectors <- lapply(n.range, function(n) hill(tib$time, k, n))

    ks.vals <- unlist(lapply(vectors, function(b) wrapper.model(b, tib, censor, num.pts.mono) ))

    ks.matrix[,idx] <- ks.vals

    progress <- idx / dimensions
    print((progress*100))
  }

  print("heavy computation done")

  col.idx <- which.min(unname(apply(ks.matrix,MARGIN=2,min)))
  row.idx <- which.min(unname(apply(ks.matrix,MARGIN=1,min)))

  k.stable <- k.range[col.idx]
  n.stable <- n.range[row.idx]

  placebo <- hill(tib$time, k.stable, n.stable)

  ks.final <- wrapper.model(placebo, tib, censor, num.pts.mono)

  guess.comb <- placebo + tib$drugA - (placebo*tib$drugA)

  tib$placebo <- placebo
  tib$y.model <- guess.comb

  cutpt <- max(tib$time)*censor
  temp <- tib %>% filter(time < cutpt)

  model.pvalue <- (sum(null > ks.final) + 1) / length(null)

  print("p value generated")

  attr(tib, "p.model") <- model.pvalue
  attr(tib, "name") <- name

  out.file <- file.path(case, "raw.results.rds")
  saveRDS(tib, out.file)

  plot.file <- file.path(case, "fit.pdf")
  commentary <- paste0("K:", round(k.stable, 3), ", N:", round(n.stable, 3), ", p-val:", round(model.pvalue, 3))
  dat.long <- pivot_longer(tib, cols=-time,names_to="drug", values_to="pfs")
  pdf(plot.file)
  plot <- ggplot(dat.long, aes(x=time, y=pfs, color=drug)) + geom_line() + theme_bw() + xlim(0,max(dat.long$time)+3) + ggtitle(paste0(name, ", ", commentary))
  print(plot)
  dev.off()

  print("done")

  return(data.frame(name=name, p.value=model.pvalue, k=k.stable, n=n.stable))
}

trials <- list.dirs(path.to.data, recursive=FALSE, full.names=FALSE)
full.path.trials <- list.dirs(path.to.data, recursive=FALSE)
trial.number <- as.numeric(sub("-.*","",trials))

placebo.combinations <- full.path.trials[trial.number %in% placebo.trial.numbers]

for (value in placebo.combinations){
  placebo.data.assembly(value, interval, nullsamples, path.to.pt, censor, outdir)
  print("done with one")
}

trials.subset <- trials[trial.number %in% placebo.trial.numbers]

print(trials.subset)

track.results <- file.path(outdir, "placebo.summary.csv")

for (item in trials.subset){

  outkey <- compute.placebo(item, limit, spacing.of.param, censor, outdir)

  if (!file.exists(track.results)){
      write.table(outkey, file=track.results, append=TRUE, sep=",", row.names=FALSE, col.names=TRUE)
  } else{
      write.table(outkey, file=track.results, append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
  }

  print("end a trial")

}