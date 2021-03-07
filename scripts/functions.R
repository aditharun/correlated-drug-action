#Functions that process a clinical trial

#does preprocessing and parsing of data and creates null distribution
assemble.data <- function(data, path.to.pt, interval, null.samples, censor, dir.out){

  sort.therapies <- grepl("Both", names(data))

  true.AB <- data[sort.therapies]
  true.A <- data[!sort.therapies][1]
  true.B <- data[!sort.therapies][2]

  basename <- sub("/.*", "", names(true.B))
  sheet <- sub(".*-", "", basename)

  A <- linear.approx(true.A[[1]], interval)
  B <- linear.approx(true.B[[1]], interval)
  AB <- linear.approx(true.AB[[1]], interval)

  val <- min(tail(A$x, n=1), tail(B$x, n=1), tail(AB$x, n=1))

  A$y <- A$y[A$x < val]
  B$y <- B$y[B$x < val]
  AB$y <- AB$y[AB$x < val]
  A$x <- A$x[A$x < val]
  tib <- as_tibble(data.frame(A$x, A$y, B$y, AB$y))
  tib <- tib %>% mutate_at(vars(-A.x), ~./100)
  colnames(tib) <- c("time", "drugA", "drugB", "drugAB")

  print("finished assembling data from raw data")

  patient.numbers <- read_excel(path.to.pt, sheet=sheet)
  num.pts.AB <- patient.numbers %>% pull(Both)

  #mono intrinsically has the minimum of the two monotherapy curves
  num.pts.mono <- patient.numbers %>% pull(mono)

  num.pts.AB <- ceiling(num.pts.AB * censor)
  num.pts.mono <- ceiling(num.pts.mono * censor)

  dir.out <- file.path(dir.out, basename)

  null.file <- file.path(dir.out, "null.rds")
  data.file <- file.path(dir.out, "data.rds")

  if (!dir.exists(dir.out)){
    dir.create(dir.out, recursive = TRUE)
  }

  #for simulation proper
  timesA <- sample2.from.cdf(tib$time, tib$drugA, n=(length(tib$time)*4), interval=0.001, sim.prep=TRUE)

  timesB <- sample2.from.cdf(tib$time, tib$drugB, n=(length(tib$time)*4), interval=0.001, sim.prep=TRUE)

  times <- data.frame(drugA=sort(timesA, decreasing=TRUE), drugB=sort(timesB, decreasing=TRUE))

  print("computed data for simulation estimation")

  attr(tib, "num.pts.mono") <- num.pts.mono
  attr(tib, "num.pts.ab") <- num.pts.AB

  out <- list(data=tib, sim.data=times)

  saveRDS(out, data.file)

  print("saved assembled data")


  if (!file.exists(null.file)){
    create.null(null.file, num.pts.AB, null.samples, censor, tib, num.pts.mono)
  }

  print("saved null distribution")
}

#where name is the full file path
#takes an already preprocessed clinical trial (tib assembled, and null distr created) and computes fit
compute <- function(name, interval, censor, possible_rho, run_sim){

  basename <- basename(name)

  data.path <- file.path(name, "data.rds")
  null <- file.path(name, "null.rds")

  data <- readRDS(data.path)
  null <- readRDS(null)

  tib <- data$data
  times <- data$sim.data

  cutpt <- max(tib$time)*censor
  num.mono.pts <- attributes(tib)$num.pts.mono
  num.ab.pts <- attributes(tib)$num.pts.ab

  model.estimate <- model.choose(tib, possible_rho)
  rmsd.model <- get.best.rho("model", cutpt, model.estimate, num.mono.pts, tib, interval)
  rho.model <- possible_rho[which.min(rmsd.model)]
  p.model <- p.value.generator(model.estimate, null, "model", num.mono.pts, cutpt, tib)
  model.chart <- data.frame(model.param=possible_rho, p=p.model)

  if (run_sim){
    simulated.estimate <- simulate.choose(times, possible_rho, tib)
    rmsd.sim <- get.best.rho("sim", cutpt, simulated.estimate, num.mono.pts, tib, interval)
    rho.sim <- possible_rho[which.min(rmsd.sim)]
    p.sim <- p.value.generator(simulated.estimate, null, "sim", num.mono.pts, cutpt, tib)
    sim.chart <- data.frame(sim.cor=possible_rho, p=p.sim)
  }

  plot.file <- file.path(name, "p_charts.pdf")
  pdf(plot.file)

  if (run_sim){
    p.chart.sim <- ggplot(sim.chart, aes(sim.cor, p)) + geom_point() + ggtitle(basename) + theme_bw()
  }

  p.chart.model <- ggplot(model.chart, aes(model.param, p)) + geom_point() + ggtitle(basename) + theme_bw()
  print(p.chart.model)

  if (run_sim){
    print(p.chart.sim)
  }

  dev.off()


  y.model <- generate.estimate.curve(tib, times, rho.model, cutpt, "model")
  tib.out <- tib %>% filter(time < cutpt)
  p.model <- p.model[which.min(rmsd.model)]

  if (run_sim){
    y.sim <- generate.estimate.curve(tib, times, rho.sim, cutpt, "sim")
    p.sim <- p.sim[which.min(rmsd.sim)]
    results <- list(y.sim=y.sim, y.model=y.model, tib=tib.out)
    output <- data.frame(name=basename, rho.sim=rho.sim, rho.model=rho.model, p.sim=p.sim, p.model=p.model)
  } else{
    results <- list(y.model=y.model, tib=tib.out)
    output <- data.frame(name=basename,  rho.model=rho.model, p.model=p.model)
  }

  out.file <- file.path(name, "results.rds")
  saveRDS(results, out.file)

  return(output)

}

#path to results.rds for trial
#does some preliminary plotting of results
plotter <- function(name, run_sim){

  basename <- basename(name)

  file <- file.path(name, "results.rds")
  data <- readRDS(file)
  model <- data$y.model
  tib <- data$tib
  model$drug <- "model"

  if (run_sim){
    sim <- data$y.sim
    sim$drug <- "sim"

    if (identical(which(sim$time==1), integer(0))){
      sim <- sim %>% add_row(time=0, pfs=1, drug="sim", .before=1)
    }
  }

  long <- pivot_longer(tib, cols=-time,names_to="drug", values_to="pfs")

  out.path <- file.path(name, "fit.pdf")
  pdf(out.path)
  if (run_sim){
    out <- ggplot(long, aes(x=time, y=pfs, color=drug)) + geom_line() + theme_bw() + xlim(0,max(long$time)+1.5) + ggtitle(basename) + geom_line(aes(x=time, y=pfs, color=drug), sim) + geom_line(aes(x=time, y=pfs, color=drug), model)
  } else{
      out <- ggplot(long, aes(x=time, y=pfs, color=drug)) + geom_line() + theme_bw() + xlim(0,max(long$time)+1.5) + ggtitle(basename) + geom_line(aes(x=time, y=pfs, color=drug), model)
  }
  print(out)
  dev.off()
}
