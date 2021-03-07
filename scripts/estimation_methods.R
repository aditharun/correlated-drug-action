#helper functions that help functions.R work


#generates survival curves for a set of rho's under model method
model.choose <- function(data, possible_rho){

  orig.data <- data %>% mutate(better=pmax(drugA,drugB),worse=pmin(drugA,drugB)) %>% select(better, worse)

  positive <- possible_rho[possible_rho >=0]
  negative <- possible_rho[possible_rho < 0]

  y.pos <- lapply(positive, function(j) orig.data$better + orig.data$worse*(1-orig.data$better)*(1-j))

  tstar <- data %>% mutate(neg=abs(drugA+drugB-1)) %>% slice(which.min(neg)) %>% pull(time)

  y.neg <- lapply(negative, function(j) data %>% mutate(comb=ifelse(time <= tstar, ((drugA + drugB - (drugA*drugB)) * (1+j)) - j, drugA + drugB - (drugA*drugB) * (1 + j))) %>% pull(comb))

  y.pos <- do.call(cbind, y.pos)
  y.neg <- do.call(cbind, y.neg)
  y <- cbind(y.neg, y.pos)
  y <- as_tibble(y)

  return(y)
}

#performs randomization with correlation
monitor.vector <- function(data, rho, error){

  out <- in.the.looking.glass(data$drugB, data, rho, error)

  progress <- rho*100
  print(progress)

  out
}

#generates survival curves for a set of rho's under simulation method
simulate.choose <- function(data, possible_rho, tib){

  zero.case <- in.the.looking.glass(data$drugB, data, 0, 0.01)

  positive <- possible_rho[possible_rho > 0]

  results <- lapply(positive, function(t) monitor.vector(data=data, rho=t, error=0.01))

  y.pos <- as_tibble(do.call(cbind, results))

  y.neg <- y.pos %>% arrange(desc(row_number()))

  y.neg <- as_tibble(t(y.neg))
  y.neg <- y.neg %>% arrange(desc(row_number()))
  y.neg <- as_tibble(t(y.neg))

  zero.case <- as_tibble(data.frame(z=zero.case))
  y <- bind_cols(y.neg, zero.case, y.pos)

  y.pmax <- map_dfr(y, pmax, data$drugA)

  y.sort <- map_dfr(y.pmax, sort)

  return(y.sort)
}

#computes rmsd for simulated survival curve and true combination
rmsd.helper.sim <- function(x, interval, tib, cutpt){

  test <- approx(x$time, x$pfs, seq(0, max(tib$time), interval) )
  test3 <- data.frame(time=test$x, pfs=test$y)
  test3 <- as_tibble(test3)
  na.vals <- test3 %>% mutate(k=row_number()) %>% filter(!is.na(pfs)) %>% slice(1, n()) %>% pull(k)

  end <- na.vals[2]
  dim.type <- dim(test3)[1]

  if (end!=dim.type){
    test3$pfs[ (na.vals[2]+1): (dim(test3)[1]) ] <- test3$pfs[na.vals[2]]
  }

  test3$pfs[1:(na.vals[1]-1)] <- 1

  tib$y <- test3$pfs
  tib <- tib %>% filter(time < cutpt)
  y.sq <- sum((tib$y - tib$drugAB)^2)

  y.sq

}

#computes p values using null as reference
p.value.generator <- function(y, null, method, mono.pts, cutpt, tib){

  if (method=="model"){
      y.subsample <- apply(y, 2, function(x) sample2.from.cdf(tib$time, x, n=mono.pts, interval=(1/mono.pts)))
      colnames(tib)[which(colnames(tib)=="drugAB")] <- "pfs"
      ks.vals <- unlist(lapply(y.subsample, function(x) ks.statistic(tib, x, cutpt)))
  }

  if (method=="sim"){
    best.mono.tail <- max(tail(tib$drugA, 1), tail(tib$drugB, 1))
    y.surv <- apply(y, 2, function(x) transform.vector.to.cdf(x, best.mono.tail))
    colnames(tib)[which(colnames(tib)=="drugAB")] <- "pfs"
    y.subsample <- lapply(y.surv, function(x) sample2.from.cdf(x$time, x$pfs, n=mono.pts, interval=(1/mono.pts), sim.prep=FALSE))
    ks.vals <- unlist(lapply(y.subsample, function(x) ks.statistic(tib, x,cutpt)))
  }

  ks.vals <- unname(ks.vals)

  p <- unlist(lapply(ks.vals, function(x) p.cal(null, x)))

  return(p)

}

#sample from the cdf - proper method!
sample2.from.cdf <- function(x, y, n, interval=interval, sim.prep=FALSE){
  #test.x <- x
  #test.y <- y
  #upper.limit <- min(y)
  y <- 1-y
  finer.res <- approx(x, y, seq(0, max(x), interval))
  y <- finer.res$y
  x <- finer.res$x
  y[is.na(y)] <- 0
  pts <- runif(n, min(y), max(y))
  x.pts <- unlist(lapply(pts, function(t) which.min(abs(t-y))))
  scale <- max(y)
  samples <- x[x.pts]
  f <- ecdf(samples)
  cdf <- f(samples)

  samples.sort <- sort(samples)
  surv <- (length(samples.sort) - seq_along(samples.sort)) / length(samples.sort)
  #surv.scale <- surv * scale
  cdf.scale <- cdf * scale
  out <- data.frame(time=samples, pfs=cdf.scale)
  out <- out %>% arrange(time)

  out$pfs <- 1 - out$pfs

  out <- as_tibble(out)

  #pdf(paste0("test_",n,".pdf"))
  #plot(out$time, out$pfs, type="l")
  #lines(test.x, test.y, col="red")
  #dev.off()
  # number of points excluded
  #print(length(out$time) - sum(out$pfs >= upper.limit))
  #out <- out %>% filter(pfs >= upper.limit)
  if (sim.prep){
    return(out$time)
  }

  return(out)

}

#turn a vector into a dataframe of values and cdf values
transform.vector.to.cdf <- function(x, best.mono.tail){
  f <- ecdf(x)
  cdf <- f(x)
  g <- data.frame(x=x, y=1-cdf)
  g <- g %>% arrange(x)
  g$y[g$y < best.mono.tail] <- best.mono.tail
  colnames(g)  <- c("time", "pfs")

  g
}

#computes ks statistic between two survival curves
ks.statistic <- function(object, comparisons, cutpt){

  object <- object %>% filter(time < cutpt)
  comparisons <- comparisons %>% filter(time < cutpt)
  object.time <- object$time
  object.pfs <- object$pfs
  comparisons.pfs <- comparisons$pfs
  comparisons.time <- comparisons$time
  closest.pts <- object.pfs[unlist(lapply(comparisons.time, function(x) which.min(abs(x-object.time))))]
  #print(which.max(abs(comparisons.pfs - closest.pts)))
  max.diff <- max(abs(comparisons.pfs - closest.pts))

  return(max.diff)
}

#creates null distribution
create.null <- function(file, pts.AB, n, censor, tib, pts.mono){

  interval.ab <- 1/pts.AB
  cutpt <- max(tib$time)*censor
  interval.mono <- 1/pts.mono

  resamp.ab <- lapply(1:n, function(x) sample2.from.cdf(tib$time, tib$drugAB, pts.AB, interval=interval.ab))
  resamp.mono <- lapply(1:n, function(x) sample2.from.cdf(tib$time, tib$drugAB, pts.mono, interval=interval.mono))

  dist.model <- unlist(lapply(1:n, function(x) ks.statistic(resamp.ab[[x]], resamp.mono[[x]], cutpt)))

  #return(dist.model)

  saveRDS(dist.model, file)

}

#get raw data and read it in
get.data <- function(path){
  dat <- read.csv(file=path, header=FALSE)
  colnames(dat) <- c("time","pfs")
  dat <- as_tibble(dat)
  dat %>% arrange(time)
}

#linear approximation of data
linear.approx <- function(data, interval){
  res <- approx(data$time, data$pfs, seq(0, tail(data$time, n=1), interval))
  res$y[which(is.na(res$y))] <- 100
  return(list(x=res$x, y=res$y))
}

# p value computation
p.cal <- function(null, x){
  out <- (sum(null > x) + 1) / length(null)
  out
}

# create final best rho curve
generate.estimate.curve <- function(tib, times, rho, cutpt, method){

  if (method=="sim"){
    out.sim <- in.the.looking.glass(times$drugB, times, rho, 0.006)
    if (rho < 0){
      out.sim <- rev(out.sim)
    }
    times.ab <- sort(pmax(out.sim, times$drugA))
    best.mono.tail <- max(tail(tib$drugA, 1), tail(tib$drugB, 1))
    out <- transform.vector.to.cdf(times.ab, best.mono.tail)
  }

  if (method=="model"){
    orig.data <- tib %>% mutate(better=pmax(drugA,drugB),worse=pmin(drugA,drugB)) %>% select(better, worse)
    if (rho < 0){
        tstar <- tib %>% mutate(neg=abs(drugA+drugB-1)) %>% slice(which.min(neg)) %>% pull(time)
        out.model <- tib %>% mutate(comb=ifelse(time <= tstar, ((drugA + drugB - (drugA*drugB)) * (1+rho)) - rho, drugA + drugB - (drugA*drugB) * (1+rho))) %>% pull(comb)
    } else{
        out.model <- orig.data$better + orig.data$worse*(1-orig.data$better)*(1-rho)
    }
    out <- data.frame(time=tib$time, pfs=out.model)
  }

  out <- out %>% filter(time < cutpt)

  return(out)

}

#wrapper around CI method for model
wrapper.approx.ci <- function(x, interval, time, drug){
  cut.pt <- max(time)
  #mval <- max(drug)
  out <- approx(x$time, x$pfs, seq(0, cut.pt , interval))
  out <- data.frame(time=out$x, pfs=out$y)
  out <- as_tibble(out)


  if (!identical(which(is.na(out)), integer(0))){

    out2 <- out %>% mutate(k=row_number()) %>% filter(is.na(pfs))

    out2 <- out2 %>% mutate(greater=(time > max(x$time))) %>% mutate(less=(time < min(x$time)))

    na.last <- out2 %>% filter(greater==TRUE) %>% pull(k)

    na.1 <- out2 %>% filter(less==TRUE) %>% pull(k)

    if (length(na.1) > 0){
      out$pfs[na.1] <- 1
    }
    if (length(na.last) > 0){
      out$pfs[na.last] <- min(x$pfs)
    }

    print("na values detected and fixed")


  }


  #na.vals <- out %>% mutate(k=row_number()) %>% filter(!is.na(pfs)) %>% slice(1, n()) %>% pull(k)
  #out$pfs[1:(na.vals[1]-1)] <- 1

    #end <- na.vals[2]
    #dim.type <- dim(out)[1]

  #  if (end!=dim.type){
  #    out$pfs[ (na.vals[2]+1): (dim(out)[1]) ] <- out$pfs[na.vals[2]]
  #  }

  #out$pfs[out$pfs < mval] <- mval
  return(out$pfs)
}

#find the rho with lowest rmsd
get.best.rho <- function(method, cutpt, y, mono.pts, tib, interval){

  if (method=="model"){
    tib <- tib %>% filter(time < cutpt)
    y <- y %>% filter(row_number() < (dim(tib)[1]))
    y.cs <- apply(y, 2, function(x) sqrt(sum((x-tib$drugAB)^2)))
    y.out <- unname(y.cs)
  }

  if (method=="sim"){
    best.mono.tail <- max(tail(tib$drugA, 1), tail(tib$drugB, 1))
    y.surv <- apply(y, 2, function(x) transform.vector.to.cdf(x, best.mono.tail))
    y.out <- unlist(lapply(y.surv, function(x) rmsd.helper.sim(x, interval, tib, cutpt)))
    y.out <- unname(y.out)
  }

  y.out

}

#boostrap data and find rho for that data
refit <- function(data, interval, possible_rho, progress, total, censor){

  n <- ceiling(length(data$time)*1.3)

  num.mono.pts <- attributes(data)$num.pts.mono
  num.ab.pts <- attributes(data)$num.pts.ab

  boot.A  <- sample2.from.cdf(data$time, data$drugA, n=n, interval=interval)
  boot.AB <- sample2.from.cdf(data$time, data$drugAB, n=n, interval=interval)
  boot.B <- sample2.from.cdf(data$time, data$drugB, n=n, interval=interval)

  bt.A <- wrapper.approx.ci(boot.A, interval, data$time, data$drugA)
  bt.B <- wrapper.approx.ci(boot.B, interval, data$time, data$drugB)
  bt.AB <- wrapper.approx.ci(boot.AB, interval, data$time, data$drugAB)

  time <- seq(0, max(data$time), interval)

  bt <- data.frame(drugA=bt.A, drugB=bt.B, drugAB=bt.AB, time=time)
  cutpt <- max(bt$time)*censor

  model.estimate <- model.choose(bt, possible_rho)

  rmsd.model <- get.best.rho("model", cutpt, model.estimate, num.mono.pts, bt, interval)
  rho.model <- possible_rho[which.min(rmsd.model)]

  progress.bar <- (progress/total)*100
  print(progress.bar)

  return(rho.model)
}

#get a certain quantile for an empirical distribution
get.quantile <- function(vector, quantile){
  return(quantile(vector, quantile))
}

#compute confidence interval master function
confidence.interval <- function(name, n.boot, interval, ci, censor=censor, path.to.results){
  dir <- path.to.results
  dirs <- list.dirs(dir, recursive=FALSE)

  case <- dirs[grepl(name, dirs)]
  basename <- tail(unlist(str_split(case, "/")), 1)
  data.path <- file.path(case, "data.rds")
  null.path <- file.path(case, "null.rds")

  data <- readRDS(data.path)
  data <- data$data
  null <- readRDS(null.path)

  print("loaded data successfully")

  possible_rho <- seq(-1,1,0.005)

  estimates <- unlist(lapply(1:n.boot, function(k) refit(data=data, interval=interval, possible_rho=possible_rho, progress=k, total=n.boot, censor=censor)))

  #no problems before here, inspect from here on

  print("finished boot")

  distr.model <- sort(estimates)
  bounds.ci <- (100-ci)/200

  lower <- get.quantile(distr.model, bounds.ci)
  upper <- get.quantile(distr.model, (1-bounds.ci))

  return(data.frame(name=name, lower=lower, upper=upper))
}
