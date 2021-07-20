#helper functions for dose space analyses

hill <- function(t, k, n){
  out <- 1 / (1 + ((t/k)^n))
  return(out)
}

loewe.f <- function(doseA, doseB, kB, nB, kA, nA, p, A, B){
    f <-  ( kA * ((doseB / kB)^(nB / nA)))
    f.inv <-  ( ((doseA / kA)^(nA / nB))*(kB) )

    #uncomment for eq65 
    #modified_conc_A <- (f*p + doseA)
    #modified_conc_B <- (f.inv*p + doseB)

    #for eq66
    modified_conc_A <- (f + doseA)
    modified_conc_B <- (f.inv + doseB)

    #0.5 is the interval on which we will have resolution
    interpolateA <- approx(A$dose, A$surv, n=max(A$dose)/0.5)
    interpolateB <- approx(B$dose, B$surv, n=max(B$dose)/0.5)


    VA <- interpolateA$y[which.min(abs(interpolateA$x - modified_conc_A))]
    VB <- interpolateB$y[which.min(abs(interpolateB$x - modified_conc_B))]

    b <-  VA > VB
    return(b)
}

f.i <- function(doseA, kA, nA, nB, kB){
  f.inv <-  ( ((doseA / kA)^(nA / nB))*(kB) )
  return(f.inv)
}

f <- function(doseB, kB, nB, nA, kA){
  f <-  ( kA * ((doseB / kB)^(nB / nA)))
  return(f)
}

negative_B <- function(doseA, kB, nB, kA, nA){
  return ( ((kA / doseA)^(nA / nB))*(kB) )
}
corresponding_B <- function(doseA, kB, nB, kA, nA){
  return ( ((doseA / kA)^(nA / nB))*(kB) )
}

outlier.finder <- function(item.lm, n){
  cutoff <- qt(1 - 0.05 / (2*n), (n - 3))
  outliers <- which(abs(unname(rstudent(item.lm))) > cutoff)
  return(outliers)
}


optimize.wrapper <- function(dosingdrugB, dosingdrugA, true, kB, kA, nA, nB, x, y, z, a, A, B, conc, LOO=FALSE){
  res <- optimize(dosingdrugB, dosingdrugA, true, kB, kA, nA, nB, x, y, z, a, A, B, LOO)
  if (conc){
    return(list(res$eqconc, res$label))
  } else{
    return(res$out)
  }
}





######################## legacy functions - not in use anymore
chisq <- function(observed, expected, n){
  value <- sum(  (  ((observed - expected)^2) / expected)  )
  p <- pchisq(value, df=(n-1), lower.tail=FALSE)
  return(p)
}

ostt <- function(x, y){

  res <- x-y
  p <- t.test(res, mu=0, alternative="two.sided")$p.value

  return(list(res=res, p=p))

}


#problem here is that the entire distirbution could be skewed where 0 is the minimal value and the p value will come up huge!
thresh <- function(x, y){
  res <- x-y
  f <- ecdf(res)
  y <- f(res)
  finer <- approx(res, y, n=100)
  samples <- unlist(lapply(runif(1000), function(x) which.min(abs(x-finer$y))))
  samples <- finer$x[samples]
  p.value <- sum(samples >= 0) / length(samples)
  return(p.value)
}

#we can use a two sample z test to determine if the two distributions are significantly differnt from each other

twosample.ztest <- function(x,y,z, t=FALSE){
  res1 <- x-y
  res2 <- x-z
  z.statistic <- abs ( ( mean(res1) - mean(res2) ) /  sqrt( var(res1) + var(res2) ) )
  z.statistic <- z.statistic * -1
  p.value <- pnorm(z.statistic)*2
  if (t){
    p.value <- t.test(res1, res2, alternative="two.sided")$p.value
  }
  return(p.value)
}

onesample.z.test <- function(x,y){
  res <- x-y
  n <- length(res)
  statistic <- abs ( mean(res) / ( sd(res) / sqrt(n) ) )
  statistic <- statistic*-1
  p.value <- pnorm(statistic)*2
  return(p.value)
}

aic <- function(x,y,p){
  param <- p*2
  SSE <- sum((x-y)^2)
  n <- length(x)
  out <- ((log(SSE, base=10) - log(n, base=10))*(-n/2)) + param
  return(out)
}


#testing hypothesis finding p value
#two sided t test because we are testing H_0 => B = 1 and B not = 1
p.value <- function(item.lm, n, slope=TRUE){
  t.bounds <- qt( 1 - 0.05/2, (n - 2))
  summary <- summary(item.lm)$coefficients

  if (slope){
    se <-  summary[2,"Std. Error"]
    coef <- summary[2,"Estimate"]
    true <- 1
  } else{
    se <-  summary[1,"Std. Error"]
    coef <- summary[1,"Estimate"]
    true <- 0
  }
  t <- abs (  (coef - true) / se )
  p <- 2*pt(t, (n-2), lower=FALSE)

  return(p)

}
