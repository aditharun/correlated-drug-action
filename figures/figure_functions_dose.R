#dose space figure helper functions
#I think i only need outlier.finder and ostt but can't be wrong to have it all 


outlier.finder <- function(item.lm, n){
  cutoff <- qt(1 - 0.05 / (2*n), (n - 3))
  outliers <- which(abs(unname(rstudent(item.lm))) > cutoff)
  return(outliers)
}
hill <- function(t, k, n){
  out <- 1 / (1 + ((t/k)^n))
  return(out)
}
fit <- function(dB, dA, kB, kA, nA, nB, D_A, D_B, p){
  if (p >= 0 ){
    if (dB < corresponding_B(dA, kB, nB, kA, nA)){
        val <-  ((D_A*D_B) + ((D_A)*(1 - D_B)*(p)))
      }
      else {
        val <-  ((D_A*D_B) + ((D_B)*(1 - D_A)*(p)))
      }
  } else{

    if (dB < negative_B(dA, kB, nB,kA, nA)){
      val <- ((D_A + D_B - 1)*(-1*p))+((D_A*D_B)*(1+p))
    } else{
      val <- (D_A*D_B)*(p + 1)
    }
  }

  return(val)
}
negative_B <- function(doseA, kB, nB, kA, nA){
  return ( ((kA / doseA)^(nA / nB))*(kB) )
}
corresponding_B <- function(doseA, kB, nB, kA, nA){
  return ( ((doseA / kA)^(nA / nB))*(kB) )
}
optimize <- function(dosingdrugB, dosingdrugA, true, kB, kA, nA, nB, rho, optimal, w.outliers, outliers){
  out.matrix <- as.data.frame(matrix(1, nrow=dosingdrugA, ncol=dosingdrugB))
  B.names <- true %>% select(-row.id) %>% colnames() %>% as.numeric()
  A.names <- true$row.id %>% as.numeric()
  for (vert in 1:dosingdrugB){
    for (horiz in 1:dosingdrugA){

      D_B <- true[nrow(true), (horiz+1)]
      D_A <- true[vert, 2]
      dA <- A.names[horiz]
      dB <- B.names[vert]

      out.matrix[vert, horiz] <- fit(dB, dA, kB, kA, nA, nB, D_A, D_B, rho)
    }
  }
  out.matrix <- out.matrix / 100

  true.values <- true %>% select(-row.id)

  if (optimal){

    if (w.outliers){
      if (!identical(outliers, integer(0))){
        cells.to.remove <- which(!is.na(true.values), arr.ind=TRUE) %>% as_tibble() %>% slice(outliers)
        for (r in 1:length(outliers)){
          true.values[unlist(unname(as.list(cells.to.remove[r,1]))), unlist(unname(as.list(cells.to.remove[r,2])))] <- out.matrix[unlist(unname(as.list(cells.to.remove[r,1]))), unlist(unname(as.list(cells.to.remove[r,2])))]
        }
      }
    }

    rmsd <- colSums((true.values - out.matrix)^2) %>% unname() %>% sum() %>% sqrt()
    return(rmsd)
  }
  return(out.matrix)
}
outlier.finder <- function(item.lm, n){
  cutoff <- qt(1 - 0.05 / (2*n), (n - 3))
  outliers <- which(abs(unname(rstudent(item.lm))) > cutoff)
  return(outliers)
}
aic <- function(x,y,p){
  param <- p*2
  SSE <- sum((x-y)^2)
  n <- length(x)
  out <- ((log(SSE, base=10) - log(n, base=10))*(-n/2)) + param
  return(out)
}
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
