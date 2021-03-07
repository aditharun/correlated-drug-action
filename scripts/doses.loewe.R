library(tidyverse)
library(readxl)
library(plot.matrix)
source("functions_dose_space.R")

#################################################################

#raw data directory
path.to.data <- file.path("..", "raw-data", "cell_line")

#output directory
outdir <- file.path("..", "results.cell.line", "doses.results.Bliss.Loewe")

#################################################################

#residual errors of p=0, p=p_opt - test for mean error suff different from 0 (t or z test)
#residual errors of p=0 and p=p_opt are suff different - AIC?
#determine if p=0 and p=p_opt are suff diff based on concordance/discordance of p-values

filenames <- list.files(path.to.data, full.names=TRUE, pattern="*.xlsx")
indir <- getwd()

fit <- function(dB, dA, kB, kA, nA, nB, D_A, D_B, p, A, B){
  if (p >= 0 ){
    if (loewe.f(dA, dB, kB, nB, kA, nA, p)){
        #val <-  ((D_A*D_B) + ((D_A)*(1 - D_B)*(p)))
        new.conc <- dB + p*f.i(dA, kA, nA, nB, kB)
        extra <- B$surv[which.min(abs(B$dose-new.conc))]
        val <- ((D_A * D_B) * (1-p) ) + extra
        label.conc <- "B"
      }
      else {
        new.conc <- dA + p*f(dB, kB, nB, nA, kA)
        label.conc <- "A"
        extra <- A$surv[which.min(abs(A$dose-new.conc))]
        val <- ((D_A*D_B)*(1-p) ) + extra
      }
  } else{

    if (dB < negative_B(dA, kB, nB,kA, nA)){
      val <- ((D_A + D_B - 1)*(-1*p))+((D_A*D_B)*(1+p))
    } else{
      val <- (D_A*D_B)*(p + 1)
    }

    new.conc <- NA
    label.conc <- NA
  }

  return(list(val=val, new.conc=new.conc, label.conc=label.conc))
}

optimize <- function(dosingdrugB, dosingdrugA, true, kB, kA, nA, nB, rho, optimal, w.outliers, outliers, A, B){
  out.matrix <- as.data.frame(matrix(1, nrow=dosingdrugA, ncol=dosingdrugB))
  eq.conc.matrix <- as.data.frame(matrix(1, nrow=dosingdrugA, ncol=dosingdrugB))
  label.matrix <- as.data.frame(matrix(1, nrow=dosingdrugA, ncol=dosingdrugB))

  B.names <- true %>% select(-row.id) %>% colnames() %>% as.numeric()
  A.names <- true$row.id %>% as.numeric()
  for (vert in 1:dosingdrugB){
    for (horiz in 1:dosingdrugA){

      D_B <- true[nrow(true), (horiz+1)]
      D_A <- true[vert, 2]
      dA <- A.names[vert]
      dB <- B.names[horiz]

      items <- fit(dB, dA, kB, kA, nA, nB, D_A, D_B, rho, A, B)

      out.matrix[vert, horiz] <- items$val
      eq.conc.matrix[vert, horiz] <- items$new.conc
      label.matrix[vert, horiz] <- items$label.conc
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

    return(list(out=rmsd, eqconc=eq.conc.matrix, label=label.matrix))
  }
  return(list(out=out.matrix, eqconc=eq.conc.matrix, label=label.matrix))
}


each.combination <- function(i, j, filenames, outdir2){
    name <- filenames[i]
    out.folder <- gsub(".xlsx", "", basename(name))
    outdir <- file.path(outdir2, out.folder)

    if (!dir.exists(outdir)){
      dir.create(outdir, recursive=TRUE)
    }

    x <- read_excel(name, sheet=j)
    x <- x  %>% drop_na()
    if ((dim(x)[1]) < 10){
      return(NA)
    }
    x <- x %>% arrange(desc(Drug1Conc), desc(Drug2Conc))

    dosingdrugA <- length(unique(x$Drug1Conc))
    dosingdrugB <- length(unique(x$Drug2Conc))

    drugA.name <- unique(x$Drug1)
    drugB.name <- unique(x$Drug2)

    combinedname <- paste0(drugA.name, " + ", drugB.name)

    values <- x$ObservedInhibition
    true <- t(matrix(values, ncol=dosingdrugB, nrow=dosingdrugA))
    class(true) <- "numeric"
    true <- 100 - true
    colnames(true) <- as.character(unique(x$Drug2Conc))
    true <- as_tibble(true)
    true <- true %>% select(rev(colnames(true)))
    true <- true %>% mutate(row.id=as.character(unique(x$Drug1Conc))) %>% select(row.id, everything())

    respB <- true %>% slice(n()) %>% select(-row.id) %>% unlist(., use.names=FALSE)
    respA <- true %>% select(-row.id) %>% pull(1)

    concB <- true %>% select(-row.id) %>% colnames() %>% as.numeric()
    concA <- true$row.id %>% as.numeric()

    A <- approx(concA, respA, seq(min(concA), max(concA), 0.2))
    B <- approx(concB, respB, rev(seq(min(concB), max(concB), 0.2)))
    lm.A <- lm(y~x, A)$coefficients
    lm.B <- lm(y~x, B)$coefficients

    #EC50
    kA <- unname((50 - lm.A[1]) / (lm.A[2]))
    kB <- unname((50 - lm.B[1]) / (lm.B[2]))

    A <- data.frame(dose=A$x, surv=A$y) %>% mutate_at(vars(surv), ~./100) %>% as_tibble()
    B <- data.frame(dose=B$x, surv=B$y) %>% arrange(dose) %>% mutate_at(vars(surv), ~./100) %>% as_tibble()

    print("preprocess done")

    n.range <- seq(0.001, 20, 0.1)

    fit.A <- lapply(n.range, function(y) hill(A$dose, kA, y))
    rmsd.A <- unlist(lapply(fit.A , function(x) sum((x - A$surv)^2)))
    nA <- n.range[which.min(rmsd.A)]

    fit.B <- lapply(n.range, function(y) hill(B$dose, kB, y))
    rmsd.B <- unlist(lapply(fit.B , function(x) sum((x - B$surv)^2)))
    nB <- n.range[which.min(rmsd.B)]

    possible_rho <- seq(-1, 1, 0.01)

    rmsds.init <- unlist(lapply(possible_rho, function(x) optimize.wrapper(dosingdrugB, dosingdrugA, true, kB, kA, nA, nB, x, TRUE, FALSE, "NA", A, B, FALSE)))

    best.rho.init <- possible_rho[which.min(rmsds.init)]
    guess.init <- optimize.wrapper(dosingdrugB, dosingdrugA, true, kB, kA, nA, nB, best.rho.init, FALSE, FALSE, "NA", A, B, FALSE)


    print("init optimization done")

    true.values <- true %>% select(-row.id)
    unwrap.true <- unlist(unname(as.list(true.values)))
    unwrap.guess.init <- unlist(unname(as.list(guess.init)))

    n <- length(unwrap.true)
    unwrap.lm <- lm(unwrap.true ~ unwrap.guess.init)
    outlier <- outlier.finder(unwrap.lm , n)
    num.outliers <- length(outlier)

    print("outlier done")

    rmsds.final <- unlist(lapply(possible_rho, function(x) optimize.wrapper(dosingdrugB, dosingdrugA, true, kB, kA, nA, nB, x, TRUE, TRUE, outlier, A, B, FALSE)))

    best.rho.final <- possible_rho[which.min(rmsds.final)]

    guess.final <- optimize.wrapper(dosingdrugB, dosingdrugA, true, kB, kA, nA, nB, best.rho.final, FALSE, FALSE, "NA", A, B, FALSE)
    guess.0 <- optimize.wrapper(dosingdrugB, dosingdrugA, true, kB, kA, nA, nB, 0, FALSE, FALSE, "NA", A, B, FALSE)

    eqconc.matrix <- optimize.wrapper(dosingdrugB, dosingdrugA, true, kB, kA, nA, nB, best.rho.final, TRUE, TRUE, outlier, A, B, TRUE)

    concentrations <- eqconc.matrix[[1]]
    #rownames(concentrations) <- rownames(true)
    rownames(concentrations) <- true %>% pull(row.id)
    colnames(concentrations) <- as.numeric(colnames(true)[-1])

    id.concs <- eqconc.matrix[[2]]

    unwrap.guess0 <- unlist(unname(as.list(guess.0)))
    unwrap.guess.final <- unlist(unname(as.list(guess.final)))

    print("final optimization done")

    if (!identical(outlier, integer(0))){
      unwrap.guess0 <- unwrap.guess0[-outlier]
      unwrap.guess.final <- unwrap.guess.final[-outlier]
      unwrap.true <- unwrap.true[-outlier]
      n <- length(unwrap.true)
    }

    out.0 <- ostt(unwrap.guess0, unwrap.true)
    out.final <- ostt(unwrap.guess.final, unwrap.true)

    comparator <- sum(out.0$p < 0.05) +  sum(out.final$p < 0.05)

    if (comparator==1){
      concordance <- FALSE
    } else{
      concordance <- TRUE
    }

    thresh.0 <- thresh(unwrap.guess0, unwrap.true)
    thresh.final <- thresh(unwrap.guess.final, unwrap.true)

    eob.z <- twosample.ztest(unwrap.true, unwrap.guess0, unwrap.guess.final)
    eob.t <- twosample.ztest(unwrap.true, unwrap.guess0, unwrap.guess.final, TRUE)


    z.0 <- onesample.z.test(unwrap.guess0,unwrap.true)
    z.final <- onesample.z.test(unwrap.guess.final, unwrap.true)



    if (aic(unwrap.true,unwrap.guess0,1) < aic(unwrap.true,unwrap.guess.final,1)){
      print("EOB is better model")
    }


    p.optrho <- chisq(unwrap.true, unwrap.guess.final, n)
    p.zerorho <- chisq(unwrap.true, unwrap.guess0, n)

    eob.lm <- lm(unwrap.guess0 ~ unwrap.guess.final)
    cda.lm <- lm(unwrap.true ~ unwrap.guess.final)

    eob.slope.p <- p.value(eob.lm, n, slope=TRUE)
    eob.int.p <- p.value(eob.lm, n, slope=FALSE)
    cda.slope.p <- p.value(cda.lm, n, slope=TRUE)
    cda.int.p <- p.value(cda.lm, n, slope=FALSE)
    eob.coeff <- coefficients(eob.lm)
    eob.eq <- paste0("y = ", round(eob.coeff[2],2), "x ", round(eob.coeff[1],2))
    cda.coeff <- coefficients(cda.lm)
    cda.eq <- paste0("y = ", round(cda.coeff[2],2), "x ", round(cda.coeff[1],2))

    chart.file <- file.path(outdir, paste0(j, ".plot.pdf"))
    pdf(chart.file)

    plot(hist(out.0$res, 15), xlab="residual p=0 case w/ true values",main=combinedname )
    plot(hist(out.final$res, 15), xlab="residual p=p_opt w/ true values", main=combinedname)

    plot(density(out.0$res), xlab="residual p=0 case w/ true values",main=combinedname )
    plot(density(out.final$res), xlab="residual p=p_opt w/ true values", main=combinedname)

    plot(unwrap.guess0, unwrap.true,pch=16,col="blue", main=combinedname, xlab="rho=0", ylab="experimental values")
    abline(a=0, b=1)

    plot(unwrap.guess.final, unwrap.true,pch=16,col="red", main=combinedname, xlab="optimal rho", ylab="experimental values")
    abline(a=0, b=1)
    dev.off()

    error <- guess.final - true.values
    heatmap <- true
    coln <- dim(heatmap)[2]

    #notice that if rho < 0, not able to find equivalent dose
    if (best.rho.final >= 0){
      conc.map <- true
      id.map <- true
      conc.map[,2:coln] <- concentrations
      id.map[,2:coln] <- id.concs
      long.conc <- pivot_longer(conc.map, -row.id, names_to="drugB", values_to="eq.conc")
      long.id <- pivot_longer(id.map, -row.id, names_to="drugB", values_to="id")
      long.conc <- inner_join(long.conc, long.id)
      colnames(long.conc)[1] <- "drugA"

      commentary.conc <- "loewe: drug + equivalent amount of drug"

      plot.file <- file.path(outdir, paste0(j, ".conc.pdf"))
      pdf(plot.file)

      long.conc$eq.conc <- round(long.conc$eq.conc, 2)

      mylevels.x <- unique(long.conc$drugB)
      mylevels.y <- unique(long.conc$drugA)
      long.conc$drugB2 <- factor(long.conc$drugB, mylevels.x[order(as.numeric(mylevels.x))])
      long.conc$drugA2 <- factor(long.conc$drugA, mylevels.y[order(as.numeric(mylevels.y))])

      plot <- ggplot(long.conc, aes(drugB2, drugA2)) +  geom_tile(aes(fill = id, width=0.95, height=0.95), size=2) + theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA)
        )  + geom_text(aes(label=eq.conc))
      plot <- plot + xlab(drugB.name) + ylab(drugA.name) + ggtitle(commentary.conc) + scale_fill_discrete(name = "Reference", labels = c(drugA.name, drugB.name))
      print(plot)
      dev.off()
    }

    colormap <- heatmap
    colormap[,2:coln] <- "non-outlier"

    if (!identical(outlier, integer(0))){
      cells.to.remove <- which(!is.na(true.values), arr.ind=TRUE) %>% as_tibble() %>% slice(outlier)
      for (r in 1:length(outlier)){
        colormap[unlist(unname(as.list(cells.to.remove[r,1]))), unlist(unname(as.list(cells.to.remove[r,2])))+1] <- "outlier"
      }
    }

    heatmap[,2:coln] <- error
    long.color <- pivot_longer(colormap, -row.id, names_to="drugB", values_to="outlier_status")
    colnames(long.color)[1] <- "drugA"
    long <- pivot_longer(heatmap, -row.id, names_to="drugB", values_to="viability")
    colnames(long)[1] <- "drugA"
    long <- inner_join(long, long.color)
    colnames(long)[3] <- "error"
    commentary <- paste0("optimal rho: ", round(best.rho.final, 3))

    plot.file <- file.path(outdir, paste0(j, ".pdf"))
    pdf(plot.file)

    coloritems <- c("non-outlier" = "transparent", "outlier" = "black")
    #make things numerically ordered on both axes - could be reordering the data wrong!
    mylevels.x <- unique(long$drugB)
    mylevels.y <- unique(long$drugA)
    long$drugB2 <- factor(long$drugB, mylevels.x[order(as.numeric(mylevels.x))])
    long$drugA2 <- factor(long$drugA, mylevels.y[order(as.numeric(mylevels.y))])

    plot <- ggplot(long, aes(drugB2, drugA2)) + scale_fill_gradient2(name="Excess over CDA", low="darkblue", high="red", guide="colorbar") +  geom_tile(aes(fill = error, color=outlier_status, width=0.95, height=0.95), size=2) + theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA)
      ) + scale_color_manual(name="Outlier Status", values=coloritems)
    plot <- plot + xlab(drugB.name) + ylab(drugA.name) + ggtitle(commentary)
    print(plot)
    dev.off()

    print("heatmap done")

    important <- data.frame(final.rho=best.rho.final, init.rho=best.rho.init, outliers=num.outliers, t.final=out.final$p, t.0=out.0$p, concordance=concordance, thresh.0=thresh.0, thresh.final=thresh.final, eob.z=eob.z, eob.t=eob.t, z.0=z.0, z.final=z.final)

    container <- list(drugA=drugA.name, drugB=drugB.name, guess0=guess.0, guessfinal=guess.final, true=true, guessinit=guess.init, conc=concentrations, id.concs=id.concs)

    #concentrations matrix is NA if p < 0 b/c not posible to find equivalent concentration
    saveRDS(container, file.path(outdir, paste0("results", j, ".rds")))

    #if opt.p > 0.05 then good fit to data - aka cda

    return(important)
}


track.results <- file.path(outdir, "doses.results.csv")

for (k in 1:length(filenames)){
  num.sheets <- length( excel_sheets( filenames[k] ) )
  for (j in 1:num.sheets){
    out <- each.combination(k, j, filenames, outdir)

    if (!file.exists(track.results)){
        write.table(out, file=track.results, append=TRUE, sep=",", row.names=FALSE, col.names=TRUE)
    } else{
        write.table(out, file=track.results, append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
    }

    print("done with a combination")

  }
  progress <- k/length(filenames)
  print((progress*100))

}
