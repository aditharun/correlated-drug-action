### Simulation Results #####

library(ggpubr)
library(tidyverse)
library(cowplot)
library(readxl)


#simulate w/ jitter and fit w/ model
#simulate w/ coin and fit w/ model
#show that rho is fraction of points on p=1 line in coin method


gen.curve <- function(params, y=y){
	  k <- params[1]
	  n <- params[2]
	  x <- k*((1/y) - 1)^(1/n)
	  return(x)
}

get_sim_data <- function(kA, nA, kB, nB){
	drugA_params <- c(kA, nA)
	drugB_params <- c(kB, nB)

	y <- seq(0,1,0.0005)

	dat <- as_tibble(data.frame(drugA=gen.curve(drugA_params, y=y), drugB=gen.curve(drugB_params, y=y), pfs=y))
	dat <- dat[-1,]

	dat
}

#window.swap or coin
	
sim_corr <- function(dat, rho, method){
	B <- find.correlated.vector(dat, rho=rho, param=0.2,stall=40,type=method,error=0.005)

	if (rho < 0){
	  B <- rev(B)
	}

	times <- data.frame(B_jumbled=B, A=dat$drugA) %>% as_tibble()


	#fraction of points on p=1 line
 	n.line <- cbind(dat, times) %>% mutate(b=B_jumbled==drugB) %>% pull(b) %>% sum()	
 	frac <- n.line / dim(dat)[1]

	AB <- rev(sort(pmax(times$B_jumbled, dat$drugA)))

	dat$drugAB <- AB

	list(dat=dat, frac=frac, rho=rho, method=method, jumbledB=B)
}


approximate <- function(dat, coln, interval){
	name <- which(colnames(dat)==coln)
	v <- approx(dat[,name] %>% pull(colnames(dat)[name]), dat$pfs, seq(0, max(dat[,name]), interval))
	v$y[which(is.na(v$y))] <- 1
	times  <- list(x=v$x, y=v$y)

	times
}

fit <- function(fdat, interval){
	A <- approximate(fdat, "drugA", interval) 
	B <- approximate(fdat, "drugB", interval) 
	AB <- approximate(fdat, "drugAB", interval) 

	val <- min(tail(A$x, n=1), tail(B$x, n=1), tail(AB$x, n=1))

	A$y <- A$y[A$x < val]
	B$y <- B$y[B$x < val]
	AB$y <- AB$y[AB$x < val]
	A$x <- A$x[A$x < val]
	tib <- as_tibble(data.frame(A$x, A$y, B$y, AB$y))
	colnames(tib) <- c("time", "drugA", "drugB", "drugAB")

	poss_rho <- seq(-1, 1, 0.01)
	model.estimate <- model.choose(tib, poss_rho)

	y.cs <- apply(model.estimate, 2, function(x) sqrt(sum((x-tib$drugAB)^2)))
	y.out <- unname(y.cs)
	rho.model <- poss_rho[which.min(y.out)]

	guess <- model.estimate[,which.min(y.out)]

	colnames(guess) <- "model"

	tib <- cbind(tib, guess)

	list(tib=tib, rho=rho.model)
}


#-----------------------------------------------------------------------------------------