#location of helper functions

helper.dir <- "../scripts"

library(ggpubr)
library(tidyverse)
library(cowplot)
library(readxl)

source(file.path(helper.dir, "randomization.R"))
source(file.path(helper.dir, "estimation_methods.R"))
source(file.path(helper.dir, "functions.R"))
source(file.path(helper.dir, "simulation.R"))
source("figure_functions_clinical.R")


###########################################
axis_text_size = 15
legend_size = 14
title_size = 18
facet_text_size = 16
label_size = 25

dir.to.clin <- "../raw-data"
dir.to.clin.results <- "../results.clinical"

xlab.wrap <- 60
ylab.wrap <- 34
#############################################




#create sim data from ka kb na nb
dat <- get_sim_data(6, 3, 5, 0.9)

#figure for fraction lies on line
fraction_figure <- function(dat, axis_text_size, legend_size, title_size, facet_text_size, label_size){
	rhos <- seq(0, 1, 0.01)
	fracs <- unlist(lapply(rhos, function(x) sim_corr(dat, x, "coin")$frac))

	g <- data.frame(rho=rhos, fraction=fracs) %>% as_tibble()

	pval <- summary(lm(rho~fraction, g))$coefficients[2,4]
	pval <- format(pval, scientific=TRUE, digits=3)
	rsq <- round(summary(lm(rho~fraction, g))$r.squared, 3)

	panel_rhos <- c(0, 0.25, 0.5, 0.75)

	panels <- lapply(panel_rhos, function(x) sim_corr(dat, x, "coin"))
	int.panels <- lapply(panels, function(x) data.frame(B=x$jumbledB, A=x$dat$drugA, unalteredB=x$dat$drugB, rho=x$rho) %>% as_tibble() %>% mutate(online=B==unalteredB) %>% select(-unalteredB))
	int.panels <- do.call(rbind, int.panels)

	int.panels$rho <- paste0("Spearman Correlation: ", int.panels$rho)


	mytheme <- theme(
		panel.background = element_rect(fill = "transparent"), # bg of the panel
		panel.grid.major = element_blank(), # get rid of major grid
		panel.grid.minor = element_blank(), # get rid of minor grid
	) + theme(axis.line = element_line(color="black", size = 0.75)) + theme(axis.text.x = element_text(size = axis_text_size),
				axis.text.y = element_text(size = axis_text_size))  + theme(
	axis.title.x = element_text(size = title_size),
		axis.title.y = element_text(size = title_size)) + theme(plot.title = element_text(size = title_size))
	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	legend <- theme(legend.justification = 'left', legend.position="bottom", legend.title = element_blank(), legend.key = element_rect(colour = "transparent", fill = "white"), legend.text=element_text(size=legend_size), )
	
	plot <- ggplot(int.panels, aes(x=A, y=B, color=online)) + geom_point() + facet_wrap(~rho, scales="free") + scale_color_manual(name="", labels=c("not on corr = 1 curve", "on corr = 1 curve"), values=cbbPalette) + guides(colour = guide_legend(nrow = 1)) + legend + mytheme + theme(strip.text=element_text(size=facet_text_size), strip.background = element_blank())
	plot <- plot + xlab("Simulated Survival Times Under Drug A") + ylab("Simulated Survival Times Under Drug B")

	plot2 <- ggplot(g, aes(x=rho, y=fraction)) + geom_point() + mytheme + xlab("Set Value of \u03c1 in Coin Simulation") + ylab("Fraction of points on the \u03c1 = 1 curve") + geom_abline(slope=1, intercept=0, color="blue", linetype="dashed", size=1)
	plot2 <- plot2 + annotate(geom = 'text', label = paste0("   R-squared: ", rsq, ", p-value: ", pval) , x = -Inf, y = -0.15, hjust = 0, vjust = 1, size=4.7)
	
	fig <- plot_grid(plot, plot2, labels="AUTO", label_size=label_size, ncol=2)

	ggsave(fig, device=cairo_pdf, filename="supplement/coin_rho_interpretation.pdf", width=20, height=10, units="in", dpi=320)
}

fraction_figure(dat, axis_text_size, legend_size, title_size, facet_text_size, label_size)

#figure for predict with jitter fit with model and predict with coin and fit with model
fit.test <- function(dat, rho2, method){
	res <- sim_corr(dat, rho2, method)
	fdat <- res$dat
	guess <- fit(fdat, 0.01)
	p <- guess$tib %>% mutate(rho=res$rho) %>% select(time, rho, drugAB, model) %>% pivot_longer(-c(time, rho), names_to="type", values_to="val")
	list(p=p, rho=rho2, guess.rho=guess$rho)
}

plot.fit.test <- function(j, coin=TRUE, title_size, axis_text_size, legend_size, facet_text_size){
	mytheme <- theme(
		panel.background = element_rect(fill = "transparent"), # bg of the panel
		panel.grid.major = element_blank(), # get rid of major grid
		panel.grid.minor = element_blank(), # get rid of minor grid
	) + theme(axis.line = element_line(color="black", size = 0.75)) + theme(axis.text.x = element_text(size = axis_text_size),
				axis.text.y = element_text(size = axis_text_size))  + theme(
	axis.title.x = element_text(size = title_size),
		axis.title.y = element_text(size = title_size)) + theme(plot.title = element_text(size = title_size))
	cbbPalette <- c("#000000",  "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	legend <- theme(legend.justification = 'left', legend.position="bottom", legend.title = element_blank(), legend.key = element_rect(colour = "transparent", fill = "white"), legend.text=element_text(size=legend_size), )
	
	l <- do.call(rbind, lapply(j, function(x) x$p))
	l$rho <- paste0("Spearman Correlation: ", l$rho)
	p <- ggplot(l, aes(x=time, y=val, color=type, linetype=type)) + facet_wrap(~rho, scales="free") + geom_line(size=1) + mytheme + legend  + guides(colour = guide_legend(nrow = 1), linetype=FALSE) + theme(strip.text=element_text(size=facet_text_size), strip.background = element_blank()) 
	p <- p + xlab("Survival Time") + ylab("Progression Free Survival")

	if (coin){
		p <- p + scale_color_manual(name="", labels=c(paste0("'Coin' Simulated Combination"), paste0("Optimal CDA Model Estimate", j$guess.rho)), values=cbbPalette)	
	} else{
		p <- p + scale_color_manual(name="", labels=c(paste0("'Window Swap' Simulated Combination"), paste0("Optimal CDA Model Estimate", j$guess.rho)), values=cbbPalette)
	}	
	p
}

master_sim_comp_figure <- function(method, dat, title_size, axis_text_size, legend_size, facet_text_size, label_size){

	mytheme <- theme(
		panel.background = element_rect(fill = "transparent"), # bg of the panel
		panel.grid.major = element_blank(), # get rid of major grid
		panel.grid.minor = element_blank(), # get rid of minor grid
	) + theme(axis.line = element_line(color="black", size = 0.75)) + theme(axis.text.x = element_text(size = axis_text_size),
				axis.text.y = element_text(size = axis_text_size))  + theme(
	axis.title.x = element_text(size = title_size),
		axis.title.y = element_text(size = title_size)) + theme(plot.title = element_text(size = title_size))
	cbbPalette <- c("#000000",  "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	legend <- theme(legend.justification = 'left', legend.position="bottom", legend.title = element_blank(), legend.key = element_rect(colour = "transparent", fill = "white"), legend.text=element_text(size=legend_size), )
	

	rho_list <- c(-0.5, 0, 0.5, 0.75)
	itemize <- lapply(rho_list, function(x) fit.test(dat, x, method))

	if (method=="coin"){
		p <- plot.fit.test(itemize, TRUE,title_size, axis_text_size, legend_size, facet_text_size )
	} else{
		p <- plot.fit.test(itemize, FALSE, title_size, axis_text_size, legend_size, facet_text_size)
	}

	rhos <- seq(-1, 1, 0.02)
	vals <- (lapply(rhos, function(x) fit.test(dat, x, method)[c(2,3)] %>% unlist() %>% unname()))
	v <- do.call(rbind, vals) %>% as_tibble()
	colnames(v) <- c("true.rho", "guess.rho")

	pval <- summary(lm(true.rho~guess.rho, v))$coefficients[2,4]
	pval <- format(pval, scientific=TRUE, digits=3)
	rsq <- round(summary(lm(true.rho~guess.rho, v))$r.squared, 3)

	if (rsq == 1){
		rsq.int <- summary(lm(true.rho~guess.rho, v))$r.squared
		rsq <- format(rsq.int, scientific=TRUE, digits=4)
	} 

	if (method=="coin"){
		plot2 <- ggplot(v, aes(x=true.rho, y=guess.rho)) + geom_point(size=1.3) + mytheme + xlab(paste0("Set Value of \u03c1 in Coin Simulated Combination")) + ylab("Optimal Estimate of \u03b1 with Temporal CDA Model") + geom_abline(slope=1, intercept=0, color="blue", linetype="dashed", size=1)
	} else{
		plot2 <- ggplot(v, aes(x=true.rho, y=guess.rho)) + geom_point(size=1.3) + mytheme + xlab(paste0("Set Value of \u03c1 in Window Swap Simulated Combination")) + ylab("Optimal Estimate of \u03b1 with Temporal CDA Model") + geom_abline(slope=1, intercept=0, color="blue", linetype="dashed", size=1)
	}
	plot2 <- plot2 + annotate(geom = 'text', label = paste0("     R-squared = ", rsq, ", p-value = ", pval) , x = -Inf, y = -1.02, hjust = 0, vjust = 1, size=5)
	fig <- plot_grid(p, plot2, labels="AUTO", label_size=label_size, ncol=2)

	ggsave(fig, device=cairo_pdf, filename=paste0("supplement/model_is_", method, ".pdf"), width=20, height=10, units="in", dpi=320)
}


master_sim_comp_figure("coin", dat,title_size, axis_text_size, legend_size, facet_text_size, label_size)
master_sim_comp_figure("window.swap", dat,title_size, axis_text_size, legend_size, facet_text_size, label_size)



sim.results <- function(kB, nB, kA, nA, rho, axis_text_size, label_size, facet_text_size, title_size, legend_size, xlab.wrap, ylab.wrap){

	dat <- get_sim_data(kA, nA, kB, nB)

	mytheme <- theme(
		panel.background = element_rect(fill = "transparent"), # bg of the panel
		panel.grid.major = element_blank(), # get rid of major grid
		panel.grid.minor = element_blank(), # get rid of minor grid
	) + theme(axis.line = element_line(color="black", size = 0.75)) + theme(axis.text.x = element_text(size = axis_text_size),
				axis.text.y = element_text(size = axis_text_size))  + theme(
	axis.title.x = element_text(size = title_size),
		axis.title.y = element_text(size = title_size)) + theme(plot.title = element_text(size = title_size))

	B.coin <- find.correlated.vector(dat, rho=rho, param=0.2,stall=40,type="coin",error=0.019)
	#B.double.walk <- find.correlated.vector(dat, rho=rho, param=0.2,stall=40,type="double_walk",error=0.015)
	B.window.swap <- find.correlated.vector(dat, rho=rho, param=0.2,stall=40,type="window.swap",error=0.019)
	#B.jitter <- find.correlated.vector(dat, rho=rho, param=0.2,stall=40,type="jitter",error=0.015)

	#print(rho - cor(B.coin, dat$drugA, method="spearman"))

	if (rho < 0){
	  B.coin <- rev(B.coin)
	  #B.jitter <- rev(B.jitter)
	  #B.double.walk <- rev(B.double.walk)
		B.window.swap <- rev(B.window.swap)
	}

	#tA vs. tB curve
	#times.actual <- data.frame(Coin=B.coin, Double.Walk=B.double.walk, Jitter=B.jitter, Window.Swap=B.window.swap, A=dat$drugA)
	times.actual <- data.frame(Coin=B.coin, Window.Swap=B.window.swap, A=dat$drugA)
	times.actual <- pivot_longer(times.actual, -A, names_to="method", values_to="B")

	drugs.labs <- c("Coin", "Window Swap")
	names(drugs.labs) <- c("Coin",  "Window.Swap")
	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	#need more separation between different tyoes, should be visually nice

	tatb <- ggplot(times.actual, aes(x=A, y=B, color=method)) + geom_point()+ xlim(0, 30) + ylim(0 ,30) + scale_color_manual(values=cbbPalette) + ggtitle(paste0("Spearman Correlation Value: ", rho))  + facet_wrap(~method, scales="free", labeller=labeller(method=drugs.labs))
	tatb <- tatb + mytheme + xlab(stringr::str_wrap("Simulated Pt. Survival Times under Drug A", width=xlab.wrap)) + ylab(stringr::str_wrap("Simulated Pt. Survival Times under Drug B", width=ylab.wrap)) + theme(strip.background =element_rect(fill="transparent")) +  theme(strip.text = element_text(size=facet_text_size)) + theme(legend.position = "none")


	coin <- rev(sort(pmax(B.coin, dat$drugA)))
	#jitter <- rev(sort(pmax(B.jitter, dat$drugA)))
	#dwalk <- rev(sort(pmax(B.double.walk, dat$drugA)))
	wswap <- rev(sort(pmax(B.window.swap, dat$drugA)))


	dat$coin <- coin
	#dat$jitter <- jitter
	#dat$dwalk <- dwalk
	dat$wswap <- wswap

	dat <- dat %>% as_tibble() %>% select(-c(drugA, drugB))

	dat.long <- pivot_longer(dat, cols=-pfs,names_to="drug", values_to="time")
	legend <- theme(legend.justification = 'left', legend.position="bottom", legend.title = element_blank(), legend.key = element_rect(colour = "transparent", fill = "white"), legend.text=element_text(size=legend_size), )


	curves <- ggplot(dat.long, aes(x=time, y=pfs, color=drug)) + geom_line(size=1.2) + mytheme + xlab("Time") + ylab("PFS") + ggtitle("Survival Curves Based On Simulated Data") + legend + scale_color_manual(values=cbbPalette, labels=c("Coin", "Window Swap")) + xlim(0, 30)

	out <- plot_grid(tatb, curves, labels="AUTO", label_size=label_size)

	out
}

part1 <- sim.results(6, 2,5, 1, 0.25, axis_text_size, label_size, facet_text_size, title_size, legend_size, xlab.wrap, ylab.wrap)

rho <- 0.25
k <- 3

real.methods.comp <- function(rho, k, dir.to.clin, dir.to.clin.results,axis_text_size, label_size, facet_text_size, title_size, legend_size, xlab.wrap, ylab.wrap){

		mytheme <- theme(
		panel.background = element_rect(fill = "transparent"), # bg of the panel
		panel.grid.major = element_blank(), # get rid of major grid
		panel.grid.minor = element_blank(), # get rid of minor grid
	) + theme(axis.line = element_line(color="black", size = 0.75)) + theme(axis.text.x = element_text(size = axis_text_size),
				axis.text.y = element_text(size = axis_text_size))  + theme(
	axis.title.x = element_text(size = title_size),
		axis.title.y = element_text(size = title_size)) + theme(plot.title = element_text(size = title_size))
	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	legend <- theme(legend.justification = 'left', legend.position="bottom", legend.title = element_blank(), legend.key = element_rect(colour = "transparent", fill = "white"), legend.text=element_text(size=legend_size), )
	
	metares <- readRDS(file.path(dir.to.clin.results, "summary.rds"))

	censor <- 0.98
	newd <- file.path(dir.to.clin, "clinical_trials")
	dirs <- list.dirs(newd, recursive=FALSE, full.name=FALSE)
	idx <- dirs[k]


	file <- file.path(dir.to.clin.results, idx, "data.rds")
	fullname <- metares %>% filter(name==idx) %>% pull(fullname)

	A.name <- metares %>% filter(name==idx) %>% pull(A_name)
	B.name <- metares %>% filter(name==idx) %>% pull(B_name)

	xlab.A <- paste0("Pt. Survival Times Under ", A.name)
	ylab.B <- paste0("Pt. Survival Times Under ", B.name)

	data <- readRDS(file)
	tib <- data$data
	times <- data$sim.data

	cutpt <- max(tib$time)*censor
	num.mono.pts <- attributes(tib)$num.pts.mono
	num.ab.pts <- attributes(tib)$num.pts.ab

	model <- bands(tib, rho)
	b.model <- find.correlated.vector(times, rho=rho, param=0.2,stall=40,type="coin",error=0.019)

	b.sim <- in.the.looking.glass(times$drugB, times, rho, error=0.019)
	y.sim <- sort(pmax(b.sim, times$drugA))

	times.actual <- data.frame(model=b.model, sim=b.sim, A=times$drugA)
	times.actual <- pivot_longer(times.actual, -A, names_to="method", values_to="B")
	drugs.labs <- c("Model", "Simulation")
	names(drugs.labs) <- c("model",  "sim")
	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	tatb <- ggplot(times.actual, aes(x=A, y=B, color=method)) + geom_point() + scale_color_manual(values=cbbPalette) + ggtitle(paste0("Spearman Correlation Value: ", rho))  + facet_wrap(~method, scales="free", labeller=labeller(method=drugs.labs)) 
	tatb <- tatb + mytheme + xlab(stringr::str_wrap(xlab.A, width=xlab.wrap)) + ylab(stringr::str_wrap(ylab.B, width=ylab.wrap)) + theme(strip.background =element_rect(fill="transparent")) +  theme(strip.text = element_text(size=facet_text_size)) + theme(legend.position = "none")

	best.mono.tail <- max(tail(tib$drugA, 1), tail(tib$drugB, 1))
	sim.out <- transform.vector.to.cdf(y.sim, best.mono.tail)

	test <- approx(sim.out$time, sim.out$pfs, tib$time)
	test2 <- data.frame(time=test$x, 'Simulated Estimate'=test$y) %>% as_tibble()
	test2$idx <- 1:dim(test2)[1]
	intervals.na <- test2 %>% filter(!is.na(Simulated.Estimate)) %>% slice(1, n()) %>% pull(idx)

	if (intervals.na[1]!=1){
	test2$Simulated.Estimate[1:(intervals.na[1]-1)] <- test2$Simulated.Estimate[intervals.na[1]]
	}
	if (intervals.na[2]!=dim(test2)[1]){
	test2$Simulated.Estimate[(intervals.na[2]+1):dim(test2)[1]] <- test2$Simulated.Estimate[intervals.na[2]]
	}
	tib$sim <- test2$Simulated.Estimate
	tib$model <- model

	tib2 <- tib %>% select(c(model, sim, time))

	tib.long <- pivot_longer(tib2, cols=-time,names_to="drug", values_to="pfs")
	legend <- theme(legend.justification = 'left', legend.position="bottom", legend.title = element_blank(), legend.key = element_rect(colour = "transparent", fill = "white"), legend.text=element_text(size=legend_size), )

	library(stringi)
	title <- str_wrap(paste0("Survival Curves Based On Clinical Trial Data: ", fullname), 65)

	curves <- ggplot(tib.long, aes(x=time, y=pfs, color=drug)) + geom_line(size=1.2) + mytheme + xlab("Time") + ylab("PFS")  + legend + ggtitle(title) + scale_color_manual(values=cbbPalette, labels=c("Coin", "Window Swap"))

	out2 <- plot_grid(tatb, curves, labels=c("C", "D"), label_size=label_size)
	out2
}

#these sim methods are similar but real world noise makes them different considerably!


part2 <- real.methods.comp(rho, k, dir.to.clin, dir.to.clin.results, axis_text_size, label_size, facet_text_size, title_size, legend_size, xlab.wrap, ylab.wrap)

t <- plot_grid(part1, part2, nrow=2)
ggsave(t, device=cairo_pdf, filename="supplement/compare_randomize_methods.pdf", width=18, height=10, units="in", dpi=320)







