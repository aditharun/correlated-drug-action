#code for clinical trial and placebo figures 

library(ggpubr)
library(tidyverse)
library(cowplot)
library(readxl)
source("figure_functions_clinical.R")

renderSim <- FALSE

#for figure 1
trial.to.follow <- 15

axis_text_size = 16
title_font_size = 20
legend_font_size = 20
facet_font_size = 16
pvalue_size_cex = 6
label_font_size = 27

width.fig = 12
height.fig = 15



###############################################


#dimensions are too big to reasonably apply the wrap title rule of str_wrap(x, width=50) rule
#maybe have to drop the dimensions, and if so implement ^^
walkthru <- function(k, dir.to.clin.results, dir.to.helper.functions, axis_text_size, title_font_size, legend_font_size, facet_font_size, label_font_size, pvalue_size_cex, width.fig, height.fig){

	metares <- readRDS(file.path(dir.to.clin.results, "summary.rds"))
	dirs <- list.dirs(dir.to.clin.results, recursive=FALSE, full.name=FALSE)
	idx <- dirs[k]
	filtered.meta <- metares %>% filter(name==idx)


	pvalue <- round(filtered.meta$p.model, 3)
	if ((pvalue==0)|(pvalue==1)){
		pvalue <- format(filtered.meta$p.model, scientific=TRUE)
	}

	file <- file.path(dir.to.clin.results, idx, "results.rds")
	res <- readRDS(file)
	tib <- res$tib
	model <- res$y.model$pfs
	tib$model <- model

	opt.rho <- filtered.meta$rho.model
	upper.rho <- filtered.meta$upper
	lower.rho <- filtered.meta$lower

	plusminus <- round((upper.rho - lower.rho)/2, 3)

	lower <- bands(tib, filtered.meta$lower)
	upper <- bands(tib, filtered.meta$upper)

	model.phrase <- (paste0("Model (","\u03C1 = ", opt.rho, " ± ", plusminus, ")"))

	poss_rho <- seq(-1,1,0.2)
	cone <- lapply(poss_rho, function(x) data.frame(bands(tib, x)))
	cones <- do.call(cbind, cone)
	colnames(cones) <- poss_rho
	cones <- cones %>% as_tibble()

	random.path <- file.path(dir.to.helper.functions, "randomization.R")
	vector.path <- file.path(dir.to.helper.functions, "functions.R")
	estimate.path <- file.path(dir.to.helper.functions, "estimation_methods.R")
	source(random.path)
	source(estimate.path)
	source(vector.path)
	rhos <- c(-1, -0.5, 0, 0.5, 1)
	data <- tib %>% select(drugA, drugB)
	

	timesA <- sample2.from.cdf(tib$time, tib$drugA, n=(length(tib$time)*4), interval=0.001, sim.prep=TRUE)
	timesB <- sample2.from.cdf(tib$time, tib$drugB, n=(length(tib$time)*4), interval=0.001, sim.prep=TRUE)
	times <- data.frame(drugA=sort(timesA, decreasing=TRUE), drugB=sort(timesB, decreasing=TRUE))
	ta.tb <- lapply(rhos, function(x) find.correlated.vector(times, x, 40, 20, "coin", 0.01))
	ta.tb2 <- lapply(seq_along(rhos), function(x) stitch.ta.tb(times, ta.tb, x, rhos))
	roll.ta.tb <- do.call(rbind, ta.tb2)
	roll.ta.tb <- roll.ta.tb %>% as_tibble()


	cones$time <- tib$time
	cones$actual <- tib$drugAB
	cones$id <- "True Combination"
	cone.wide <- pivot_longer(cones, -c(time, actual, id), names_to ="correlation", values_to="survival")


	lower <- mono(lower)
	upper <- mono(upper)
	tib$drugAB <- mono(tib$drugAB)
	tib$drugA <- mono(tib$drugA)
	tib$drugB <- mono(tib$drugB)
	tib$model <- mono(tib$model)
	colnames(tib)[which(colnames(tib)=="drugA")] <- filtered.meta$A_name
	colnames(tib)[which(colnames(tib)=="drugB")] <- filtered.meta$B_name
	colnames(tib)[which(colnames(tib)=="drugAB")] <- "Combination"
	colnames(tib)[which(colnames(tib)=="model")] <- model.phrase

	tib <- pivot_longer(tib, -time, names_to="drug", values_to="survival")

	tib$lower <- NA
	tib$lower[tib$drug==model.phrase] <- lower
	tib$upper <- NA
	tib$upper[tib$drug==model.phrase] <- upper

	tib.naive <- tib %>% filter(drug!=model.phrase)

	tib$drug <- as.factor(tib$drug)
	tib.naive$drug <- as.factor(tib.naive$drug)

	orig.levels <- levels(tib$drug)
	comb <- which(orig.levels=="Combination")
	mod <- which(orig.levels==model.phrase)
	others <- setdiff(1:4, c(comb, mod))

	naive.levels <- levels(tib.naive$drug)
	comb.naive <- which(naive.levels=="Combination")
	others.naive <- setdiff(1:3, comb.naive)

	tib.naive$drug <- factor(tib.naive$drug, levels = c(naive.levels[comb.naive], naive.levels[others.naive[1]], naive.levels[others.naive[2]]), ordered = TRUE)
	tib$drug <- factor(tib$drug, levels = c(orig.levels[mod], orig.levels[comb], orig.levels[others[1]], orig.levels[others[2]]), ordered = TRUE)

	###start plotting


	mytheme <- theme(
		panel.background = element_rect(fill = "transparent"), # bg of the panel
		panel.grid.major = element_blank(), # get rid of major grid
		panel.grid.minor = element_blank(), # get rid of minor grid
	) + theme(axis.line = element_line(color="black", size = 0.75)) + theme(axis.text.x = element_text(size = axis_text_size),
				axis.text.y = element_text(size = axis_text_size))  + theme(
	axis.title.x = element_text(size = title_font_size),
		axis.title.y = element_text(size = title_font_size)) + theme(plot.title = element_text(size = title_font_size))

	plt.all <- ggplot(tib, aes(x=time, y=survival)) + geom_smooth(aes(ymin=lower, ymax=upper, color=drug), stat="identity", size=1.4) + mytheme + xlab("Time (months)") + ylab("PFS") 

	#+ ggtitle(filtered.meta$fullname)

	legend <- theme(legend.justification = 'left', legend.position="bottom", legend.title = element_blank(), legend.key = element_rect(colour = "transparent", fill = "white"), legend.text=element_text(size=legend_font_size), )

	plt.naive <- ggplot(tib.naive, aes(x=time, y=survival)) + geom_smooth(aes(ymin=lower, ymax=upper, color=drug), size=1.4, stat="identity") + mytheme + xlab("Time (months)") + ylab("PFS") 

	#+ ggtitle(filtered.meta$fullname)

	cbbPalette.all <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	cbbPalette.naive <- cbbPalette.all[-1]

	plt.naive <- plt.naive + legend + scale_color_manual(values=cbbPalette.naive)+ guides(colour = guide_legend(nrow = 1)) + theme(legend.key = element_rect(colour = "transparent", fill = "white"))
	plt.all <- plt.all + legend + scale_color_manual(values=cbbPalette.all)+ guides(colour = guide_legend(nrow = 1)) + theme(legend.key = element_rect(colour = "transparent", fill = "white"))

	cone.wide$correlation <- as.numeric(cone.wide$correlation)
	legend <- theme(legend.key = element_rect(colour = "transparent", fill = "white"), legend.text=element_text(size=legend_font_size), legend.title=element_text(size=(legend_font_size+1)) )

	cone.plt <- ggplot(cone.wide) + geom_line(aes(x=time, y=survival, group=correlation, color=correlation), size=1.4) + mytheme + scale_color_gradientn(name="Correlation", colours=rainbow(length(poss_rho))) + xlab("Time (months)") + ylab("PFS") + geom_line(aes(x=time, y=actual, linetype=id), size=1.4) +  scale_linetype_manual(values=c("dashed"),guide="legend", name="", labels=c("Obs.\nComb.")) + legend 
	#+ ggtitle("Cone of Possible Estimated Combination Survival Curves")


	#active development

	roll.ta.tb2 <- roll.ta.tb %>% mutate(bin=1:n())
	negative_ref <- roll.ta.tb %>% filter(id==-1)
	positive_ref <- roll.ta.tb %>% filter(id==1)

	nbin <- roll.ta.tb2 %>% left_join(negative_ref, by=c("tA"="tA", "tB"="tB")) %>% drop_na() %>% pull(bin)
	pbin <- roll.ta.tb2 %>% left_join(positive_ref, by=c("tA"="tA", "tB"="tB")) %>% drop_na() %>% pull(bin)

	roll.ta.tb2 <- roll.ta.tb2 %>% mutate(neg_line=FALSE, pos_line=FALSE)
	roll.ta.tb2$neg_line[nbin] <- TRUE
	roll.ta.tb2$pos_line[pbin] <- TRUE

	roll.ta.tb2 <- roll.ta.tb2 %>% mutate(color=ifelse(pos_line, "a", ifelse(neg_line, "b", "c")))

	plt.all <- plt.all + annotate(geom = 'text', label = paste0('   p-value = ', pvalue) , x = -Inf, y = 0.1, hjust = 0, vjust = 1, size=pvalue_size_cex)

	cbbPalette <- c( "#0072B2", "#D55E00", "#000000","#009E73", "#F0E442", "#CC79A7",  "#E69F00", "#56B4E9")

	roll.ta.tb2$id <- paste0("Correlation: ", roll.ta.tb2$id)

	tatb.plt <- ggplot(roll.ta.tb2, aes(x=tA, y=tB, color=color)) + geom_point() + facet_wrap(~id, scales="free") + mytheme + xlab(paste0("Simulated Survival Times \nfor ",filtered.meta$A_name, " Patients")) + ylab(paste0("Simulated Survival Times \nfor ",filtered.meta$B_name, " Patients")) + theme(strip.text=element_text(size=facet_font_size)) + theme(axis.title.x = element_text(size = title_font_size), axis.title.y = element_text(size = title_font_size)) + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) 

	legend <- theme(legend.justification = 'left', legend.position="bottom", legend.title = element_blank(), legend.key = element_rect(colour = "transparent", fill = "white"), legend.text=element_text(size=legend_font_size), )

	tatb.plt <- tatb.plt + scale_color_manual(name="", labels=c("on corr = 1 curve", "on corr = -1 curve", "neither"), values=cbbPalette) + guides(colour = guide_legend(nrow = 2)) + legend

	middle_row <- plot_grid(cone.plt, tatb.plt, labels=c("B", "C"), nrow=1, label_size=label_font_size)

	fig1 <- plot_grid(plt.naive, middle_row, plt.all, labels=c("A", "", "D"), label_size=label_font_size, ncol=1)

	ggsave(fig1, device=cairo_pdf, filename="figure1.pdf", width=width.fig, height=height.fig, units="in", dpi=320)

	return(fig1)
}

walkthru(trial.to.follow, "../results.clinical", "../scripts",axis_text_size, title_font_size, legend_font_size, facet_font_size, label_font_size, pvalue_size_cex, width.fig, height.fig)

#FIGURE RENDER ------------- panels of clinical trials figure ---------------
panels_clinical_trials <- function(dir.to.clin.results, sim.figs=FALSE, title_font_size, axis_text_size, pvalue_size_cex, legend_font_size=14, nwidth=40){

	metares <- readRDS(file.path(dir.to.clin.results, "summary.rds"))

	dirs <- list.dirs(dir.to.clin.results, recursive=FALSE, full.name=FALSE)

	#distribution of rho model values for clinical trials
	#it looks really shitty either as histogram or density

	#multiple panels clinical trial results
	store.plts <- vector("list", length(dirs))
	names(store.plts) <- dirs



	for (k in 1:length(dirs)){

		idx <- dirs[k]
		filtered.meta <- metares %>% filter(name==idx)

		reformat_comb1 <- "Irinotecan, Bevacizumab and Panitumumab in Advanced Colorectal Cancer"
		reformat_comb2 <- "Oxaliplatin, Bevacizumab and Panitumumab in Advanced Colorectal Cancer"

		if (filtered.meta$fullname == reformat_comb1 ){

			filtered.meta$fullname <- "Irinotecan, Bevacizumab & Panitumumab in Advanced Colorectal Cancer"
		}

		if (filtered.meta$fullname == reformat_comb2 ){
			filtered.meta$fullname <- "Oxaliplatin, Bevacizumab & Panitumumab in Advanced Colorectal Cancer"
		}

		rho.model <- round(filtered.meta$rho.model, 3)
		rho.upper <- round(filtered.meta$upper, 3)
		rho.lower <- round(filtered.meta$lower, 3)

		file <- file.path(dir.to.clin.results, idx, "results.rds")
		res <- readRDS(file)
		tib <- res$tib

		if (sim.figs){
			rho.sim <- round(filtered.meta$rho.sim, 3)
			pvalue <- round(filtered.meta$p.sim, 3)
			if ((pvalue==0)|(pvalue==1)){
				pvalue <- format(filtered.meta$p.sim, scientific=TRUE)
			}

			sim <- res$y.sim
			tib$drugAB <- mono(tib$drugAB)
			tib$drugA <- mono(tib$drugA)
			tib$drugB <- mono(tib$drugB)
			sim$pfs <- mono(sim$pfs)
			colnames(tib)[which(colnames(tib)=="drugA")] <- filtered.meta$A_name
			colnames(tib)[which(colnames(tib)=="drugB")] <- filtered.meta$B_name
			colnames(tib)[which(colnames(tib)=="drugAB")] <- "Combination"

			test <- approx(sim$time, sim$pfs, tib$time)
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

			model.phrase <- paste0("Simulation (\u03C1 = ", rho.sim, ")")
			colnames(tib)[which(colnames(tib)=="sim")] <- model.phrase

			tib <- pivot_longer(tib, -time, names_to="drug", values_to="survival")
			tib$drug <- as.factor(tib$drug)
			orig.levels <- levels(tib$drug)
			comb <- which(orig.levels=="Combination")
			mod <- which(orig.levels==model.phrase)
			others <- setdiff(1:4, c(comb, mod))
			tib$drug <- factor(tib$drug, levels = c(orig.levels[mod], orig.levels[comb], orig.levels[others[1]], orig.levels[others[2]]), ordered = TRUE)
			mytheme <- theme(
				panel.background = element_rect(fill = "transparent"), # bg of the panel
				panel.grid.major = element_blank(), # get rid of major grid
				panel.grid.minor = element_blank(), # get rid of minor grid
			) + theme(axis.line = element_line(color="black", size = 0.75)) + theme(axis.text.x = element_text(size = 11),
						axis.text.y = element_text(size = 11))  + theme(
			axis.title.x = element_text(size = 14),
				axis.title.y = element_text(size = 14)) + theme(plot.title = element_text(size = 14))

			plt <- ggplot(tib) + geom_line(aes(x=time, y=survival, color=drug), size=1.3) + mytheme + xlab("Time (months)") + ylab("PFS") + ggtitle(stringr::str_wrap(filtered.meta$fullname, width=nwidth))

			legend <- theme(legend.justification = 'left', legend.position="bottom", legend.title = element_blank(), legend.key = element_rect(colour = "transparent", fill = "white"), legend.text=element_text(size=10), )

			cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

			plt <- plt + legend + scale_color_manual(values=cbbPalette) + guides(colour = guide_legend(nrow = 2))
			plt <- plt + annotate(geom = 'text', label = paste0('   p-value = ', pvalue) , x = -Inf, y = 0.15, hjust = 0, vjust = 1, size=4.7)
			store.plts[[k]] <- plt

			next
		}

		pvalue <- round(filtered.meta$p.model, 3)
		if ((pvalue==0)|(pvalue==1)){
			pvalue <- format(filtered.meta$p.model, scientific=TRUE, digits=3)
		}

		model <- res$y.model$pfs
		tib$model <- model

		lower <- bands(tib, filtered.meta$lower)
		upper <- bands(tib, filtered.meta$upper)

		lower <- mono(lower)
		upper <- mono(upper)
		tib$drugAB <- mono(tib$drugAB)
		tib$drugA <- mono(tib$drugA)
		tib$drugB <- mono(tib$drugB)
		tib$model <- mono(tib$model)
		model.phrase <- (paste0("Model: ","\u03C1 = (", rho.lower, ", ", rho.upper, ")"))

		colnames(tib)[which(colnames(tib)=="drugA")] <- filtered.meta$A_name
		colnames(tib)[which(colnames(tib)=="drugB")] <- filtered.meta$B_name
		colnames(tib)[which(colnames(tib)=="drugAB")] <- "Combination"
		colnames(tib)[which(colnames(tib)=="model")] <- model.phrase

		tib <- pivot_longer(tib, -time, names_to="drug", values_to="survival")

		tib$lower <- NA
		tib$lower[tib$drug==model.phrase] <- lower
		tib$upper <- NA
		tib$upper[tib$drug==model.phrase] <- upper

		tib$drug <- as.factor(tib$drug)

		orig.levels <- levels(tib$drug)
		comb <- which(orig.levels=="Combination")

		mod <- which(orig.levels==model.phrase)
		others <- setdiff(1:4, c(comb, mod))

		tib$drug <- factor(tib$drug, levels = c(orig.levels[mod], orig.levels[comb], orig.levels[others[1]], orig.levels[others[2]]), ordered = TRUE)


		mytheme <- theme(
			panel.background = element_rect(fill = "transparent"), # bg of the panel
			panel.grid.major = element_blank(), # get rid of major grid
			panel.grid.minor = element_blank(), # get rid of minor grid
		) + theme(axis.line = element_line(color="black", size = 0.75)) + theme(axis.text.x = element_text(size = axis_text_size),
					axis.text.y = element_text(size = axis_text_size))  + theme(
		axis.title.x = element_text(size = title_font_size),
			axis.title.y = element_text(size = title_font_size)) + theme(plot.title = element_text(size = (title_font_size-1)))
		plt <- ggplot(tib, aes(x=time, y=survival)) + geom_smooth(aes(ymin=lower, ymax=upper, color=drug), stat="identity") + mytheme + xlab("Time (months)") + ylab("PFS") + ggtitle(stringr::str_wrap(filtered.meta$fullname, width=nwidth))
		legend <- theme(legend.justification = 'left', legend.position="bottom", legend.title = element_blank(), legend.key = element_rect(colour = "transparent", fill = "white"), legend.text=element_text(size=legend_font_size), )

		cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

		#x label tick marks should be integers not floats! SO USE THE BELOW CODE APPROPRIATELY
		#breaks.x <-
		#scale_x_continuous(breaks=c(1,3,7,10))
		plt <- plt + legend + scale_color_manual(values=cbbPalette) + guides(colour = guide_legend(nrow = 2))
		plt <- plt + annotate(geom = 'text', label = paste0('   p-value = ', pvalue) , x = -Inf, y = 0.05, hjust = 0, vjust = 1, size=pvalue_size_cex)

		store.plts[[k]] <- plt


	}

	if (sim.figs){
		saveRDS(store.plts, "SIM_panels_clinical_trials.rds")
	}  else{
		saveRDS(store.plts, "MODEL_panels_clinical_trials.rds")
		print("done")
	}

}


panels_clinical_trials("../results.clinical", FALSE, title_font_size, axis_text_size, pvalue_size_cex=4.9)

if (renderSim){
	panels_clinical_trials("../results.clinical", sim.figs=TRUE, title_font_size, axis_text_size, pvalue_size_cex, legend_font_size)
}

if (renderSim){

	store.sim.plts <- readRDS("SIM_panels_clinical_trials.rds")

	ranges <- c(1, seq(0, length(store.sim.plts), length(store.sim.plts)/3)[-1])

	for (x in 1:(length(ranges)-1)){
		temp <- store.sim.plts[ranges[x]:ranges[x+1]]
		a <- temp[[1]]
		b <- temp[[2]]
		c <- temp[[3]]
		d <- temp[[4]]
		e <- temp[[5]]
		f <- temp[[6]]
		temp.fig <- plot_grid(a,b,c,d,e,f, label_size=label_font_size, labels="AUTO", ncol=2)
		fname <- paste0("supplement/simulated_curves_", x, ".pdf")
		ggsave(temp.fig, filename=fname, width=13.5, height=20, units="in", dpi=320, device=cairo_pdf)
	}
}


### ------------------------------------------------------------------------------------------------------

#Already made with summary_clintrial.R in final_version_of_scripts folder
#### create csv file with model and simulation summary statistics for supplement

#dir.to.clin.results <- "../scripts/final_version_of_scripts/results.clinical.test"
#summ <- readRDS(file.path(dir.to.clin.results,"summary.rds"))
#summ$id <- as.numeric(gsub("([0-9]+).*$", "\\1", summ$name))
#summ <- summ %>% as_tibble() %>% select(-c(name, A_name, B_name)) %>% rename(trial=fullname, model.lower.ci=lower, model.upper.ci=upper)
#summ <- summ %>% arrange(id)
#summ <- summ %>% dplyr::rename(`Trial ID`=id)
#summ <- summ %>% select(`Trial ID`, trial, model.lower.ci, rho.model, model.upper.ci, p.model, rho.sim, p.sim)
#write.csv(summ, "supplement/trials_results.csv", row.names=FALSE)


#Dot plot w/ CI showing rho model values
dotplot <- function(dir.to.clin.results, axis_text_size, title_font_size, legend_font_size){

	p.cutoff <- .01

	metares <- readRDS(file.path(dir.to.clin.results, "summary.rds"))
	metares$id <- as.numeric(gsub("([0-9]+).*$", "\\1", metares$name))
	metares <- metares %>% arrange(desc(id))

	mytheme <- theme(
		panel.background = element_rect(fill = "transparent"), # bg of the panel
		panel.grid.major.x = element_blank(),
		panel.grid.major.y = element_line(colour="grey70", size=0.15, linetype="dashed"),# get rid of major grid
		panel.grid.minor = element_blank(), # get rid of minor grid
	) + theme(axis.line = element_line(color="black", size = 0.75)) + theme(axis.text.x = element_text(size = axis_text_size),
				axis.text.y = element_text(size = axis_text_size))  + theme(
	axis.title.x = element_text(size = title_font_size),
		axis.title.y = element_text(size = title_font_size)) + theme(plot.title = element_text(size = title_font_size))

	cbbPalette <- c("#000000", "#F0E442", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

	metares <- metares %>% arrange(desc(rho.model)) %>% mutate(id=1:n())

	plot <- ggplot(metares, aes(x=id, y=rho.model)) + geom_hline(yintercept=0, size=1, color="grey", linetype="dashed") + geom_point(aes(color="c"), size=1.5) + geom_errorbar(aes(x=id, ymin=lower, ymax=upper), color="black")  + ylab("Correlation of Model Estimate") + xlab("Trial ID") + mytheme  + scale_x_continuous(breaks=seq(1,18,1)) + coord_flip()
	legend <- theme(legend.justification = 'left', legend.position="bottom", legend.title = element_blank(), legend.key = element_rect(colour = "transparent", fill = "white"), legend.text=element_text(size=legend_font_size), )

	nonsig.rows <- metares %>% filter(p.model > p.cutoff) %>% pull(id)
	plot <- plot + geom_point(data = metares %>% filter(id %in% nonsig.rows), aes(x=id, y=rho.model, color="b"), size=2.5) + scale_color_manual(name="", labels=c("CDA (p-value > 0.01)", "Non-CDA (p-value < 0.01)"), values=c("purple", "black")) + legend + guides(colour = guide_legend(nrow = 2))
	
	#highlight y labels another option for labeling non-sig rows
	#plot + theme(axis.text.y = element_text(color=rep("red", 18)))
	#ggsave(plot, device="pdf", filename="../../figures/dotplot.pdf", width=7, height=9, units="in", dpi=320)
	return(plot)
}

#i is the # of the trial - either 14, 2, 4, 8 - maybe not 8 unclear whether label is gc or gemcitabine
# we could compute some sort of CI for placebo by shuffling A and AP in the script that generates results of placebo...
placebo.plot <- function(i, dir.to.placebo, off=FALSE, axis_text_size, legend_font_size, title_font_size, nwidth=43){
	dirs <- list.dirs(dir.to.placebo, recursive=FALSE)
	match <- grepl(as.character(i), dirs)
	if (sum(match) <= 0){
			return("DID NOT ENTER VALID TRIAL NUMBER")
	}

	summary <- file.path(dir.to.placebo, "placebo.summary.csv")
	s <- read_csv(summary)

	case <- dirs[match]

	if (length(case)>1){
		case <- case[2]
	}

	file <- paste0(case, "/raw.results.rds")
	raw.results <- readRDS(file)

	s.f <- s %>% filter(name==basename(case))

	colnames(raw.results)[which(colnames(raw.results)=="drugA")] <- s.f$drugA_name
	colnames(raw.results)[which(colnames(raw.results)=="drugAB")] <- "Combination"
	colnames(raw.results)[which(colnames(raw.results)=="placebo")] <- "Placebo"
	colnames(raw.results)[which(colnames(raw.results)=="y.model")] <- "Model"

	tib <- pivot_longer(raw.results, -time, names_to="drug", values_to="survival")
	tib$drug <- as.factor(tib$drug)
	orig.levels <- levels(tib$drug)
	comb <- which(orig.levels=="Combination")
	mod <- which(orig.levels=="Model")
	others <- setdiff(1:4, c(comb, mod))

	tib$drug <- factor(tib$drug, levels = c(orig.levels[mod], orig.levels[comb], orig.levels[others[1]], orig.levels[others[2]]), ordered = TRUE)

	mytheme <- theme(
		panel.background = element_rect(fill = "transparent"), # bg of the panel
		panel.grid.major = element_blank(), # get rid of major grid
		panel.grid.minor = element_blank(), # get rid of minor grid
	) + theme(axis.line = element_line(color="black", size = 0.75)) + theme(axis.text.x = element_text(size = axis_text_size),
				axis.text.y = element_text(size = axis_text_size))  + theme(
	axis.title.x = element_text(size = title_font_size),
		axis.title.y = element_text(size = title_font_size)) + theme(plot.title = element_text(size = title_font_size))

	plt <- ggplot(tib, aes(x=time, y=survival, color=drug))  + geom_line(size=1.3) + mytheme + xlab("Time (months)") + ylab("PFS") + ggtitle(stringr::str_wrap(s.f$fullname, width=nwidth))

	legend <- theme(legend.justification = 'left', legend.position="bottom", legend.title = element_blank(), legend.key = element_rect(colour = "transparent", fill = "white"), legend.text=element_text(size=legend_font_size), )

	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	#x label tick marks should be integers not floats! SO USE THE BELOW CODE APPROPRIATELY
	#breaks.x <-
	#scale_x_continuous(breaks=c(1,3,7,10))
	plt2 <- plt + legend + scale_color_manual(values=cbbPalette) + guides(colour = guide_legend(nrow = 2))
	if (off){
			plt2 <- plt + legend + scale_color_manual(values=cbbPalette)
	}
	return(plt2)

}



f <- dotplot("../results.clinical", axis_text_size, title_font_size, legend_font_size-4)

store.plts <- readRDS("MODEL_panels_clinical_trials.rds")

model.trials.to.show <- c(2, 4, 16, 12, 17)

a <- store.plts[[model.trials.to.show[1]]]
b <- store.plts[[model.trials.to.show[2]]]
c <- store.plts[[model.trials.to.show[3]]]
d <- store.plts[[model.trials.to.show[4]]]
e <- store.plts[[model.trials.to.show[5]]]

fig2 <- plot_grid(f, a, b, c, d, e, label_size=label_font_size, labels="AUTO", ncol=2)

ggsave(fig2, device=cairo_pdf, filename=paste0("figure2.pdf"), width=12, height=18, units="in", dpi=320)

### Supplemental figures for placebo

p.list <- c(2,4,8,14)

store.p <- vector("list", 4L)

for (y in 1:length(p.list)){
	x <- p.list[y]
	temp <- placebo.plot(x, "../placebo.results", FALSE, axis_text_size, legend_font_size, title_font_size)
	store.p[[y]] <- temp
}

a <- store.p[[1]]
b <- store.p[[2]]
c <- store.p[[3]]
d <- store.p[[4]]

placebo.others <- plot_grid(a,b,c,d, labels="AUTO", label_size=label_font_size)

ggsave(placebo.others, device=cairo_pdf, filename="supplement/placebo.pdf", width=13, height=13, units="in", dpi=320)



### Supplemental Figures SHOW OTHER MODEL TRIALS DATA
store.plts <- readRDS("MODEL_panels_clinical_trials.rds")

n.others <- setdiff(1:length(store.plts), c(model.trials.to.show, trial.to.follow ))

fig.idxs <- list(s1=n.others[1:4], s2=n.others[5:8], s3=n.others[9:12])

for (y in 1:length(fig.idxs)){
	x <- fig.idxs[[y]]

	a <- store.plts[[x[1]]]
	b <- store.plts[[x[2]]]
	c <- store.plts[[x[3]]]
	d <- store.plts[[x[4]]]

	
	model.sup <- plot_grid(a,b,c,d,label_size=label_font_size, labels="AUTO", ncol=2)
	ggsave(model.sup, device=cairo_pdf, filename=paste0("supplement/model_results_", names(fig.idxs)[y], ".pdf"), width=12.65, height=12.65, units="in", dpi=320)
	
}

print("executed completely")






