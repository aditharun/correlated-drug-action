##Doses Figures
library(tidyverse)
library(readxl)
library(ggpubr)
library(cowplot)
source("figure_functions_dose.R")


###############################################
#title_size <- 18
#axis_text_size <- 15
#label_size <- 25
#legend_size <- 14
#facet_size <- 16
#pvalue_cex <- 6.4

facet_size = 16

axis_text_size = 16
title_size = 21
legend_size = 18
pvalue_cex = 7
label_size = 27

size_geom_text = 5

dose.dir <- "../results.cell.line/doses.results.Bliss.Loewe"
sp.dose.dir <- "../results.cell.line/sp.doses.results.Bliss.Loewe"
#############################################

#format dose.results.csv into supplement form 
organize.doses <- function(x){
	r <- readRDS(x)

	guess.0 <- r$guess0
	drugB.name <- r$drugB
	drugA.name <- r$drugA
	guess.final <- r$guessfinal
	true <- r$true
	concentrations <- r$conc
	id.concs <- r$id.concs
	guess.init <- r$guessinit

	true.values <- true %>% select(-row.id)
	unwrap.true <- unlist(unname(as.list(true.values)))
	unwrap.guess.init <- unlist(unname(as.list(guess.init)))
	n <- length(unwrap.true)
	unwrap.lm <- lm(unwrap.true ~ unwrap.guess.init)
	outlier <- outlier.finder(unwrap.lm , n)

	unwrap.guess0 <- unlist(unname(as.list(guess.0)))
	unwrap.guess.final <- unlist(unname(as.list(guess.final)))

	if (!identical(outlier, integer(0))){
	    unwrap.guess0 <- unwrap.guess0[-outlier]
	    unwrap.guess.final <- unwrap.guess.final[-outlier]
	    unwrap.true <- unwrap.true[-outlier]
	    n <- length(unwrap.true)
	}

	pvalue <- ostt(unwrap.guess.final, unwrap.true)$p

	sp <- grepl("sp", x)

	if (sp){
		list(pvalue=pvalue, name=paste0(r$hour, " Hr ", drugA.name, " and ", drugB.name))
	} else{
		list(pvalue=pvalue, name=paste0(drugA.name, " and ", drugB.name))
	}
}




barplot.doses <- function(dose.dir, sp.dose.dir, axis_text_size, title_size, size_geom_text){

	results <- read_csv(list.files(dose.dir, full.names=TRUE, pattern="*.csv")) %>% drop_na()
	results2 <- read_csv(list.files(sp.dose.dir, full.names=TRUE, pattern="*.csv")) %>% drop_na() %>% select(-c(hours,name))

	results <- rbind(results, results2)

	files <- list.files(dose.dir, recursive=TRUE, full.names=TRUE, pattern="*.rds")
	files2 <- list.files(sp.dose.dir, recursive=TRUE, full.names=TRUE, pattern="*.rds")

	files3 <- c(files, files2)

	names <- do.call(rbind, lapply(files3, function(x) organize.doses(x) %>% unlist() %>% unname())) %>% as_tibble()
	colnames(names) <- c("pvalue", "name")

	names$pvalue <- as.numeric(names$pvalue)

	final <- results %>% left_join(names, by=c("t.final"="pvalue")) %>% mutate(trial_id=1:n()) %>% select(trial_id, name, final.rho, init.rho, outliers, t.final, eob.t)

	write.csv(final, "supplement/doses.results.csv")

	mylevels.x <- unique(final$trial_id)
	final$trial_id <- factor(final$trial_id, mylevels.x[order(as.numeric(mylevels.x))])

	k <- ggplot(final, aes(x=trial_id, y=final.rho+0.001)) + geom_bar(stat="identity") 

	goodfit <- final %>% filter(t.final > 0.05) %>% mutate(p=ifelse(final.rho < 0, final.rho-0.007, final.rho+0.007))
	label.df <- data.frame(trial_id = goodfit %>% pull(trial_id),
	                       final.rho = goodfit %>% pull(p))

	eobfit <- final %>% filter(eob.t > 0.05) %>% mutate(p=ifelse(final.rho < 0, final.rho-0.025, final.rho+0.02))
	label.df2 <- data.frame(trial_id = eobfit %>% pull(trial_id),
	                       final.rho = eobfit %>% pull(p))

	k <- k + geom_text(data = label.df, label = "\u2021", size=(size_geom_text+1)) + geom_text(data = label.df2, label = "**", size=(size_geom_text+1))

		mytheme <- theme(
			panel.background = element_rect(fill = "transparent"), # bg of the panel
			panel.grid.major = element_blank(), # get rid of major grid
			panel.grid.minor = element_blank(), # get rid of minor grid
		) + theme(axis.line = element_line(color="black", size = 0.75)) + theme(axis.text.x = element_text(size = axis_text_size),
					axis.text.y = element_text(size = axis_text_size))  + theme(
		axis.title.x = element_text(size = title_size),
			axis.title.y = element_text(size = title_size)) + theme(plot.title = element_text(size = title_size))

	k <- k + mytheme + xlab("Trial ID") + ylab("Estimate for \u03c1") 
	k <- k + geom_hline(yintercept=0)
	k <- k + annotate(geom = 'text', label = "   ** = EOB condition is valid, \u2021 = Good fit to observed combination" , x = -Inf, y = -0.17, hjust = 0, vjust = 1, size=(size_geom_text+1.5))

	k
}



#-------------------------------------------------------------------------------------------------#

###Include special dose results too

#location for results
loc <- dose.dir
loc2 <- sp.dose.dir

pvalues.scatter <- function(loc, loc2, legend_size, title_size, axis_text_size){
	library(readxl)
	library(ggExtra)

	doses.results <- read_csv(file.path(loc, "doses.results.csv"))
	sp.doses.results <- read_csv(file.path(loc2, "sp.doses.results.csv"))

	threshold <- -log(0.05, base=10)
	doses.results <- doses.results %>% mutate(transformed.fit.p=-log(t.final, base=10), transformed.eob.p=-log(eob.t, base=10))
	sp.doses.results <- sp.doses.results %>% mutate(transformed.fit.p=-log(t.final, base=10), transformed.eob.p=-log(eob.t, base=10))

	doses.results <- doses.results %>% drop_na()
	sp.doses.results <- sp.doses.results %>% drop_na()

	doses.results <- rbind(doses.results, sp.doses.results %>% select(-c(name, hours)))
	doses.results$guide <- "none"
	for (k in 1:dim(doses.results)[1]){
		bool1 <- doses.results$transformed.fit.p[k] <= threshold
		bool2 <- doses.results$transformed.eob.p[k] <= threshold


		if (bool1 & bool2){
			item <- "Fit and EOB is good"
		}
		if (bool1 & !bool2){
			item <- "Fit is good, EOB is not"
		}
		if (!bool1 & bool2){
			item <- "EOB is good, Fit is not good"
		}
		if (!bool1 & !bool2){
			item <- "EOB and Fit both not good"
		}

		doses.results$guide[k] <- item
	}

	#sum(doses.results$transformed.fit.p <= threshold) -- 20 / 26 is good for Fit
	#sum(doses.results$transformed.eob.p <= threshold) -- 17 / 26 is good for EOB
	#sum(doses.results$guide=="Fit and EOB is good") --- 11 / 26 is good for both


	mytheme <- theme(
		panel.background = element_rect(fill = "transparent"), # bg of the panel
		panel.grid.major = element_blank(), # get rid of major grid
		panel.grid.minor = element_blank(), # get rid of minor grid
	) + theme(axis.line = element_line(color="black", size = 0.75)) + theme(axis.text.x = element_text(size = axis_text_size),
				axis.text.y = element_text(size = axis_text_size))  + theme(
	axis.title.x = element_text(size = title_size),
		axis.title.y = element_text(size = title_size)) + theme(plot.title = element_text(size = title_size))


	legend <- theme(legend.justification = 'left', legend.position="bottom", legend.title = element_blank(), legend.key = element_rect(colour = "transparent", fill = "white"), legend.text=element_text(size=legend_size), ) 
	cbbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	plot <- ggplot(doses.results, aes(x=transformed.fit.p, y=transformed.eob.p, color=guide)) + geom_point(size = 2.4) + geom_hline(yintercept=threshold, linetype="dashed", color="grey", size=1.6) + geom_vline(xintercept=threshold, linetype="dashed", color="grey", size=1.6) + mytheme + xlab("-log10 p-value for GoF") + ylab("-log10 p-value for EOB") + legend + scale_color_manual(values=cbbPalette) + guides(colour = guide_legend(nrow = 2))
	
	#for 2 rows of legend instead of 1
	#plot <- plot +  guides(color = guide_legend(nrow = 2))

	final <- ggMarginal(plot, type="densigram", color="darkgrey", bins=20)

	#ggsave(final, device="pdf", filename="../../figures/scatterp.pdf", width=10, height=10, units="in", dpi=320)

	return(final)
}



#iterate through eaceh file in files
#heatmap, regression and concentration map

files <- list.files(dose.dir, recursive=TRUE, full.names=TRUE, pattern="*.rds")
files2 <- list.files(sp.dose.dir, recursive=TRUE, full.names=TRUE, pattern="*.rds")

files3 <- c(files, files2)


dose.maps.plots <- function(dose.dir, sp.dose.dir, x, item_labels, title_size, axis_text_size, legend_size, label_size, pvalue_cex){



	mytheme <- theme(
		panel.background = element_rect(fill = "transparent"), # bg of the panel
		panel.grid.major = element_blank(), # get rid of major grid
		panel.grid.minor = element_blank(), # get rid of minor grid
	) + theme(axis.line = element_line(color="black", size = 0.75)) + theme(axis.text.x = element_text(size = axis_text_size),
				axis.text.y = element_text(size = axis_text_size))  + theme(
	axis.title.x = element_text(size = title_size),
		axis.title.y = element_text(size = title_size)) + theme(plot.title = element_text(size = title_size+5))
	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	legend <- theme(legend.justification = 'left', legend.position="bottom", legend.title = element_blank(), legend.key = element_rect(colour = "transparent", fill = "white"), legend.text=element_text(size=(legend_size-2)))
	

	r <- readRDS(x)

	guess.0 <- r$guess0

	sp <- grepl("sp", x)

	if (sp){
		drug.hr <- r$hour
	}

	drugB.name <- r$drugB
	drugA.name <- r$drugA
	guess.final <- r$guessfinal
	true <- r$true
	concentrations <- r$conc
	id.concs <- r$id.concs
	guess.init <- r$guessinit


	true.values <- true %>% select(-row.id)
	unwrap.true <- unlist(unname(as.list(true.values)))
	unwrap.guess.init <- unlist(unname(as.list(guess.init)))
	n <- length(unwrap.true)
	unwrap.lm <- lm(unwrap.true ~ unwrap.guess.init)
	outlier <- outlier.finder(unwrap.lm , n)

	unwrap.guess0 <- unlist(unname(as.list(guess.0)))
	unwrap.guess.final <- unlist(unname(as.list(guess.final)))

	if (!identical(outlier, integer(0))){
	    unwrap.guess0 <- unwrap.guess0[-outlier]
	    unwrap.guess.final <- unwrap.guess.final[-outlier]
	    unwrap.true <- unwrap.true[-outlier]
	    n <- length(unwrap.true)
	}

	pvalue <- ostt(unwrap.guess.final, unwrap.true)$p

	sp.doses.results <- read_csv(paste0(sp.dose.dir, "/sp.doses.results.csv")) %>% as_tibble() %>% drop_na()
	doses.results <- read_csv(paste0(dose.dir, "/doses.results.csv")) %>% as_tibble() %>% drop_na()

	ps <- c(sp.doses.results$t.final, doses.results$t.final)
	rs <- c(sp.doses.results$final.rho, doses.results$final.rho)

	rho <- rs[which.min(abs(ps-pvalue))]
	rho <- round(rho, 2)

	pround <- round(pvalue, 3)
	if ((pround==0)|(pround==1)){
		pround <- format(pvalue, scientific=TRUE, digits=3)
	}
	#add p value to regression bottom corner like fig 2 b-e

	data <- data.frame(gf=unwrap.guess.final, tf=unwrap.true)

	combination.name <- paste0(drugA.name, " and ", drugB.name)

	title.combination <- ggdraw() + draw_label(paste0(combination.name," Combination"), size = (title_size+6), fontface='bold') #to left align use x = 0, hjust = 0


	if (sp){
	title <- paste0("Time: ", drug.hr, " Hours")
	}
	#regression line
	p <- ggplot(data, aes(x=gf, y=tf)) + geom_point(size=2.4) + mytheme + xlab(paste0("Estimate for ", combination.name, " Viability")) + ylab(paste0("Observed ", combination.name," Viability")) + geom_abline(slope=1, intercept=0, color=cbbPalette[8], linetype="dashed", size=1.6)
	p <- p + annotate(geom = 'text', label = paste0('   GoF p-value = ', pround, '; \u03C1 = ', rho) , x = -Inf, y = 5, hjust = 0, vjust = 1, size=(pvalue_cex+1))

	error <- guess.final - true.values
	heatmap <- true
	coln <- dim(heatmap)[2]
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

	    coloritems <- c("non-outlier" = "transparent", "outlier" = "black")
	    #make things numerically ordered on both axes - could be reordering the data wrong!
	    mylevels.x <- unique(long$drugB)
	    mylevels.y <- unique(long$drugA)
	    long$drugB2 <- factor(long$drugB, mylevels.x[order(as.numeric(mylevels.x))])
	    long$drugA2 <- factor(long$drugA, mylevels.y[order(as.numeric(mylevels.y))])

	    b.title <- paste0(drugB.name, " Concentration (\u00B5M)") 
	    a.title <- paste0(drugA.name, " Concentration (\u00B5M)")

	    #heatmap
	    plot <- ggplot(long, aes(drugB2, drugA2)) + scale_fill_gradient2(name="Excess over CDA", low="darkblue", high="red", guide="colorbar") +  geom_tile(aes(fill = error, color=outlier_status, width=0.95, height=0.95), size=2) + theme(
	        panel.grid.major = element_blank(),
	        panel.grid.minor = element_blank(),
	        panel.background = element_rect(fill = "transparent",colour = NA),
	        plot.background = element_rect(fill = "transparent",colour = NA)
	      ) + scale_color_manual(name="Outlier Status", values=coloritems)
	    plot <- plot + xlab(b.title) + ylab(a.title) + guides(color=FALSE) 
		plot <- plot + theme(axis.text.x = element_text(size = axis_text_size),
					axis.text.y = element_text(size = axis_text_size))  + theme(
		axis.title.x = element_text(size = title_size),
			axis.title.y = element_text(size = title_size), legend.title=element_text(size=legend_size), legend.text=element_text(size=(legend_size)))
		plot <- plot + theme(panel.background = element_rect(colour = "black", size=1.7))


		if (sum(is.na(id.concs))==0){

		  conc.map <- true
	      id.map <- true
	      conc.map[,2:coln] <- concentrations
	      id.map[,2:coln] <- id.concs
	      long.conc <- pivot_longer(conc.map, -row.id, names_to="drugB", values_to="eq.conc")
	      long.id <- pivot_longer(id.map, -row.id, names_to="drugB", values_to="id")
	      long.conc <- inner_join(long.conc, long.id)
	      colnames(long.conc)[1] <- "drugA"

	      long.conc$eq.conc <- as.character(round(long.conc$eq.conc, 2))
	      long.conc$eq.conc <- paste0(long.conc$eq.conc, "\n\u00B5M")
	      mylevels.x <- unique(long.conc$drugB)
	      mylevels.y <- unique(long.conc$drugA)
	      long.conc$drugB2 <- factor(long.conc$drugB, mylevels.x[order(as.numeric(mylevels.x))])
	      long.conc$drugA2 <- factor(long.conc$drugA, mylevels.y[order(as.numeric(mylevels.y))])

	      plot2 <- ggplot(long.conc, aes(drugB2, drugA2)) +  geom_tile(aes(fill = id, width=0.95, height=0.95)) + theme(
	          panel.grid.major = element_blank(),
	          panel.grid.minor = element_blank(),
	          panel.background = element_rect(fill = "transparent",colour = NA),
	          plot.background = element_rect(fill = "transparent",colour = NA)
	        )  + geom_text(aes(label=eq.conc), size=(size_geom_text+1))


	      plot2 <- plot2 + xlab(b.title) + ylab(a.title)  + scale_fill_manual(name = "", labels = c(paste0(drugA.name, " Equivalent Dose"), paste0(drugB.name, " Equivalent Dose")), values=c(cbbPalette[2], cbbPalette[3]))
		  plot2 <- plot2 + guides(fill = guide_legend(nrow = 2)) + legend
		  plot2 <- plot2 + theme(axis.text.x = element_text(size = axis_text_size),
					axis.text.y = element_text(size = axis_text_size))  + theme(
		axis.title.x = element_text(size = title_size),
			axis.title.y = element_text(size = title_size), legend.title=element_text(size=legend_size))

			
		is.conc <- "YES"
		} else{
			is.conc <- "NO"
		}

		if (is.conc=="YES"){

			if (sp){
				p <- p + ggtitle(title)
				plot <- plot + ggtitle(title)
				plot2 <- plot2 + ggtitle(title)
			}
			c("1", "2", "3")

			fp <- plot_grid(plot, p, plot2, labels=item_labels,label_size=label_size, nrow=1)
		} else{

			if (sp){
				p <- p + ggtitle(title)
				plot <- plot + ggtitle(title)
			}

			fp <- plot_grid(plot, p, labels=item_labels, label_size=label_size, nrow=1)
		}

		plot_grid(title.combination, fp, ncol = 1, rel_heights = c(0.1, 1))
}



a <- pvalues.scatter(loc, loc2, legend_size, title_size, axis_text_size)
d <- barplot.doses(dose.dir, sp.dose.dir, axis_text_size-2, title_size, size_geom_text)

ex1 <- files[grepl("combo3/results5", files3)]
ex2 <- files[grepl("combo2/results1", files3)]


b <- dose.maps.plots(dose.dir, sp.dose.dir, ex1, c("C", "D", "E"),  title_size, axis_text_size, legend_size, label_size, pvalue_cex)
c <-  dose.maps.plots(dose.dir, sp.dose.dir, ex2, c("F", "G", "H"),  title_size, axis_text_size, legend_size, label_size, pvalue_cex)

figure_dose <- plot_grid(plot_grid(d, a, labels="AUTO", label_size=label_size), b, c, labels=c("", "", ""), label_size=label_size, ncol=1)
ggsave(figure_dose, device=cairo_pdf, filename="figure3.pdf", width=20, height=30, units="in", dpi=320)


remaining <- setdiff(files3, c(ex1, ex2))


assign_priority <- function(x, dose.dir, sp.dose.dir){


	r <- readRDS(x)
	sp <- grepl("sp", x)
	guess.0 <- r$guess0
	guess.final <- r$guessfinal
	true <- r$true
	concentrations <- r$conc
	id.concs <- r$id.concs
	guess.init <- r$guessinit
	true.values <- true %>% select(-row.id)
	unwrap.true <- unlist(unname(as.list(true.values)))
	unwrap.guess.init <- unlist(unname(as.list(guess.init)))
	n <- length(unwrap.true)
	unwrap.lm <- lm(unwrap.true ~ unwrap.guess.init)
	outlier <- outlier.finder(unwrap.lm , n)
	unwrap.guess0 <- unlist(unname(as.list(guess.0)))
	unwrap.guess.final <- unlist(unname(as.list(guess.final)))
	if (!identical(outlier, integer(0))){
	    unwrap.guess0 <- unwrap.guess0[-outlier]
	    unwrap.guess.final <- unwrap.guess.final[-outlier]
	    unwrap.true <- unwrap.true[-outlier]
	    n <- length(unwrap.true)
	}
	pvalue <- ostt(unwrap.guess.final, unwrap.true)$p
	sp.doses.results <- read_csv(paste0(sp.dose.dir, "/sp.doses.results.csv")) %>% as_tibble() %>% drop_na()
	doses.results <- read_csv(paste0(dose.dir, "/doses.results.csv")) %>% as_tibble() %>% drop_na()
	ps <- c(sp.doses.results$t.final, doses.results$t.final)
	rs <- c(sp.doses.results$final.rho, doses.results$final.rho)
	rho <- rs[which.min(abs(ps-pvalue))]
	rho <- round(rho, 2)

	if (rho < 0){
		n.panels <- 2
	} else{
		n.panels <- 3
	}

	c(n.panels, sp)
}

lookup <- do.call(rbind, lapply(remaining, function(x) assign_priority(x, dose.dir, sp.dose.dir))) 
colnames(lookup) <- c("n.panels", "special.dose")
lookup <- lookup %>% as_tibble() %>% mutate(file=remaining) %>% arrange(special.dose, n.panels)

#I want groups of 3 and there is 24 rows right now so we want 8 groups
num_groups = 8

final.lookup <- lookup %>% group_by((row_number()-1) %/% (n()/num_groups)) %>% nest %>% pull(data)


remaining.plots <- list()
for (x in 1:length(final.lookup)){
		table <- final.lookup[[x]]
		item_labels <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")
	for (j in 1:3){
		used <- item_labels[1:table$n.panels[j]]
		remaining.plots[[(3*(x-1)+j)]] <- dose.maps.plots(dose.dir, sp.dose.dir, table$file[j], used ,title_size, axis_text_size, legend_size, label_size, pvalue_cex )
		item_labels <- item_labels[!item_labels %in% used]
	}
}


for (l in 1:8){
	dr <- plot_grid(remaining.plots[[ (3*l - 2) ]], remaining.plots[[ (3*l - 1)]], remaining.plots[[ (3*l) ]], labels=NULL, nrow=3)
	ggsave(dr, device=cairo_pdf, filename=paste0("supplement/doses_results_", l ,".pdf"), width=24, height=24, units="in", dpi=320)
}







