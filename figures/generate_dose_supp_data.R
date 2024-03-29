##Doses Figures
library(tidyverse)
library(readxl)
library(ggpubr)
library(cowplot)
library(writexl)
source("figure_functions_dose.R")


#############################################
title_size <- 18
axis_text_size <- 15
label_size <- 25
legend_size <- 14
facet_size <- 16
pvalue_cex <- 6.4


results.dir <- "../results.cell.line"


bliss.loewe.dose.dir <- file.path(results.dir, "doses.results.Bliss.Loewe")
bliss.loewe.sp.dose.dir <- file.path(results.dir, "sp.doses.results.Bliss.Loewe")

#bliss.hsa.dose.dir <- file.path(results.dir, "doses.results.Bliss.HSA")
#bliss.hsa.sp.dose.dir <- file.path(results.dir, "sp.doses.results.Bliss.HSA")

##############################################################


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

	pvalue <- t.test(unwrap.guess.final, unwrap.true, paired=TRUE, alternative="two.sided")$p.value

	sp <- grepl("sp", x)

	if (sp){
		list(pvalue=pvalue, name=paste0(r$hour, " Hr ", drugA.name, " and ", drugB.name))
	} else{
		list(pvalue=pvalue, name=paste0(drugA.name, " and ", drugB.name))
	}
}


format.data <- function(dose.dir, sp.dose.dir, LOEWE=FALSE){
	results <- read_csv(list.files(dose.dir, full.names=TRUE, pattern="*.csv")) %>% drop_na()
	results2 <- read_csv(list.files(sp.dose.dir, full.names=TRUE, pattern="*.csv")) %>% drop_na() %>% select(-c(hours,name))

	results <- rbind(results, results2)

	files <- list.files(dose.dir, recursive=TRUE, full.names=TRUE, pattern="*.rds")
	files2 <- list.files(sp.dose.dir, recursive=TRUE, full.names=TRUE, pattern="*.rds")

	files3 <- c(files, files2)

	files3 <- files3[!grepl("regression_*", files3)]

	names <- do.call(rbind, lapply(files3, function(x) organize.doses(x) %>% unlist() %>% unname())) %>% as_tibble()
	colnames(names) <- c("pvalue", "name")


	results$name <- unlist(lapply(as.numeric(results$gof.p), function(x) names$name[which.min(abs(x - as.numeric(names$pvalue)))]))

	results <- results %>% select(name, init.rho, final.rho, outliers, gof.p, eob.p) %>% rename(n.outliers=outliers, gof.p=gof.p, eob.p=eob.p, pre.outlier.rho=init.rho, post.outlier.rho=final.rho)

	if (LOEWE){
		return (	results %>% arrange(desc(post.outlier.rho)) %>% mutate(trial_id=1:n() ) )
	} else{
		return (	results %>% arrange(desc(post.outlier.rho)) )
	}

}

bl <- format.data(bliss.loewe.dose.dir, bliss.loewe.sp.dose.dir, LOEWE=TRUE)
#bh <- format.data(bliss.hsa.dose.dir, bliss.hsa.sp.dose.dir, LOEWE=FALSE)


write_xlsx(list(`dCDA Results`=bl), file.path("supplement", "Supplemental Table 2.xlsx"))


