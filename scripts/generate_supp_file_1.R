setwd("~/Documents/correlated-drug-action/raw-data/clinical_trials/")

library(readxl)
library(writexl)
library(tidyverse)

#non-placebo files
files <- list.files(recursive = TRUE, pattern="*.csv")
f <- files[!grepl("PL",files)]

n.trials <- length(f) / 3

linear.approx <- function(data, interval){
  res <- approx(data$time, data$pfs, seq(0, max(data$time), interval))
  res$y[which(is.na(res$y))] <- 100
  return(list(x=res$x, y=res$y))
}


interval <- 0.01

store <- list()

for (i in 1:n.trials){

	tmp <- list()
	names <- list()
	
	for (k in 0:2){

		ind <- 3*i - k
		data <- read_csv(f[ind], col_names = FALSE)
		colnames(data) <- c("time", "pfs")

		data <- linear.approx(data, interval)

		tmp[[(k+1)]] <- data
		names[[(k+1)]] <- str_sub(basename(f[ind]), 1, -5)

	}

	val <- min(tail(tmp[[1]]$x, n=1), tail(tmp[[2]]$x, n=1), tail(tmp[[3]]$x, n=1))

  	tmp[[1]]$y <- tmp[[1]]$y[tmp[[1]]$x < val]
  	tmp[[2]]$y <- tmp[[2]]$y[tmp[[2]]$x < val]
  	tmp[[3]]$y <- tmp[[3]]$y[tmp[[3]]$x < val]
  	tmp[[1]]$x <- tmp[[1]]$x[tmp[[1]]$x < val]
  	tib <- as_tibble(data.frame(A.x=tmp[[1]]$x, A.y=tmp[[1]]$y, B.y=tmp[[2]]$y, AB.y=tmp[[3]]$y))

  	tib <- tib %>% mutate_at(vars(-A.x), ~./100)

  	colnames(tib) <- c("time", unlist(names))

  	store[[i]] <- tib

}

trials <- list.dirs(full.names = FALSE)[-1]

trial.num <- lapply(str_split( trials, "-"), function(x) x[1]) %>% unlist() %>% as.numeric()

correctorder <- read_excel("../../figures/supplement/Supplemental Table 1.xlsx")


names(store) <- lapply(trial.num, function(j) correctorder %>% filter(`Trial ID` == j) %>% pull(`Combination Name`) %>% paste0(j, "-" ,.)) %>% unlist()

store <- c(store[1], store[11:18], store[2:10])

writexl::write_xlsx(store, col_names=TRUE, "../../figures/supplement/Supplemental File 1.xlsx")
