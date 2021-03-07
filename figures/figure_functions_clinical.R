#helper functions for rendering figures


#monotonicity imposed after the fact!
mono <- function(y){
	n <- length(y)
	for (i in 1:(n-1)){
		if (y[i] < y[i+1]){
			y[i+1] <- y[i]
		}
	}
	return(y)
}

bands <- function(data, p){
		orig.data <- data %>% mutate(better=pmax(drugA,drugB),worse=pmin(drugA,drugB)) %>% select(better, worse)

		if (p >= 0){
			y <- orig.data$better + orig.data$worse*(1-orig.data$better)*(1-p) #rho >= 0
		} else{
		  tstar <- data %>% mutate(neg=abs(drugA+drugB-1)) %>% slice(which.min(neg)) %>% pull(time)
			y <- data %>% mutate(comb=ifelse(time <= tstar, ((drugA + drugB - (drugA*drugB)) * (1+p)) - p, drugA + drugB - (drugA*drugB) * (1 + p))) %>% pull(comb) #rho is negative
		}

		return(y)
}

stitch.ta.tb <- function(data, y, idx, rhos){
		vec <- y[[idx]]
		if (rhos[[idx]] < 0){
			vec <- rev(vec)
		}
		i <- (data.frame(tA=data$drugA, tB=vec))
		corr <- cor(i$tA, i$tB, method="spearman")
		i$id <- rhos[which.min(abs(corr-rhos))]

		return(i)

	}