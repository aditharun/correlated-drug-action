library(tidyverse)
library(readxl)

#########################################
ref <- readRDS("ref.rds")

results.dir <- "../results.clinical"
placebo.dir <- "../placebo.results"
fig.dir <- "../figures"
supp.dir <- file.path(fig.dir, "supplement")

ci <- read_csv(file.path(results.dir, "ci.csv"))
res <- read_csv(file.path(results.dir, "results.csv"))


dir.create(fig.dir, showWarnings = FALSE)
dir.create(supp.dir, showWarnings = FALSE)

########################################
ci$basename <- basename(ci$name)

df <- left_join(res, ci, by=c("name"="basename")) %>% select(-name.y)

summary <- left_join(df, ref, by=c("name"="name"))

saveRDS(summary, file.path(results.dir, "summary.rds"))

supp.table <- summary %>% select(-c(A_name, B_name)) %>% rename(lower.95.ci=lower, upper.95.ci=upper, trial=fullname, p.value=p.model, sp.cor.estimate=rho.model)

supp.table$`Trial ID` <- as.numeric(unlist(lapply(supp.table$name,function(x) unlist(str_split(x, "-"))[1])))

supp.table <- supp.table %>% select(`Trial ID`, trial, lower.95.ci, sp.cor.estimate, upper.95.ci, p.value) %>% rename(`Combination`=trial)
supp.table <- supp.table %>% arrange(`Trial ID`)
write_csv(supp.table, file.path(supp.dir, "Supplemental File 1.csv"))

summary2 <- file.path(placebo.dir, "placebo.summary.csv")
s <- read_csv(summary2)

s$drugA_name <- c("Gemcitabine, Cisplatin", "Gemcitabine, Carboplatin", "Dabrafenib", "Gemcitabine")
s$fullname <- c("Estimate of Placebo Survival in Advanced Non-Small Cell Lung Cancer", "Estimate of Placebo Survival in Recurrent Ovarian Cancer", "Estimate of Placebo Survival in Metastatic BRAF-Mutant Cutaneous Melanoma", "Estimate of Placebo Survival in Advanced Pancreatic Cancer")

write_csv(s, summary2)

placeboresults <- left_join(s, summary, by=c("name"="name")) %>% rename(parent_trial=fullname.y, hill.param.k=k, hill.param.n=n, reference_drug=drugA_name) %>% select(-c(name, fullname.x, rho.model, p.model, lower, upper, A_name, B_name)) %>% select(parent_trial, reference_drug, hill.param.k, hill.param.n, p.value)

write_csv(placeboresults, file.path(supp.dir, "Supplemental File 3.csv"))