# Correlated Drug Action as a Baseline Additivity Model for Combination Cancer Therapy in Patient Cohorts and Cell Cultures

Adith S. Arun, Sung-Cheol Kim, Mehmet Eren Ahsen, Gustavo Stolovitzky

March 7, 2021

For questions, comment, or any correspondence, please email adith.3.arun@gmail.com

# Overview 

We present additivity models at the level of patient populations (e.g. clinical trials) and cells in culture (cell lines). Namely, the temporal correlated drug action (CDA) model describes a baseline additivtity model for clinical trials. And, the dose-space CDA model acts on cell lines. The raw data, scripts, and code for generating figures are presented in this repository.  

# Data 

The raw data is presented in the `raw-data` directory. 

The clinical trial data aggregated from many different trials can be found within the `raw-data/clinical_trial` directory. Namely, each combination is stored as a folder here, and within each folder are `.csv` files with progression-free survival (PFS) curve data (time and probability of progression free survival). The guide that describes the clinical trials from which we collected data can be found at `raw-data/guide_to_clinical_trials.xlsx`. 

The MCF7 cell line tested combinations are stored within the `raw-data/cell_line` directory. There are four Excel files, each with a number of tested combinations for which many variables are reported including doses, observed viabilities, and Excess over Bliss values. We used data from a time-series experiment (from this [paper](https://elifesciences.org/articles/52707)) where viability data for MCF7 cells for 2 different combinations were collected at different time points (12, 24, 48 hrs). 

# Usage and Workflow

All the code necessary to perform analyses are located within the `scripts` directory. 

To generate Temporal CDA Model results and figures, run the R script `metascript.clintrial.R` which analyzes clinical trials and imputes placebo performance. To generate Dose-space CDA Model results and figures, run the R script `metascript.doses.R` which analyzes cell line experiments. 

### Temporal CDA Model for Clinical Trials

Run the R script `scripts/clinicaltrial.R` to generate the directory `results.clinical.trial` with estimates for the combination PFS cuurve, bootstrapped 95% confidence intervals, and other summary measures. Within the script, there are changeable parameters including the number of bootstrap samples, and interpolation frequency. 

To impute a plausible placebo PFS curve for collected clinical trial combinations for which the relevant data is provided (Trial ID's 2, 4, 8, 14), run the R script `scripts/placebo.R`. 

Subsequently, after running `clinicaltrial.R`, run the R script `scripts/summary_clintrial.R` to generate summary tables with all relevant output parameters. 

### Dose-space CDA Model for Cell Lines

There are two different formulations of the dose-space CDA Model. 

##### The method we use is the one that incorporates the Bliss and Loewe additivity models: 

Run the R script `scripts/doses.loewe.R` to generate the results directory `results.cell.line/doses.results.Bliss.Loewe`. Similarly, for the time series expeirments, run the R script `scripts/sp.doses.loewe.R` to generate the results directory `results.cell.line/sp.doses.results.Bliss.Loewe`. 

##### The alternative method is the one that incorporates the Bliss and Highest Single Agent (HSA) additivity models:

Run the R script `scripts/doses.loewe.R` to generate the results directory `results.cell.line/doses.results.Bliss.Loewe`. Similarly, for the time series expeirments, run the R script `scripts/sp.doses.hsa.R` to generate the results directory `results.cell.line/sp.doses.results.Bliss.HSA`. 

In our tested combinations, both methods perform very similarly. The Bliss-Loewe is favored because we can apply the principle of dose-equivalence from the theory of Loewe additivity to determine the minimum monotherapy dose necessary to recreate the effect of the combination. 

# Figures

All figures and supplemental figures are shown here. To generate them, simply run the `figures/render_figures_clintrial.R` script for clinical trial and placebo figures; `figures/render_figures_sim.R` for simulation based figures; `figures/render_figures_doses.R` and `figures/generate_dose_supp_data.R` for cell line figures.
