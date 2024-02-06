## DNNMP: Discrete nearest-neighbor mixture process

This repository contains the R package [**dnnmp**](https://github.com/xzheng42/dnnmp-examples-env-2023/tree/main/dnnmp) (currently developer's version) 
and R scripts to reproduce the numerical results in

Xiaotian Zheng. Athanasios Kottas. Bruno Sansó. 
Zheng, X., Kottas, A., & Sansó, B. (2023). Bayesian geostatistical modeling for discrete-valued processes. Environmetrics, 34(7), e2805. https://doi.org/10.1002/env.2805

### Installing and using the **dnnmp** package

You can install the package with **devtools**
```
devtools::install_github("xzheng42/dnnmp-examples-env-2023", subdir = "dnnmp")
library(dnnmp)
```

Main functions of the package are `dnnmp` and `predict.dnnmp`:

- `dnnmp` fits a DNNMP model via Markov chain Monte Carlo (MCMC).
- `predict.dnnmp` (or simply `predict`) generates posterior predictive samples for a set of new locations given a dnnmp object.

Notes: The current version was tested on macOS 10.15.7 under R version 4.2.2 and on Fedora Linux 38 under R version 4.3.2.

### Workflow to reproduce numerical results

R scripts to reproduce results in the paper are available in 
[*data-examples/*](https://github.com/xzheng42/dnnmp-examples-env-2023/tree/main/data-examples),
and [*data/*](https://github.com/xzheng42/dnnmp-examples-env-2023/tree/main/data) contains the survey routes (names, numbers, location information, status, etc.)
available in the [North American Breeding Bird Survey (BBS) Dataset](https://www.sciencebase.gov/catalog/item/52b1dfa8e4b0d9b325230cd9)
and the northern cardinal (*Cardinalis cardinalis*) count data.

- Run all simulation experiments: `run_all_sim_rscripts.R`.
- Run all BBS data examples: `run_all_bbs_rscripts.R`.

- Simulation experiments in Section 5.1 and Section D.1: 
  `sim1_mcmc.R` and `sim1_results` (first simulation experiment); 
  `sim2_mcmc.R`, `sim2_model_comp`, `sim2_results` (second simulation experiment);
  
- Data analyses in Section 5.3: `bbs_data_analysis.R`.
- Data analyses in Section D.2: 
  `bbs_gaus_sa_1_mcmc.R`, `bbs_gaus_sa_2_mcmc.R`, `bbs_gaus_sa_1_results.R`, and `bbs_gaus_sa_2_results.R` (Section D.2.1); 
  `bbs_comp_1_mcmc.R` and `bbs_comp_1_results.R` (Section D.2.2); 
  `bbs_comp_2_mcmc.R` and `bbs_comp_2_results.R` (Section D.2.3);

- Randomized quantile residual analysis for model checking:
  `sim_model_checking.R` and `bbs_model_checking.R`.