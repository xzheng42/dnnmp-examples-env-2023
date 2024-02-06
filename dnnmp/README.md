## Discrete nearest-neighbor mixture process (DNNMP)

This is the R package **dnnmp** (currently developer's version) for the paper:

Xiaotian Zheng. Athanasios Kottas. Bruno Sansó. 
Zheng, X., Kottas, A., & Sansó, B. (2023). Bayesian geostatistical modeling for discrete-valued processes. Environmetrics, 34(7), e2805. https://doi.org/10.1002/env.2805

You can install the **dnnmp** package with **devtools**
```
devtools::install_github("xzheng42/dnnmp-examples-env-2023", subdir = "dnnmp")
```

Main functions of the package are `dnnmp` and `predict.dnnmp`:

- `dnnmp` fits an NNMP model via Markov chain Monte Carlo (MCMC).
- `predict.dnnmp` (or simply `predict`) generates posterior predictive samples for a set of new locations given a dnnmp object.

Detailed guidelines for using the functions are referred to their help pages in R. 
R scripts to reproduce results in the paper are available in [*data-examples/*](https://github.com/xzheng42/dnnmp-examples-env-2023/tree/main/data-examples)
with instructions available in [*dnnmp-examples-env-2003*](https://github.com/xzheng42/dnnmp-examples-env-2023/).

Notes: The current version was tested on macOS 10.15.7 under R version 4.2.2 and on Fedora Linux 38 under R version 4.3.2.

We are currently working on combining the two R packages [**nnmp**](https://github.com/xzheng42/nnmp-examples-ba-2023/tree/main/nnmp) and **dnnmp** into a single one.