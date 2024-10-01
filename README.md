# Panel with Interactive Fixed Effects

This package implements the estimation and inference procedure for panel data with interactive fixed effects. This package provides two different method for estimation: one is the commonly used linear panel data model estimation procedure, and another one is the bias-aware estimation procedure that allows weak factors.

### Installation

The latest development version can be installed from the Github source using `devtools` package in R:

```R
devtools::install_github("chenweihsiang/PanelIFE")
```

### Examples

First, the following code generates the sample data:

```R
library(PanelIFE)
# Generate sample data with number of factors equals 2
dt <- sample_data(N = 100, T = 20, R = 2, kappa = c(0.5, 0.3))
```

The following example considers the method of least squares estimation of linear panel data models with interactive fixed effects (Bai, 2009; Moon and Weidner, 2017):

```R
# Estimate the linear panel data model with interactive fixed effects
res <- ls_factor(Y = dt$Y, X = dt$X, R = dt$R, report = "silent",
                 precision_beta = 10^(-8), method = "m1",
                 start = c(0), repMIN = 3, repMAX = 10, M1 = 2, M2 = 2)
# Summarize and print the results
sum_res <- summary(res)
print(sum_res)
```

The following example considers the method of robust estimation and inference in panels with interactive fixed effects which can accommodate weak factors (Armstrong, Weidner, and Zeleneev, 2023):

```R
# Estimate the linear panel data model with interactive fixed effects
res <- honest_weak_factors(Y = dt$Y, X = dt$X, R = dt$R,
                           Gamma_LS = NULL, alpha = 0.05)
# Summarize and print the results
sum_res <- summary(res)
print(sum_res)
```

### References

- Armstrong, T. B., Weidner, M., & Zeleneev, A. (2022). Robust estimation and inference in panels with interactive fixed effects. *arXiv preprint arXiv:2210.06639*. [[arxiv](https://doi.org/10.48550/arXiv.2210.06639)]
- Bai, J. (2009). Panel data models with interactive fixed effects. *Econometrica*, *77*(4), 1229-1279. [[paper](https://doi.org/10.3982/ECTA6135)]
- Moon, H. R., & Weidner, M. (2015). Linear regression for panel with unknown number of factors as interactive fixed effects. *Econometrica*, *83*(4), 1543-1579. [[paper](https://doi.org/10.3982/ECTA9382)]
- Moon, H. R., & Weidner, M. (2017). Dynamic linear panel regression models with interactive fixed effects. *Econometric Theory*, *33*(1), 158-195. [[paper](https://doi.org/10.1017/S0266466615000328)]
