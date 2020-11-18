# STEPWISE
An R package for STEPWISE method implementation

# Description
The package carries out STEPWISE method that performs a variable selection and consistent estimation of the model parameters in Generalized Linear Models (GLM) via stepwise regression.

# Dependency
The functions in this package involve parallel computing elements and rely on `parallel`, `doParallel` and `foreach` packages.

# Installation
To install this package, first install `devtools` package and execute the following lines:
```
devtools::install_git("https://github.com/AlexPijyan/STEPWISE")
library("STEPWISE")
```

# Example
Suppose we generate data according to `Example 1` from the file `Simulation` with `N = 400` and `p = 1000`, where the response variable follows `Gaussian` distribution. The following code illustrates the utilization of the `STEPWISE` function:

```
library("STEPWISE")
output = STEPWISE(x = x, y = y, family = "gaussian", eta1 = 1, eta2 = 2.5)
output$finalset2
```
`output$finalset2` contains variables selected by `STEPWISE`. Corresponding model coefficients can be found by extracting `output$coef2`. In addition, `STEPWISE` function provides variables selected after Stage 1 and their corresponding coefficients if one desires to use them in the final model. These elements can be accessed at `output$finalset` and `output$coef1` respectively. 

Stopping criterias `eta1` and `eta2` take non-negative values and can be specified in the function. Default values are `0` and `1` respectively.


# Authors

**Alex Pijyan** and **Hyokyoung G. Hong**
