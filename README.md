# rfrontier
## Stata package for robust stochastic frontier analysis
This facilitates estimation of stochastic frontier models with alternative distributional assumptions, including models with in which the noise terms follows a Student's t, Cauchy, or logistic distribution. The package can then be used to generate efficiency predictions and related postestimation outputs. The main purposes of this package are to enable easy implementation of some of the stochastic frontier specifications explored in publications I have co-authored, and to ease replication of some of the results reported in these publications.

## Installation
In order to install the command, enter the following commands into Stata (the first is unnecessary if you already have the `github` command installed):
```stata
net install github, from("https://haghish.github.io/github/")
github install AlexStead/rfrontier
```

## Getting started
For help with the command's syntax and options, enter (after installation):
```stata
help rfrontier
help rfrontier_postestimation
```

## Linked publications

- Wheat P, Stead AD, Greene WH. 2019. Robust stochastic frontier analysis: a Student’s t-half normal model with application to highway maintenance costs in England. Journal of Productivity Analysis. 51(1), pp. 21-38, https://doi.org/10.1007/s11123-018-0541-y
- Stead AD, Wheat P, Greene WH. 2018. Estimating efficiency in the presence of extreme outliers: A logistic-half normal stochastic frontier model with application to highway maintenance costs in England. In: Greene WH; Khalaf L; Makdissi P; Sickles RC; Veall MR; Voia M-C (eds.) Productivity and Inequality. Springer Proceedings in Business and Economics. Springer International Publishing, pp. 1-19, https://doi.org/10.1007/978-3-319-68678-3_1
