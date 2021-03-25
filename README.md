# sPLSPM R package
[Multiset sparse partial least squares path modeling for high dimensional data analysis](https://doi.org/10.1186/s12859-019-3286-3)

This R application is based on [plspm](https://cran.r-project.org/web/packages/plspm/index.html)

## Installation
Install the development version of [sPLSPM](https://github.com/acsala/sPLSPM)
```r
# install "devtools"
install.packages("devtools") 
library(devtools)

# install "sPLSPM"
install_github("acsala/sPLSPM")
```
## Minimum example

library(sPLSPM)

```r
Generated_data <- sPLSPM::generate_data(number_of_ksi = 1,
                                        number_of_patients = 150,
                                        number_of_Xs_associated_with_ksis = c(15),
                                        number_of_not_associated_Xs = 100,
                                        mean_of_the_regression_weights_of_the_associated_Xs = c(0.9),
                                        sd_of_the_regression_weights_of_the_associated_Xs = c(0.05),
                                        Xnoise_min = -0.2, Xnoise_max = 0.2,
                                        number_of_Ys_associated_with_ksis = c(15),
                                        number_of_not_associated_Ys = 85,
                                        mean_of_the_regression_weights_of_the_associated_Ys = c(0.9),
                                        sd_of_the_regression_weights_of_the_associated_Ys = c(0.05),
                                        Ynoise_min = -0.2, Ynoise_max = 0.2)
X <- Generated_data$X
Y <- Generated_data$Y
data_info1 <- Generated_data$data_info
Data <- cbind(X,Y)
dim(X)
dim(Y)

dim(Data)
## 2 DATASETS connectivity matrix
EXPL_X = c(0,0)
RESP_Y = c(1,0)
path_matrix = rbind(EXPL_X, RESP_Y)
## blocks of outer model
blocks = list(1:dim(X)[2], dim(X)[2]+1:dim(Y)[2])
dim(Data[,blocks[[1]]])
dim(Data[,blocks[[2]]])
## define vector of reflective modes
modes = c("B","A")

## fit the model
time_data <- system.time(
  s_satpls <- splspm(Data, path_matrix, blocks, modes, scheme="path",
                     scaled=T, penalization = "enet", nonzero = c(40,20),
                     cross_validate = T, lambda = 1, maxiter = 100)
)

s_satpls$outer_model
s_satpls$model$iter
s_satpls$nonzero
s_satpls$lambda
```
