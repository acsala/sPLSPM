% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splspm.r
\name{splspm}
\alias{splspm}
\title{PLS-PM: Partial Least Squares Path Modeling}
\usage{
splspm(Data, path_matrix, blocks, modes = NULL, scaling = NULL,
  scheme = "path", scaled = TRUE, tol = 1e-06, maxiter = 100,
  plscomp = NULL, boot.val = FALSE, br = NULL, dataset = TRUE,
  penalization = "none", nonzero = 0, lambda = 0, alpha = 0.5,
  cross_validate = FALSE, parallel_cv = FALSE, nr_subsets = 10,
  warning_non_convergence = TRUE)
}
\arguments{
\item{Data}{matrix or data frame containing the manifest variables.}

\item{path_matrix}{A square (lower triangular) boolean matrix representing
the inner model (i.e. the path relationships between latent variables).}

\item{blocks}{list of vectors with column indices or column names
from \code{Data} indicating the sets of manifest variables forming
each block (i.e. which manifest variables correspond to each block).}

\item{modes}{character vector indicating the type of measurement for each
block. Possible values are: \code{"A", "B", "newA", "PLScore", "PLScow"}.
The length of \code{modes} must be equal to the length of \code{blocks}.}

\item{scaling}{optional argument for runing the non-metric approach;
it is a list of string vectors indicating the type of
measurement scale for each manifest variable specified in \code{blocks}.
\code{scaling} must be specified when working with non-metric variables.
Possible values: \code{"num"} (linear transformation,
suitable for numerical variables), \code{"raw"} (no transformation),
\code{"nom"} (non-monotonic transformation, suitable for nominal variables),
and \code{"ord"} (monotonic transformation, suitable for ordinal variables).}

\item{scheme}{string indicating the type of inner weighting
scheme. Possible values are \code{"centroid"}, \code{"factorial"}, or
\code{"path"}.}

\item{scaled}{whether manifest variables should be standardized.
Only used when \code{scaling = NULL}. When (\code{TRUE}, data is
scaled to standardized values (mean=0 and variance=1).
The variance is calculated dividing by \code{N} instead of \code{N-1}).}

\item{tol}{decimal value indicating the tolerance criterion for the
iterations (\code{tol=0.000001}). Can be specified between 0 and 0.001.}

\item{maxiter}{integer indicating the maximum number of iterations
(\code{maxiter=100} by default). The minimum value of \code{maxiter} is 100.}

\item{plscomp}{optional vector indicating the number of PLS components
(for each block) to be used when handling non-metric data
(only used if \code{scaling} is provided)}

\item{boot.val}{whether bootstrap validation should be performed.
(\code{FALSE} by default).}

\item{br}{number bootstrap resamples. Used only
when \code{boot.val=TRUE}. When \code{boot.val=TRUE}, the default number of
re-samples is 100.}

\item{dataset}{whether the data matrix used in the computations should be
retrieved (\code{TRUE} by default).}
}
\value{
An object of class \code{"plspm"}.

\item{outer_model}{Results of the outer model. Includes:
outer weights, standardized loadings, communalities, and redundancies}

\item{inner_model}{Results of the inner (structural) model.
Includes: path coeffs and R-squared for each endogenous latent variable}

\item{scores}{Matrix of latent variables used to estimate the inner
model. If \code{scaled=FALSE} then \code{scores} are latent variables
calculated with the original data (non-stardardized).}

\item{path_coefs}{Matrix of path coefficients
(this matrix has a similar form as \code{path_matrix})}

\item{crossloadings}{Correlations between the latent variables
and the manifest variables (also called crossloadings)}

\item{inner_summary}{Summarized results of the inner model.
Includes: type of LV, type of measurement, number of indicators, R-squared,
average communality, average redundancy, and average variance
extracted}

\item{effects}{Path effects of the structural relationships.
Includes: direct, indirect, and total effects}

\item{unidim}{Results for checking the unidimensionality of blocks
(These results are only meaningful for reflective blocks)}

\item{gof}{Goodness-of-Fit index}

\item{data}{Data matrix containing the manifest variables used in the
model. Only available when \code{dataset=TRUE}}

\item{boot}{List of bootstrapping results; only available
when argument \code{boot.val=TRUE}}
}
\description{
Estimate path models with latent variables by partial least squares approach
(for both metric and non-metric data)
}
\details{
The function \code{plspm} estimates a path model by partial least squares
approach providing the full set of results. \cr

The argument \code{path_matrix} is a matrix of zeros and ones that indicates
the structural relationships between latent variables. \code{path_matrix}
must be a lower triangular matrix; it contains a 1 when column \code{j}
affects row \code{i}, 0 otherwise. \cr
}
\examples{
\dontrun{
## typical example of PLS-PM in customer satisfaction analysis
## model with six LVs and reflective indicators

# load dataset satisfaction
data(satisfaction)

# path matrix
IMAG = c(0,0,0,0,0,0)
EXPE = c(1,0,0,0,0,0)
QUAL = c(0,1,0,0,0,0)
VAL = c(0,1,1,0,0,0)
SAT = c(1,1,1,1,0,0)
LOY = c(1,0,0,0,1,0)
sat_path = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)

# plot diagram of path matrix
innerplot(sat_path)

# blocks of outer model
sat_blocks = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)

# vector of modes (reflective indicators)
sat_mod = rep("A", 6)

# apply plspm
satpls = plspm(satisfaction, sat_path, sat_blocks, modes = sat_mod,
   scaled = FALSE)

# plot diagram of the inner model
innerplot(satpls)

# plot loadings
outerplot(satpls, what = "loadings")

# plot outer weights
outerplot(satpls, what = "weights")

#NEWER DOCS

Generated_data <- generate_data(number_of_ksi = 1,
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

# Only 2 DATASETS####
EXPL_X = c(0,0)
RESP_Y = c(1,0)
path_matrix = rbind(EXPL_X, RESP_Y)

# blocks of outer model
blocks = list(1:dim(X)[2], dim(X)[2]+1:dim(Y)[2])

dim(Data[,blocks[[1]]])

dim(Data[,blocks[[2]]])

# define vector of reflective modes

modes = c("B","A")



time_data <- system.time(
s_satpls <- splspm(Data, path_matrix, blocks, modes, scheme="path",
scaled=T, penalization = "enet", nonzero = 40, lambda = 1, maxiter = 100)
)

s_satpls$outer_model
s_satpls$model$iter

s_satpls$nonzero
s_satpls$lambda




# ALL DATASETS####

Generated_data <- generate_data(number_of_ksi = 1,
number_of_patients = 150,
number_of_Xs_associated_with_ksis = c(10),
number_of_not_associated_Xs = 150,
mean_of_the_regression_weights_of_the_associated_Xs = c(0.9),
sd_of_the_regression_weights_of_the_associated_Xs = c(0.05),
Xnoise_min = -0.2, Xnoise_max = 0.2,
number_of_Ys_associated_with_ksis = c(10),
number_of_not_associated_Ys = 100,
mean_of_the_regression_weights_of_the_associated_Ys = c(0.9),
sd_of_the_regression_weights_of_the_associated_Ys = c(0.05),
Ynoise_min = -0.3, Ynoise_max = 0.3)

Z <- Generated_data$X
V <- Generated_data$Y

data_info2 <- Generated_data$data_info

Data <- cbind(X,Y,Z,V)

dim(X)
dim(Y)
dim(Z)
dim(V)

dim(Data)

EXPL_X = c(0,1,0,0)
RESP_Y = c(1,0,1,0)
EXPL_Z = c(0,1,0,0)
RESP_V = c(0,0,1,0)
path_matrix = rbind(EXPL_X, RESP_Y,EXPL_Z,RESP_V)

# blocks of outer model
blocks = list(1:dim(X)[2], dim(X)[2]+1:dim(Y)[2],dim(Y)[2]+1:dim(Z)[2],
dim(Z)[2]+1:dim(V)[2])

dim(Data[,blocks[[1]]])

# define vector of reflective modes

modes = c("B","B","B","A")



time_data <- system.time(
s_satpls <- splspm(Data, path_matrix, blocks, modes, scheme="path",
scaled=T, penalization = "enet", nonzero = c(20,40),
lambda = c(0.1,1), maxiter = 100, cross_validate = T, nr_subsets = 10)
)

s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_X",]
s_satpls$outer_model[s_satpls$outer_model[,2]=="RESP_Y",]
s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_Z",]

s_satpls$nonzero
s_satpls$lambda

}
}
\references{
Tenenhaus M., Esposito Vinzi V., Chatelin Y.M., and Lauro C.
(2005) PLS path modeling. \emph{Computational Statistics & Data Analysis},
\bold{48}, pp. 159-205.

Lohmoller J.-B. (1989) \emph{Latent variables path modeling with partial
least squares.} Heidelberg: Physica-Verlag.

Wold H. (1985) Partial Least Squares. In: Kotz, S., Johnson, N.L. (Eds.),
\emph{Encyclopedia of Statistical Sciences}, Vol. 6. Wiley, New York, pp.
581-591.

Wold H. (1982) Soft modeling: the basic design and some extensions. In: K.G.
Joreskog & H. Wold (Eds.), \emph{Systems under indirect observations:
Causality, structure, prediction}, Part 2, pp. 1-54. Amsterdam: Holland.

Russolillo, G. (2012) Non-Metric Partial Least Squares. \emph{Electronic
Journal of Statistics}, \bold{6}, pp. 1641-1669.
\url{http://projecteuclid.org/euclid.ejs/1348665231}
}
\seealso{
\code{\link{innerplot}}, \code{\link{outerplot}},
}
\author{
Gaston Sanchez, Giorgio Russolillo
}
