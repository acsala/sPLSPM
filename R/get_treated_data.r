#' @title Apply corresponding treatment to MV
#'
#' @details
#' Internal function. \code{get_treated_data} is called by \code{plspm}.
#'
#' @note
#' In non-metric case, all mvs are starndardized, and those ordinal or nominal
#' are rankified
#'
#' @param MV matrix or data frame from with manifest variables
#' @param specs list with algorithm specifications
#' @return matrix or data frame of treated MVs (scaling)
#' @keywords internal
#' @template internals
#' @export

colSdColMeans <- function(x, na.rm=TRUE) {
  if (na.rm) {
    n <- colSums(!is.na(x)) # thanks @flodel
  } else {
    n <- nrow(x)
  }
  colVar <- colMeans(x*x, na.rm=na.rm) - (colMeans(x, na.rm=na.rm))^2
  return(sqrt(colVar * n/(n-1)))
}

get_treated_data <- function(MV, specs)
{
  metric = get_metric(specs$scaling)

  if (metric) {
    # standardize if all numeric
    if (specs$scaled) {
      #sd.X = sqrt((nrow(MV)-1)/nrow(MV)) * apply(MV, 2, sd)
      sd.X = sqrt((nrow(MV)-1)/nrow(MV)) * colSdColMeans(MV)
      X = scale(MV, scale=sd.X)


    } else {
      # center if all raw
      X = scale(MV, scale=FALSE)
    }
  } else {
    # standardize
#    sd.X = sqrt((nrow(MV)-1)/nrow(MV)) * apply(MV, 2, sd)
#    X = scale(MV, scale=sd.X)
    X = scale(MV) / sqrt((nrow(MV)-1)/nrow(MV))
    # all variables quantified as "ord"/"nom" are pre-treated with get_rank
    scaling = unlist(specs$scaling)
    nominal_ordinal = which(scaling %in% c("ord", "nom"))

    for (j in nominal_ordinal) {
      X[,j] = get_rank(X[,j])
    }
  }

  # add names
  dimnames(X) = dimnames(MV)
  # output
  X
}


