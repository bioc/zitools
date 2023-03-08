#'@name median
#'@title Median Value
#'
#'@param x An Object of class "Zi", NA matrix will be extracted using this function
#'@param na.rm logical, default = TRUE, NAs are excluded
#'
#'@description computes the median of the NA matrix of an S4 "Zi" class object
#'
#'@returns median value
#'
#'@example

median.Zi <- function(x, na.rm = TRUE, ...)
{
  median(x@output, na.rm = na.rm,  ...)
}



#'@name colMedians
#'
#'@title Calculates the median for each row (column) of a matrix-like object
#'
#'@param  x         An Object of class "Zi", NA matrix will be extracted using this function
#'@param  rows,     A vector indicating the subset of rows (and/or columns) to operate
#'@param  cols      over. If NULL, no subsetting is done
#'@param  na.rm     logical, default = TRUE, NAs are excluded
#'@param  ...       Additional arguments passed to specific methods
#'@param  useNames  If NA, the default behavior of the function about naming support
#'                  is remained. If FALSE, no naming support is done. Else if TRUE, names
#'                  attributes of result are set.
#'
#'@description      Calculates the median for each row (column) of a matrix-like object.
#'                  Prior to calculating the median, the NA matrix is extracted
#'
#'@return           returns a numeric vector of length N(K)
#'
#'@examples
#'
setMethod("colMedians", "Zi" , function(x,
                                        rows = NULL,
                                        cols = NULL,
                                        na.rm = TRUE,
                                        ...,
                                        useNames = NA) {
  MatrixGenerics::colMedians(
    x = x@output,
    rows = rows,
    cols = cols,
    na.rm = na.rm,
    ...,
    useNames = useNames
  )
})

setMethod("rowMedians", "Zi", function(x,
                                       rows = NULL,
                                       cols = NULL,
                                       na.rm = TRUE,
                                       ...,
                                       useNames = NA) {
  MatrixGenerics::rowMedians(
    x = x@output,
    rows = rows,
    cols = cols,
    na.rm = na.rm,
    ...,
    useNames = useNames
  )
})

#'@name quantile
#'@title Sample Quantiles
#'
#'@param x An S4 object of class "Zi", NA matrix will be extracted using this function
#'@param na.rm logical, default = TRUE, NAs are excluded
#'
#'@description computes the quantiles of the NA matrix of an S4 "Zi" class object
#'
#'@returns quantiles
#'
#'@example
quantile.Zi <- function(x, na.rm = TRUE, ...) {
  quantile(x@output, na.rm = na.rm, ...)
}

#'@name rowQuantiles
#'@title Calculates quantiles for each row (column) of a matrix-like object
#'
#'@param  x         An Object of class "Zi", NA matrix will be extracted using this function
#'@param  rows,     A vector indicating the subset of rows (and/or columns) to operate
#'@param  cols      over. If NULL, no subsetting is done
#'@param  probs     A numeric vector of J probabilities in [0,1]
#'@param  na.rm     logical, default = TRUE, NAs are excluded
#'@param  type      An integer specifying the type of estimator
#'@param  ...       Additional arguments passed to specific methods
#'@param  useNames  If NA, the default behavior of the function about naming support
#'                  is remained. If FALSE, no naming support is done. Else if TRUE, names
#'                  attributes of result are set.
#'@param  drop      If TRUE a vector is returned if J == 1.
#'
#'@description  Calculates quantiles for each row (column) of a matrix-like object
#'
#'@example

setMethod("rowQuantiles", "Zi", function(x,
                                         rows = NULL,
                                         cols = NULL,
                                         probs = seq(from = 0, to = 1, by = 0.25),
                                         na.rm = TRUE,
                                         type = 7L,
                                         ...,
                                         useNames = NA,
                                         drop = TRUE) {
  MatrixGenerics::rowQuantiles(
    x = x@output,
    rows = rows,
    cols = cols,
    probs = probs,
    na.rm = na.rm,
    type = type,
    ...,
    useNames = useNames,
    drop = drop
  )
})

setMethod("colQuantiles", "Zi", function(x,
                                         rows = NULL,
                                         cols = NULL,
                                         probs = seq(from = 0, to = 1, by = 0.25),
                                         na.rm = TRUE,
                                         type = 7L,
                                         ...,
                                         useNames = NA,
                                         drop = TRUE) {
  MatrixGenerics::colQuantiles(
    x = x@output,
    rows = rows,
    cols = cols,
    probs = probs,
    na.rm = na.rm,
    type = type,
    ...,
    useNames = useNames,
    drop = drop
  )
})


#'@name mean
#'@title Arithmetic Mean
#'@param x  An Object of class "Zi", input matrix will be used to calculate the
#'mean taking structural zero weights into account
#'@param ...
#'
#'@description  calculate the arithmetic mean of zero inflated data taking weights
#'for structural zeros into account
#'
#'@retuns value
#'@example
#'

mean.Zi <- function(zi_result, ...) {
  x <- zi_result@inputmatrix
  w <- zi_result@weights
  mean <- weighted.mean(x, w, ...)
  return(mean)
}

#'@name colMeans
#'@title Row and Column Means
#'
#'@param  x   An Object of class "Zi", input matrix will be used to calculate the
#'mean taking structural zero weights into account
#'
#'@description calculate row and column means for matrix like objects
#'
#'@returns a vector of row/col length
#'@example


setMethod("colMeans", "Zi", function(x) {
  colmean <-
    mapply(
      weighted.mean,
      as.data.frame(x@inputmatrix),
      as.data.frame(x@weights, USE.NAMES = TRUE)
    )
  return(colmean)
})

setMethod("rowMeans", "Zi", function(x) {
  rowmean <-
    mapply(weighted.mean,
           as.data.frame(t(x@inputmatrix)),
           as.data.frame(t(x@weights), USE.NAMES = TRUE))
  return(rowmean)
})

#'@name sd
#'@title Standard Deviation
#'@param x  An Object of class "Zi", input matrix will be used to calculate the
#'standard deviation taking structural zero weights into account
#'@param ...
#'
#'@description  calculate the standard deviation of zero inflated data taking weights
#'for structural zeros into account
#'
#'@retuns value
#'@example

setMethod("sd", "Zi", function(x) {
  sd <- matrixStats::weightedSd(x = x@inputmatrix,
                                w = x@weights)
  return(sd)
})

#'@name rowSds
#'@title Row and Column Standard Deviations
#'
#'@param  x   An Object of class "Zi", input matrix will be used to calculate the
#'row or col Standard Deviation taking structural zero weights into account
#'
#'@description calculate row and column sds for matrix like objects
#'
#'@returns a vector of row/col length
#'@example

setMethod("rowSds", "Zi", function(x) {
  mapply(weightedSd,
         as.data.frame(t(x@inputmatrix)),
         as.data.frame(t(x@weights), USE.NAMES = TRUE))
})
setMethod("colSds", "Zi", function(x) {
  mapply(
    weightedSd,
    as.data.frame(x@inputmatrix),
    as.data.frame(x@weights, USE.NAMES = TRUE)
  )
})

#'@name var
#'@title Variance
#'@param x  An Object of class "Zi", input matrix will be used to calculate the
#'variance taking structural zero weights into account
#'@param ...
#'
#'@description  calculate the variance of zero inflated data taking weights
#'for structural zeros into account
#'
#'@retuns value
#'@example

setMethod("var", "Zi", function(x) {
  var <- matrixStats::weightedVar(x = x@inputmatrix,
                                  w = x@weights)
  return(var)
})

#'@name rowVars
#'@title Row and Column Variances
#'
#'@param  x   An Object of class "Zi", input matrix will be used to calculate the
#'row or col Variance  taking structural zero weights into account
#'
#'@description calculate row and column variances  for matrix like objects
#'
#'@returns a vector of row/col length
#'@example

setMethod("rowVars", "Zi", function(x) {
  mapply(weightedVar,
         as.data.frame(t(x@inputmatrix)),
         as.data.frame(t(x@weights), USE.NAMES = TRUE))
})
setMethod("colVars", "Zi", function(x) {
  mapply(
    weightedVar,
    as.data.frame(x@inputmatrix),
    as.data.frame(x@weights, USE.NAMES = TRUE)
  )
})


#'@name weighted.mean
#'@title Weighted Arithmetic Mean
#'
#'@param x An Object of class "Zi", input matrix will be used to calculate the
#'meam taking structural zero weights into account
#'@param w a numerical vector of weight the same length as x giving the weights
#'to use for elements of x
#'@param ...
#'
#'@description compute a weighted mean
#'@returns value
#'@example


setMethod("weighted.mean", "Zi", function(x, w, ...) {
  mean <- weighted.mean(x@inputmatrix, w = w * x@weights)
  return(mean)
})

?rowWeightedMeans

#'@name rowWeightedMeans
#'@title Calculates the weighted mean for each row (column) of a matrix-like object
#'
#'@param x An Object of class "Zi", input matrix will be used to calculate the
#'meam taking structural zero weights into account
#'@param w a numerical vector of weight the same length as x giving the weights
#'to use for elements of x
#'@param ...
#'
#'@description Calculates the weighted mean for each row (column) of a matrix-like object.
#'@returns a numeric vector of length N(K)
#'@example
#'
setMethod("rowWeightedMeans", "Zi", function(x, w, ...) {
  mapply(weighted.mean,
         as.data.frame(t(x@inputmatrix)),
         as.data.frame((t(x@weights) * w), USE.NAMES = TRUE))
})

setMethod("colWeightedMeans", "Zi", function(x, w, ...) {
  mapply(
    weighted.mean,
    as.data.frame(x@inputmatrix),
    as.data.frame(x@weights) * w, USE.NAMES = TRUE)
})

#'@name weightedSd
#'@title Weighted Standard Deviation
#'
#'@param x An Object of class "Zi", input matrix will be used to calculate the
#'meam taking structural zero weights into account
#'@param w a numerical vector of weight the same length as x giving the weights
#'to use for elements of x
#'@param ...
#'
#'@description compute a weighted standard deviation
#'@returns value
#'@example

setMethod("weightedSd", "Zi", function(x, w, ...) {
  sqrt(weightedVar(x=x,w=w,...))
})

#'@name rowWeightedSds
#'@title Calculates the weighted mean for each row (column) of a matrix-like object
#'
#'@param x An Object of class "Zi", input matrix will be used to calculate the
#'row weighted Standard Deviation taking structural zero weights into account
#'@param w a numerical vector of weight the same length as x giving the weights
#'to use for elements of x
#'@param ...
#'
#'@description Calculates the row weighted mean for each row (column) of a matrix-like object.
#'@returns a numeric vector of length N(K)
#'@example
setMethod("rowWeightedSds", "Zi", function(x, w, ...) {
  mapply(weightedSd,
         as.data.frame(t(x@inputmatrix)),
         as.data.frame((t(x@weights) * w), USE.NAMES = TRUE))
})
setMethod("colWeightedSds", "Zi", function(x, w, ...) {
  mapply(
    weightedSd,
    as.data.frame(x@inputmatrix),
    as.data.frame(x@weights * w, USE.NAMES = TRUE)
  )
})

#'@name weightedVar
#'@title Weighted Variance and weighted Standard Deviation
#'
#'@param x An Object of class "Zi", input matrix will be used to calculate the
#'row weighted Standard Deviation taking structural zero weights into account
#'@param w a vector of weights the same length as x giving the weights to use for
#'each element of x. Negative weights are treated as zero weights. Default value
#'is equal weight to all values.
#'@description Computes a weighted variance / standard deviation of a numeric
#'vector or across rows or columns of a matrix
#'@returns a numeric scalar
#'@example

setMethod("weightedVar", "Zi", function(x, w, ...) {
  weightedVar(x = x@inputmatrix,
              w = w * x@weights)
})

setMethod("rowWeightedVars", "Zi", function(x, w, ...) {
  mapply(weightedVar,
         as.data.frame(t(x@inputmatrix)),
         as.data.frame((t(x@weights) * w), USE.NAMES = TRUE))
})
setMethod("colWeightedVars", "Zi", function(x, w, ...) {
  mapply(
    weightedVar,
    as.data.frame(x@inputmatrix),
    as.data.frame(x@weights * w, USE.NAMES = TRUE))
})

?ZiMain

