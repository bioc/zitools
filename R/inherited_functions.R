#'@export
#'@name median
#'@title Median Value
#'
#'@param x 'Zi'-class object, output matrix will be extracted using this function
#'@param na.rm logical, default = TRUE, NAs are excluded
#'
#'@description computes the median of the NA matrix of an S4 "Zi" class object
#'
#'@returns median value
#'@importFrom stats median
#'@example

median.Zi <- function(x, na.rm = TRUE, ...)
{
  median(x@output, na.rm = na.rm,  ...)
}



#'@export
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
#'@importFrom MatrixGenerics colMedians
#'@examples
#'
setMethod("colMedians", "Zi" , function(x,
                                        rows = NULL,
                                        cols = NULL,
                                        na.rm = TRUE,
                                        ...,
                                        useNames = NA) {
  colMedians(
    x = x@output,
    rows = rows,
    cols = cols,
    na.rm = na.rm,
    ...,
    useNames = useNames
  )
})

#'@name rowMedians
#'@export
#'@rdname colMedians
#'@importFrom MatrixGenerics rowMedians
setMethod("rowMedians", "Zi", function(x,
                                       rows = NULL,
                                       cols = NULL,
                                       na.rm = TRUE,
                                       ...,
                                       useNames = NA) {
  rowMedians(
    x = x@output,
    rows = rows,
    cols = cols,
    na.rm = na.rm,
    ...,
    useNames = useNames
  )
})

#'@export
#'@name quantile
#'@title Sample Quantiles
#'
#'@param x An S4 object of class "Zi", NA matrix will be extracted using this function
#'@param na.rm logical, default = TRUE, NAs are excluded
#'
#'@description computes the quantiles of the NA matrix of an S4 "Zi" class object
#'
#'@returns quantiles
#'@importFrom stats quantile
#'@example

quantile.Zi <- function(x, na.rm = TRUE, ...) {
  quantile(x@output, na.rm = na.rm, ...)
}

#'@export
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
#'@importFrom MatrixGenerics rowQuantiles
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
  rowQuantiles(
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

#'@name colQuantiles
#'@export
#'@rdname rowQuantiles
#'@importFrom MatrixGenerics colQuantiles

setMethod("colQuantiles", "Zi", function(x,
                                         rows = NULL,
                                         cols = NULL,
                                         probs = seq(from = 0, to = 1, by = 0.25),
                                         na.rm = TRUE,
                                         type = 7L,
                                         ...,
                                         useNames = NA,
                                         drop = TRUE) {
  colQuantiles(
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


#'@export
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
  x <- zi_result@countmatrix
  w <- zi_result@weights
  mean <- weighted.mean(x, w, ...)
  return(mean)
}

#'@export
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
      as.data.frame(x@countmatrix),
      as.data.frame(x@weights, USE.NAMES = TRUE)
    )
  return(colmean)
})

#'@name rowMeans
#'@export
#'@rdname colMeans
#'
setMethod("rowMeans", "Zi", function(x) {
  rowmean <-
    mapply(weighted.mean,
           as.data.frame(t(x@countmatrix)),
           as.data.frame(t(x@weights), USE.NAMES = TRUE))
  return(rowmean)
})


#'@export
#'@name sd
#'@title Standard Deviation
#'@param x  An Object of class "Zi", input matrix will be used to calculate the
#'standard deviation taking structural zero weights into account
#'@param ...
#'
#'@description  calculate the standard deviation of zero inflated data taking weights
#'for structural zeros into account
#'@importFrom matrixStats weightedSd
#'@retuns value
#'@example

setMethod("sd", "Zi", function(x) {
  sd <- weightedSd(x = x@countmatrix,
                                w = x@weights)
  return(sd)
})

#'@export
#'@name rowSds
#'@title Row and Column Standard Deviations
#'
#'@param  x   An Object of class "Zi", input matrix will be used to calculate the
#'row or col Standard Deviation taking structural zero weights into account
#'
#'@description calculate row and column sds for matrix like objects
#'
#'@returns a vector of row/col length
#'@importFrom matrixStats weightedSd
#'@importFrom MatrixGenerics rowSds
#'@example

setMethod("rowSds", "Zi", function(x) {
  mapply(weightedSd,
         as.data.frame(t(x@countmatrix)),
         as.data.frame(t(x@weights), USE.NAMES = TRUE))
})

#'@name colSds
#'@export
#'@rdname rowSds
#'@importFrom matrixStats weightedSd
#'@importFrom MatrixGenerics colSds
#'
setMethod("colSds", "Zi", function(x) {
  mapply(
    weightedSd,
    as.data.frame(x@countmatrix),
    as.data.frame(x@weights, USE.NAMES = TRUE)
  )
})

#'@export
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
#'@importFrom matrixStats weightedVar
#'@example

setMethod("var", "Zi", function(x) {
  var <- weightedVar(x = x@countmatrix,
                                  w = x@weights)
  return(var)
})

#'@export
#'@name rowVars
#'@title Row and Column Variances
#'
#'@param  x   An Object of class "Zi", input matrix will be used to calculate the
#'row or col Variance  taking structural zero weights into account
#'
#'@description calculate row and column variances  for matrix like objects
#'
#'@returns a vector of row/col length
#'@importFrom matrixStats weightedVar
#'@importFrom MatrixGenerics rowVars
#'@example

setMethod("rowVars", "Zi", function(x) {
  mapply(weightedVar,
         as.data.frame(t(x@countmatrix)),
         as.data.frame(t(x@weights), USE.NAMES = TRUE))
})

#'@name colVars
#'@export
#'@rdname rowVars
#'@importFrom matrixStats weightedVar
#'@importFrom MatrixGenerics colVars

setMethod("colVars", "Zi", function(x) {
  mapply(
    weightedVar,
    as.data.frame(x@countmatrix),
    as.data.frame(x@weights, USE.NAMES = TRUE)
  )
})


#'@export
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
#'@importFrom stats weighted.mean
#'@example


setMethod("weighted.mean", "Zi", function(x, w, ...) {
  mean <- weighted.mean(x@countmatrix, w = w * x@weights)
  return(mean)
})


#'@export
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
#'@importFrom stats weighted.mean
#'@importFrom MatrixGenerics rowWeightedMeans
#'@example
#'
setMethod("rowWeightedMeans", "Zi", function(x, w, ...) {
  mapply(weighted.mean,
         as.data.frame(t(x@countmatrix)),
         as.data.frame((t(x@weights) * w), USE.NAMES = TRUE))
})

#'@name colWeightedMeans
#'@export
#'@rdname rowWeightedMeans
#'@importFrom stats weighted.mean
#'@importFrom MatrixGenerics colWeightedMeans

setMethod("colWeightedMeans", "Zi", function(x, w, ...) {
  mapply(
    weighted.mean,
    as.data.frame(x@countmatrix),
    as.data.frame(x@weights) * w, USE.NAMES = TRUE)
})


#setMethod("weightedSd", "Zi", function(x, w, ...) {
  #sqrt(weightedVar(x=x,w=w,...))
#})

#'@export
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
#'@importFrom matrixStats weightedSd
#'@importFrom MatrixGenerics rowWeightedSds
#'@example

setMethod("rowWeightedSds", "Zi", function(x, w, ...) {
  mapply(weightedSd,
         as.data.frame(t(x@countmatrix)),
         as.data.frame((t(x@weights) * w), USE.NAMES = TRUE))
})

#'@name colWeightedSds
#'@export
#'@rdname rowWeightedSds
#'@importFrom matrixStats weightedSd
#'@importFrom MatrixGenerics colWeightedSds

setMethod("colWeightedSds", "Zi", function(x, w, ...) {
  mapply(
    weightedSd,
    as.data.frame(x@countmatrix),
    as.data.frame(x@weights * w, USE.NAMES = TRUE)
  )
})

#'@export
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
#'@importFrom matrixStats weightedVar
#'@example

setGeneric("weightedVar", function(x, w = NULL, idxs = NULL, na.rm = FALSE, center = NULL,
                                   ...) standardGeneric("weightedVar"))

setMethod("weightedVar", "Zi", function(x, w, ...) {
  weightedVar(x = x@countmatrix,
              w = w * x@weights)
})

#'@export
#'@name rowWeightedVars
#'@title Calculates the weighted variance for each row (column) of a matrix-like object
#'
#'@param x An Object of class "Zi", input matrix will be used to calculate the
#'row weighted variances taking structural zero weights into account
#'@param w a numerical vector of weight the same length as x giving the weights
#'to use for elements of x
#'@param ...
#'
#'@description Calculates the row weighted variances for each row (column) of a matrix-like object.
#'@returns a numeric vector of length N(K)
#'@importFrom matrixStats weightedVar
#'@importFrom MatrixGenerics rowWeightedVars
#'@example
setMethod("rowWeightedVars", "Zi", function(x, w, ...) {
  mapply(weightedVar,
         as.data.frame(t(x@countmatrix)),
         as.data.frame((t(x@weights) * w), USE.NAMES = TRUE))
})

#'@name colWeightedVars
#'@export
#'@rdname rowWeightedVars
#'@importFrom matrixStats weightedVar
#'@importFrom MatrixGenerics colWeightedVars

setMethod("colWeightedVars", "Zi", function(x, w, ...) {
  mapply(
    weightedVar,
    as.data.frame(x@countmatrix),
    as.data.frame(x@weights * w, USE.NAMES = TRUE))
})

#'@name t()
#'@title Matrix Transpose
#'
#'@export
#'
setMethod(
  "t",
  signature = "Zi",
  definition = function(x) {
    countmatrix <- t(x@countmatrix)
    output <- t(x@output)
    weights <- t(x@weights)
    result <- new(
      Class = "Zi",
      datafile = x@datafile,
      countmatrix = countmatrix,
      ZINBModel = x@ZINBModel,
      output = output,
      weights = weights
    )
    return(result)
  }
)
