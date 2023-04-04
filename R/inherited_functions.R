#'@include zi_function.R
NULL

#'@export
#'@name median
#'@title Calculate the median of zero-deinflated count data
#'
#'@param x \code{\linkS4class{Zi}}-class object
#'@param na.rm \code{\link[base]{logical}} If \code{\link{TRUE}} NAs are excluded, otherwise
#'not. default = \code{\link{TRUE}}
#'@param ... see \link[stats]{median.default}
#'@description Caluclate the median of  zero-deinflated data of a 'Zi'-class
#'object. To calculate the median, the output matrix will be extracted
#'
#'@returns median value
#'@importFrom stats median
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'median(Zi)

median.Zi <- function(x, na.rm = TRUE, ...)
{
  median(x@output, na.rm = na.rm,  ...)
}



#'@export
#'@name colMedians
#'
#'@title Calculate the row or column median of zero-deinflated count data
#'@description Calculate the row or column median of  zero-deinflated data of a
#' \code{\linkS4class{Zi}}-class object. To calculate the median, the output matrix will be extracted
#'@param  x         \code{\linkS4class{Zi}}-class object
#'@param  rows,     A  \code{\link[base]{vector}} indicating the subset of rows (and/or columns) to operate
#'@param  cols      over. If  \code{\link{NULL}}, no subsetting is done
#'@param  na.rm     \code{\link[base]{logical}} If \code{\link{TRUE}} \code{\link{NA}}s are excluded, otherwise
#'not. default = \code{\link{TRUE}}
#'@param  useNames  \code{\link[base]{logical}}.  If \code{\link{TRUE}} (default), names
#'attributes of result are set. Else if \code{\link{FALSE}}, no naming support is done.
#'
#'@returns returns a numeric vector of row/column length
#'@importFrom MatrixGenerics colMedians
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'colMedians(Zi, useNames = TRUE)
#'rowMedians(Zi, useNames = TRUE)
#'
setMethod("colMedians", "Zi" , function(x,
                                        rows = NULL,
                                        cols = NULL,
                                        na.rm = TRUE,
                                        ...,
                                        useNames = TRUE) {
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
                                       useNames = TRUE) {
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
#'@title Calculate the quantiles of zero-deinflated count data
#'
#'@param x A \code{\linkS4class{Zi}}-class object
#'@param na.rm \code{\link[base]{logical}} If \code{\link{TRUE}} NAs are excluded, otherwise
#'not. default = \code{\link{TRUE}}
#'@param ... \link[stats]{quantile}
#'@description Calculate the quantiles of  zero-deinflated data of a
#'\code{\linkS4class{Zi}}-class object. To calculate the quantiles, the output
#'matrix will be extracted.
#'
#'@returns quantile values
#'@importFrom stats quantile
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'quantile(Zi)

quantile.Zi <- function(x, probs = seq(0, 1, 0.25), na.rm = TRUE, ...) {
  quantile(x@output, probs = probs, na.rm = na.rm, ...)
}

#'@export
#'@name rowQuantiles
#'@title Calculate the row or column quantiles of zero-deinflated count data
#'
#'@param  x         A \code{\linkS4class{Zi}}-class object
#'@param  rows,cols A \code{\link[base]{vector}} indicating the subset of rows
#'and/or columns to operate over. If \code{\link{NULL}} (default), no subsetting is done
#'@param  probs     A numeric \code{\link[base]{vector}} of J probabilities in [0,1]
#'@param  na.rm     \code{\link[base]{logical}} If \code{\link{TRUE}}
#'\code{\link{NA}}s are excluded, otherwise not. default = \code{\link{TRUE}}
#'@param  type      An integer specifying the type of estimator
#'@param  ...       Additional arguments passed to specific methods
#'\link[MatrixGenerics]{rowQuantiles}
#'@param  useNames  \code{\link[base]{logical}}  If \code{\link{TRUE}} (default), names
#'attributes of result are set. Else if \code{\link{FALSE}}, no naming support is done.
#'@param  drop      If \code{\link{TRUE}} a \code{\link[base]{vector}} is returned if J == 1.
#'
#'@description Calculate the row or column quantiles of  zero-deinflated data of a
#' \code{\linkS4class{Zi}}-class object. To calculate the quantiles, the output matrix will be extracted
#'@importFrom MatrixGenerics rowQuantiles
#'
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'rowQuantile(Zi, useNames = TRUE)
#'colQuantiles(Zi)

setMethod("rowQuantiles", "Zi", function(x,
                                         rows = NULL,
                                         cols = NULL,
                                         probs = seq(from = 0, to = 1, by = 0.25),
                                         na.rm = TRUE,
                                         type = 7L,
                                         ...,
                                         useNames = TRUE,
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
#'@param x  A \code{\linkS4class{Zi}}-class object
#'@param ... \link[base]{mean.default}
#'
#'@description  Calculate the arithmetic mean of zero inflated data taking weights
#'for structural zeros into account
#'
#'@returns value
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'mean(Zi)
#'

mean.Zi <- function(x, ...) {
  Zi <- x
  x <- Zi@countmatrix
  w <- Zi@weights
  mean <- weighted.mean(x, w, ...)
  return(mean)
}

#'@export
#'@name colMeans2
#'@title Calculate the row or column means of zero-inflated count data
#'
#'@param  x   A \code{\linkS4class{Zi}}-class object
#'@param  rows,cols A \code{\link[base]{vector}} indicating the subset of rows and/or columns to
#'operate over. If \code{\link{NULL}} (default), no subsetting is done
#'@param na.rm \code{\link[base]{logical}} If \code{\link{TRUE}} \code{\link{NA}}s are excluded, otherwise
#'not. default = \code{\link{FALSE}}
#'@param useNames \code{\link[base]{logical}}  If \code{\link{TRUE}} (default), names
#'attributes of result are set. Else if \code{\link{FALSE}}, no naming support is done.
#'
#'@description Calculate row and column means of zero-inflated count data taking
#'weights for structural zeros into account.
#'
#'@returns a numeric \code{\link[base]{vector}} of row/column length
#'@importFrom stats weighted.mean
#'@importFrom MatrixGenerics colMeans2
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'colMeans2(Zi)
#'rowMeans2(Zi)


setMethod("colMeans2", "Zi", function(x, rows = NULL, cols = NULL, na.rm = FALSE, useNames=TRUE) {
  if(is.null(rows)){
    rows <-  c(1:nrow(x@countmatrix))}
  if(is.null(cols)){
    cols <-  c(1:ncol(x@countmatrix))}
  colmean <-
    mapply(
      weighted.mean,
      as.data.frame(x@countmatrix)[rows,cols],
      as.data.frame(x@weights)[rows, cols],USE.NAMES = useNames, na.rm = na.rm)
  return(colmean)
})

#'@name rowMeans2
#'@export
#'@rdname colMeans2
#'@importFrom stats weighted.mean
#'@importFrom MatrixGenerics rowMeans2
#'
setMethod("rowMeans2", "Zi", function(x, rows = NULL, cols = NULL, na.rm = FALSE, useNames=TRUE) {
  if(is.null(rows)){
    rows <- c(1:nrow(t(x@countmatrix)))}
  if(is.null(cols)){
    cols <- c(1:ncol(t(x@countmatrix)))}
  rowmean <-
    mapply(weighted.mean,
           as.data.frame(t(x@countmatrix))[rows,cols],
           as.data.frame(t(x@weights))[rows,cols], USE.NAMES = useNames, na.rm=na.rm)
  return(rowmean)
})


#'@export
#'@name sd
#'@title Standard Deviation of zero inflated count data
#'@param x  A \code{\linkS4class{Zi}}-class object
#'@param na.rm \code{\link[base]{logical}} If \code{\link{TRUE}} \code{\link{NA}}s
#'are excluded, otherwise not. default = \code{\link{FALSE}}
#'
#'@description  Calculate the standard deviation of zero inflated count data
#'taking weights for structural zeros into account.
#'@importFrom matrixStats weightedSd
#'@importFrom stats sd
#'@returns value
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'sd(Zi)

setMethod("sd", "Zi", function(x, na.rm = FALSE) {
  sd <- weightedSd(x = x@countmatrix,
                                w = x@weights, na.rm = na.rm)
  return(sd)
})

#'@export
#'@name rowSds
#'@title Row and Column Standard Deviations of zero inflated count data
#'
#'@param  x   A \code{\linkS4class{Zi}}-class object
#'@param  rows,cols A \code{\link[base]{vector}} indicating the subset of rows and/or columns to
#'operate over. If \code{\link{NULL}} (default), no subsetting is done
#'@param na.rm \code{\link[base]{logical}} If \code{\link{TRUE}} \code{\link{NA}}s
#'are excluded, otherwise not. default = \code{\link{FALSE}}
#'@param useNames \code{\link[base]{logical}}  If \code{\link{TRUE}} (default), names
#'attributes of result are set. Else if \code{\link{FALSE}}, no naming support is done.
#'
#'@description Calculate row and column standard deviations of zero inflated
#'count data taking weights for structural zeros into account
#'
#'@returns a \code{\link[base]{vector}} of row/column length
#'@importFrom matrixStats weightedSd
#'@importFrom MatrixGenerics rowSds
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'rowSds(Zi)
#'colSds(Zi)

setMethod("rowSds", "Zi", function(x, rows = NULL, cols = NULL, na.rm = FALSE, useNames = TRUE) {
  if(is.null(rows)){
    rows <- c(1:nrow(t(x@countmatrix)))}
  if(is.null(cols)){
    cols <- c(1:ncol(t(x@countmatrix)))}
  mapply(weightedSd,
         as.data.frame(t(x@countmatrix))[rows,cols],
         as.data.frame(t(x@weights))[rows,cols], USE.NAMES = useNames, na.rm = na.rm)
})

#'@name colSds
#'@export
#'@rdname rowSds
#'@importFrom matrixStats weightedSd
#'@importFrom MatrixGenerics colSds

setMethod("colSds", "Zi", function(x, rows = NULL, cols = NULL, na.rm = FALSE, useNames = TRUE) {
  if(is.null(rows)){
    rows <-  c(1:nrow(x@countmatrix))}
  if(is.null(cols)){
    cols <-  c(1:ncol(x@countmatrix))}
  mapply(
    weightedSd,
    as.data.frame(x@countmatrix)[rows,cols],
    as.data.frame(x@weights)[rows,cols], USE.NAMES = useNames, na.rm = na.rm)
})

#'@export
#'@name var
#'@title Variance of zero inflated count data
#'@param x  A \code{\linkS4class{Zi}}-class object
#'@param na.rm \code{\link[base]{logical}} If \code{\link{TRUE}} \code{\link{NA}}s
#'are excluded, otherwise not. default = \code{\link{FALSE}}
#'
#'@description  Calculate the variance of zero inflated count data taking weights
#'for structural zeros into account.
#'
#'@returns value
#'@importFrom matrixStats weightedVar
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'var(Zi)

setMethod("var", "Zi", function(x, na.rm = FALSE) {
  var <- weightedVar(x = x@countmatrix,
                                  w = x@weights, na.rm = na.rm)
  return(var)
})

#'@export
#'@name rowVars
#'@title Row and Column Variances of zero inflated count data
#'
#'@param  x   A \code{\linkS4class{Zi}}-class object
#'@param  rows,cols A \code{\link[base]{vector}} indicating the subset of rows and/or columns to
#'operate over. If \code{\link{NULL}} (default), no subsetting is done
#'@param na.rm \code{\link[base]{logical}} If \code{\link{TRUE}} \code{\link{NA}}s
#'are excluded, otherwise not. default = \code{\link{FALSE}}
#'@param useNames \code{\link[base]{logical}}  If \code{\link{TRUE}} (default), names
#'attributes of result are set. Else if \code{\link{FALSE}}, no naming support is done.
#'
#'@description Calculate row and column variances of zero inflated count data
#'taking weights for structural zeros into account.
#'
#'@returns a vector of row/col length
#'@importFrom matrixStats weightedVar
#'@importFrom MatrixGenerics rowVars
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'rowVars(Zi)
#'colVars(Zi)

setMethod("rowVars", "Zi", function(x, rows = NULL, cols = NULL, na.rm = FALSE, useNames = TRUE) {
  if(is.null(rows)){
    rows <- c(1:nrow(t(x@countmatrix)))}
  if(is.null(cols)){
    cols <- c(1:ncol(t(x@countmatrix)))}
  mapply(weightedVar,
         as.data.frame(t(x@countmatrix))[rows,cols],
         as.data.frame(t(x@weights))[rows,cols], USE.NAMES = useNames, na.rm = na.rm)
})

#'@name colVars
#'@export
#'@rdname rowVars
#'@importFrom matrixStats weightedVar
#'@importFrom MatrixGenerics colVars

setMethod("colVars", "Zi", function(x, rows = NULL, cols = NULL, na.rm = FALSE, useNames = TRUE) {
  if(is.null(rows)){
    rows <-  c(1:nrow(x@countmatrix))}
  if(is.null(cols)){
    cols <-  c(1:ncol(x@countmatrix))}
  mapply(
    weightedVar,
    as.data.frame(x@countmatrix)[rows,cols],
    as.data.frame(x@weights)[rows,cols], USE.NAMES = useNames, na.rm = na.rm)
})

#'@export
#'@name weighted.mean
#'@title Weighted Arithmetic Mean of zero inflated count data
#'
#'@param x A \code{\linkS4class{Zi}}-class object
#'@param w a numerical \code{\link[base]{vector}} of weight the same length as
#'x giving the weights to use for elements of x
#'@param ... \link[stats]{weighted.mean}
#'
#'@description Calculate a weighted mean of zero inflated count data, additionally
#' taking weights for structural zeros into account
#'@returns value
#'@importFrom stats weighted.mean
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'weight <- runif(length(Zi@countmatrix), 0.1, 1)
#'weighted.mean(Zi, w= weight)


setMethod("weighted.mean", "Zi", function(x, w, ...) {
  mean <- weighted.mean(x@countmatrix, w = w * x@weights, ...)
  return(mean)
})

#'@export
#'@name rowWeightedMeans
#'@title Row and Column weighted means of zero inflated count data
#'
#'@param x A \code{\linkS4class{Zi}}-class object
#'@param w a numerical vector of weights either of length = rows or length = cols
#' giving the weights to use for elements of x
#'@param  rows,cols A \code{\link[base]{vector}} indicating the subset of rows and/or columns to
#'operate over. If \code{\link{NULL}} (default), no subsetting is done
#'@param na.rm \code{\link[base]{logical}} If \code{\link{TRUE}} \code{\link{NA}}s are excluded, otherwise
#'not. default = \code{\link{FALSE}}
#'@param useNames \code{\link[base]{logical}}  If \code{\link{TRUE}} (default), names
#'attributes of result are set. Else if \code{\link{FALSE}}, no naming support is done.
#'
#'@description Calculate row and column weighted means of zero inflated count
#'data, additionally taking weights for structural zeros into account.
#'@returns a numeric vector of length N(K)
#'@importFrom stats weighted.mean
#'@importFrom MatrixGenerics rowWeightedMeans
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'rowWeightedMeans(Zi, w = runif(ncol(Zi@countmatrix), 0.1,1))
#'colWeightedMeans(Zi, w = runif(nrow(Zi@countmatrix), 0.1,1))
#'

setMethod("rowWeightedMeans", "Zi", function(x,
                                             w,
                                             rows = NULL,
                                             cols = NULL,
                                             na.rm = FALSE,
                                             useNames = TRUE) {
  if (is.null(rows)) {
    rows <- c(1:nrow(t(x@countmatrix)))
  }
  if (is.null(cols)) {
    cols <- c(1:ncol(t(x@countmatrix)))
  }
  mapply(
    weighted.mean,
    as.data.frame(t(x@countmatrix))[rows, cols],
    as.data.frame(t(x@weights * w))[rows, cols],
    USE.NAMES = useNames,
    na.rm = na.rm
  )
})


#'@name colWeightedMeans
#'@export
#'@rdname rowWeightedMeans
#'@importFrom stats weighted.mean
#'@importFrom MatrixGenerics colWeightedMeans

setMethod("colWeightedMeans", "Zi", function(x, w, rows = NULL, cols = NULL, na.rm = FALSE, useNames = TRUE) {
  if(is.null(rows)){
    rows <-  c(1:nrow(x@countmatrix))}
  if(is.null(cols)){
    cols <-  c(1:ncol(x@countmatrix))}
  mapply(
    weighted.mean,
    as.data.frame(x@countmatrix)[rows,cols],
    as.data.frame(x@weights*w)[rows,cols], USE.NAMES = useNames, na.rm = na.rm)
})

#setMethod("weightedSd", "Zi", function(x, w, ...) {
  #sqrt(weightedVar(x=x,w=w,...))
#})

#'@export
#'@name rowWeightedSds
#'@title Row and column weighted standard deviations or variances of zero
#'inflated count data
#'
#'@param x A \code{\linkS4class{Zi}}-class object
#'@param w a numerical vector of weights either of length = rows or length = cols
#' giving the weights to use for elements of x
#'@param  rows,cols A \code{\link[base]{vector}} indicating the subset of rows and/or columns to
#'operate over. If \code{\link{NULL}} (default), no subsetting is done
#'@param na.rm \code{\link[base]{logical}} If \code{\link{TRUE}} \code{\link{NA}}s are excluded, otherwise
#'not. default = \code{\link{FALSE}}
#'@param useNames \code{\link[base]{logical}}  If \code{\link{TRUE}} (default), names
#'attributes of result are set. Else if \code{\link{FALSE}}, no naming support is done.
#'
#'@description Calculate row and column standard deviations or variances of
#'zero inflated count data, additionally taking weights for structural zeros
#'into account.
#'@returns a numeric vector of length N(K)
#'@importFrom matrixStats weightedSd
#'@importFrom MatrixGenerics rowWeightedSds
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'rowWeightedSds(Zi, w = runif(ncol(Zi@countmatrix), 0.1,1))
#'colWeightedSds(Zi, w = runif(nrow(Zi@countmatrix), 0.1,1))
#'rowWeightedVars(Zi, w = runif(ncol(Zi@countmatrix), 0.1,1))
#'colWeightedVars(Zi, w = runif(nrow(Zi@countmatrix), 0.1,1)
#'

setMethod("rowWeightedSds", "Zi", function(x,
                                             w,
                                             rows = NULL,
                                             cols = NULL,
                                             na.rm = FALSE,
                                             useNames = TRUE) {
  if (is.null(rows)) {
    rows <- c(1:nrow(t(x@countmatrix)))
  }
  if (is.null(cols)) {
    cols <- c(1:ncol(t(x@countmatrix)))
  }
  mapply(
    weightedSd,
    as.data.frame(t(x@countmatrix))[rows, cols],
    as.data.frame(t(x@weights * w))[rows, cols],
    USE.NAMES = useNames,
    na.rm = na.rm
  )
})

#'@name colWeightedSds
#'@export
#'@rdname rowWeightedSds
#'@importFrom matrixStats weightedSd
#'@importFrom MatrixGenerics colWeightedSds

setMethod("colWeightedSds", "Zi", function(x, w, rows = NULL, cols = NULL, na.rm = FALSE, useNames = TRUE) {
  if(is.null(rows)){
    rows <-  c(1:nrow(x@countmatrix))}
  if(is.null(cols)){
    cols <-  c(1:ncol(x@countmatrix))}
  mapply(
    weightedSd,
    as.data.frame(x@countmatrix)[rows,cols],
    as.data.frame(x@weights*w)[rows,cols], USE.NAMES = useNames, na.rm = na.rm)
})

#'@export
#'@name weightedVar
#'@title Weighted Variance and weighted Standard Deviation
#'
#'@param x A \code{\linkS4class{Zi}}-class object
#'@param w a numerical \code{\link[base]{vector}} of weight the same length as
#'x giving the weights to use for elements of x
#'@param idxs A \code{\link[base]{vector}} indicating subset of elements to
#'operate over. If \code{\link{NULL}}, no subsetting is done.
#'@param na.rm \code{\link[base]{logical}} If \code{\link{TRUE}} \code{\link{NA}}s are excluded, otherwise
#'not. default = \code{\link{FALSE}}
#'@param center	 \code{\link{numeric}} scalar specifying the center location of
#' the data. If \code{\link{NULL}}, it is estimated from data.
#'@description  Calculate a weighted variance of zero inflated count data,
#'additionally taking weights for structural zeros into account
#'@returns a \code{\link{numeric}} scalar
#'@importFrom matrixStats weightedVar
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'weight <- runif(length(Zi@countmatrix), 0.1, 1)
#'weightedVar(Zi, w= weight)

setGeneric("weightedVar", function(x, w = NULL, idxs = NULL, na.rm = FALSE, center = NULL,
                                   ...) standardGeneric("weightedVar"))

setMethod("weightedVar", "Zi", function(x, w, idxs = NULL, na.rm = FALSE, center = NULL, ...) {
  weightedVar(x = x@countmatrix,
              w = w * x@weights, idxs = idxs, na.rm = na.rm, center = NULL, ...)
})

#'@export
#'@name rowWeightedVars
#'@rdname rowWeightedSds
#'@importFrom matrixStats weightedVar
#'@importFrom MatrixGenerics rowWeightedVars

setMethod("rowWeightedVars", "Zi", function(x,
                                           w,
                                           rows = NULL,
                                           cols = NULL,
                                           na.rm = FALSE,
                                           useNames = TRUE) {
  if (is.null(rows)) {
    rows <- c(1:nrow(t(x@countmatrix)))
  }
  if (is.null(cols)) {
    cols <- c(1:ncol(t(x@countmatrix)))
  }
  mapply(
    weightedVar,
    as.data.frame(t(x@countmatrix))[rows, cols],
    as.data.frame(t(x@weights * w))[rows, cols],
    USE.NAMES = useNames,
    na.rm = na.rm
  )
})

#'@name colWeightedVars
#'@export
#'@rdname rowWeightedSds
#'@importFrom matrixStats weightedVar
#'@importFrom MatrixGenerics colWeightedVars


setMethod("colWeightedVars", "Zi", function(x, w, rows = NULL, cols = NULL, na.rm = FALSE, useNames = TRUE) {
  if(is.null(rows)){
    rows <-  c(1:nrow(x@countmatrix))}
  if(is.null(cols)){
    cols <-  c(1:ncol(x@countmatrix))}
  mapply(
    weightedVar,
    as.data.frame(x@countmatrix)[rows,cols],
    as.data.frame(x@weights*w)[rows,cols], USE.NAMES = useNames, na.rm = na.rm)
})


#'@name log1p
#'@title log(1+x)
#'@description Calculate log(1+x) of all 'matrix' objects of a 'Zi'-class
#'object, log calculates by default natural logarithms
#'@param x \code{\linkS4class{Zi}}-class object
#'@export
#'@seealso \link[base]{log1p}
#'
#'
setMethod("log1p", signature ="Zi", definition = function(x){
  countmatrix <- log1p(x@countmatrix)
  output <- log1p(x@output)
  weights <- log1p(x@weights)
  result <- new(
    Class = "Zi",
    inputdata = x@inputdata,
    countmatrix = countmatrix,
    ZiModel = x@ZiModel,
    output = output,
    weights = weights)
})
