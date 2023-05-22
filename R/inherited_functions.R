#'@include ziMain.R
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
#'object. To calculate the median, the deinflatedcounts matrix will be extracted
#'
#'@returns median value
#'@importFrom stats median
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'median(Zi)
#'@seealso \link[stats]{median}, \link[zitools]{colMedians}, \link[zitools]{rowMedians}

median.Zi <- function(x, na.rm = TRUE, ...)
{
  median(x@deinflatedcounts, na.rm = na.rm,  ...)
}



#'@export
#'@name colMedians
#'@aliases colMedians,Zi-method
#'@title Calculate the row or column median of zero-deinflated count data
#'@description Calculate the row or column median of  zero-deinflated data of a
#' \code{\linkS4class{Zi}}-class object. To calculate the median, the deinflatedcounts matrix will be extracted
#'@param  x         \code{\linkS4class{Zi}}-class object
#'@param  rows,     A  \code{\link[base]{vector}} indicating the subset of rows (and/or columns) to operate
#'@param  cols      over. If  \code{\link{NULL}}, no subsetting is done
#'@param  na.rm     \code{\link[base]{logical}} If \code{\link{TRUE}} \code{\link{NA}}s are excluded, otherwise
#'not. default = \code{\link{TRUE}}
#'@param  useNames  \code{\link[base]{logical}}.  If \code{\link{TRUE}} (default), names
#'attributes of result are set. Else if \code{\link{FALSE}}, no naming support is done.
#'@param ... see \code{\link[MatrixGenerics]{colMedians}}
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
    x = x@deinflatedcounts,
    rows = rows,
    cols = cols,
    na.rm = na.rm,
    ...,
    useNames = useNames
  )
})

#'@name rowMedians
#'@aliases rowMedians,Zi-method
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
    x = x@deinflatedcounts,
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
#'@param probs A numeric \code{\link[base]{vector}} of J probabilities in \[0,1\]
#'@param ... \link[stats]{quantile}
#'@description Calculate the quantiles of  zero-deinflated data of a
#'\code{\linkS4class{Zi}}-class object. To calculate the quantiles, the deinflatedcounts
#'matrix will be extracted.
#'
#'@returns quantile value
#'@importFrom stats quantile
#'@seealso \link[stats]{quantile}, \link[zitools]{rowQuantiles}, \link[zitools]{colQuantiles}
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'quantile(Zi)

quantile.Zi <- function(x, probs = seq(0, 1, 0.25), na.rm = TRUE, ...) {
  quantile(x@deinflatedcounts, probs = probs, na.rm = na.rm, ...)
}

#'@export
#'@name rowQuantiles
#'@aliases rowQuantiles,Zi-method
#'@title Calculate the row or column quantiles of zero-deinflated count data
#'
#'@param  x         A \code{\linkS4class{Zi}}-class object
#'@param  rows,cols A \code{\link[base]{vector}} indicating the subset of rows
#'and/or columns to operate over. If \code{\link{NULL}} (default), no subsetting is done
#'@param  probs     A numeric \code{\link[base]{vector}} of J probabilities in \[0,1\]
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
#' \code{\linkS4class{Zi}}-class object. To calculate the quantiles, the deinflatedcounts matrix will be extracted
#'@importFrom MatrixGenerics rowQuantiles
#'
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'rowQuantiles(Zi, useNames = TRUE)
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
    x = x@deinflatedcounts,
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
#'@aliases colQuantiles,Zi-method
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
    x = x@deinflatedcounts,
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
#'@returns mean value
#'@seealso \link[stats]{weighted.mean}, \link[zitools]{colMeans2}, \link[zitools]{rowMeans2}
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'mean(Zi)
#'

mean.Zi <- function(x, ...) {
  Zi <- x
  x <- Zi@inputcounts
  w <- Zi@weights
  mean <- weighted.mean(x, w, ...)
  return(mean)
}

#'@export
#'@name colMeans2
#'@title Calculate the row or column means of zero-inflated count data
#'@aliases colMeans2,Zi-method
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
    rows <-  c(1:nrow(x@inputcounts))}
  if(is.null(cols)){
    cols <-  c(1:ncol(x@inputcounts))}
  colmean <-
    mapply(
      weighted.mean,
      as.data.frame(x@inputcounts)[rows,cols],
      as.data.frame(x@weights)[rows, cols],USE.NAMES = useNames, na.rm = na.rm)
  return(colmean)
})

#'@name rowMeans2
#'@aliases rowMeans2,Zi-method
#'@export
#'@rdname colMeans2
#'@importFrom stats weighted.mean
#'@importFrom MatrixGenerics rowMeans2
#'
setMethod("rowMeans2", "Zi", function(x, rows = NULL, cols = NULL, na.rm = FALSE, useNames=TRUE) {
  if(is.null(rows)){
    rows <- c(1:nrow(t(x@inputcounts)))}
  if(is.null(cols)){
    cols <- c(1:ncol(t(x@inputcounts)))}
  rowmean <-
    mapply(weighted.mean,
           as.data.frame(t(x@inputcounts))[rows,cols],
           as.data.frame(t(x@weights))[rows,cols], USE.NAMES = useNames, na.rm=na.rm)
  return(rowmean)
})



#'@name sd
#'@aliases sd,Zi-method
#'@title Standard Deviation of zero inflated count data
#'@param x  A \code{\linkS4class{Zi}}-class object
#'@param na.rm \code{\link[base]{logical}} If \code{\link{TRUE}} \code{\link{NA}}s
#'are excluded, otherwise not. default = \code{\link{FALSE}}
#'
#'@description  Calculate the standard deviation of zero inflated count data
#'taking weights for structural zeros into account.
#'@importFrom matrixStats weightedSd
#'@importFrom BiocGenerics sd
#'@returns standard deviation value
#'@seealso \link[matrixStats]{weightedSd}, \link[zitools]{rowSds}, \link[zitools]{colSds}
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'sd(Zi)


setMethod("sd", "Zi", function(x, na.rm = FALSE) {
  sd <- matrixStats::weightedSd(x = x@inputcounts,
                                w = x@weights, na.rm = na.rm)
  return(sd)
})

#'@export
#'@name rowSds
#'@aliases rowSds,Zi-method
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
    rows <- c(1:nrow(t(x@inputcounts)))}
  if(is.null(cols)){
    cols <- c(1:ncol(t(x@inputcounts)))}
  mapply(weightedSd,
         as.data.frame(t(x@inputcounts))[rows,cols],
         as.data.frame(t(x@weights))[rows,cols], USE.NAMES = useNames, na.rm = na.rm)
})

#'@name colSds
#'@aliases colSds,Zi-method
#'@export
#'@rdname rowSds
#'@importFrom matrixStats weightedSd
#'@importFrom MatrixGenerics colSds

setMethod("colSds", "Zi", function(x, rows = NULL, cols = NULL, na.rm = FALSE, useNames = TRUE) {
  if(is.null(rows)){
    rows <-  c(1:nrow(x@inputcounts))}
  if(is.null(cols)){
    cols <-  c(1:ncol(x@inputcounts))}
  mapply(
    weightedSd,
    as.data.frame(x@inputcounts)[rows,cols],
    as.data.frame(x@weights)[rows,cols], USE.NAMES = useNames, na.rm = na.rm)
})

#'@export
#'@name var
#'@aliases var,Zi-method
#'@title Variance of zero inflated count data
#'@param x  A \code{\linkS4class{Zi}}-class object
#'@param na.rm \code{\link[base]{logical}} If \code{\link{TRUE}} \code{\link{NA}}s
#'are excluded, otherwise not. default = \code{\link{FALSE}}
#'
#'@description  Calculate the variance of zero inflated count data taking weights
#'for structural zeros into account.
#'
#'@returns variance value
#'@importFrom matrixStats weightedVar
#'@importFrom BiocGenerics var
#'@seealso \link[matrixStats]{weightedVar}, \link[zitools]{rowVars}, \link[zitools]{colVars}
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'var(Zi)

setMethod("var", "Zi", function(x, na.rm = FALSE) {
  var <- weightedVar(x = x@inputcounts,
                                  w = x@weights, na.rm = na.rm)
  return(var)
})

#'@export
#'@name rowVars
#'@aliases rowVars,Zi-method
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
    rows <- c(1:nrow(t(x@inputcounts)))}
  if(is.null(cols)){
    cols <- c(1:ncol(t(x@inputcounts)))}
  mapply(weightedVar,
         as.data.frame(t(x@inputcounts))[rows,cols],
         as.data.frame(t(x@weights))[rows,cols], USE.NAMES = useNames, na.rm = na.rm)
})

#'@name colVars
#'@aliases colVars,Zi-method
#'@export
#'@rdname rowVars
#'@importFrom matrixStats weightedVar
#'@importFrom MatrixGenerics colVars

setMethod("colVars", "Zi", function(x, rows = NULL, cols = NULL, na.rm = FALSE, useNames = TRUE) {
  if(is.null(rows)){
    rows <-  c(1:nrow(x@inputcounts))}
  if(is.null(cols)){
    cols <-  c(1:ncol(x@inputcounts))}
  mapply(
    weightedVar,
    as.data.frame(x@inputcounts)[rows,cols],
    as.data.frame(x@weights)[rows,cols], USE.NAMES = useNames, na.rm = na.rm)
})

#'@export
#'@name weighted.mean
#'@aliases weighted.mean,Zi-method
#'@title Weighted Arithmetic Mean of zero inflated count data
#'
#'@param x A \code{\linkS4class{Zi}}-class object
#'@param w a numerical \code{\link[base]{vector}} of weight the same length as
#'x giving the weights to use for elements of x
#'@param ... \link[stats]{weighted.mean}
#'
#'@description Calculate a weighted mean of zero inflated count data, additionally
#' taking weights for structural zeros into account
#'@returns weighted mean value
#'@importFrom stats weighted.mean
#'@seealso \link[stats]{weighted.mean}, \link[zitools]{rowWeightedMeans}, \link[zitools]{colWeightedMeans}
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'weight <- runif(length(Zi@inputcounts), 0.1, 1)
#'weighted.mean(Zi, w= weight)


setMethod("weighted.mean", "Zi", function(x, w, ...) {
  mean <- weighted.mean(x@inputcounts, w = w * x@weights, ...)
  return(mean)
})

#'@export
#'@name rowWeightedMeans
#'@aliases rowWeightedMeans,Zi-method
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
#'rowWeightedMeans(Zi, w = runif(ncol(Zi@inputcounts), 0.1,1))
#'colWeightedMeans(Zi, w = runif(nrow(Zi@inputcounts), 0.1,1))
#'

setMethod("rowWeightedMeans", "Zi", function(x,
                                             w,
                                             rows = NULL,
                                             cols = NULL,
                                             na.rm = FALSE,
                                             useNames = TRUE) {
  if (is.null(rows)) {
    rows <- c(1:nrow(t(x@inputcounts)))
  }
  if (is.null(cols)) {
    cols <- c(1:ncol(t(x@inputcounts)))
  }
  mapply(
    weighted.mean,
    as.data.frame(t(x@inputcounts))[rows, cols],
    as.data.frame(t(x@weights * w))[rows, cols],
    USE.NAMES = useNames,
    na.rm = na.rm
  )
})


#'@name colWeightedMeans
#'@aliases colWeightedMeans,Zi-method
#'@export
#'@rdname rowWeightedMeans
#'@importFrom stats weighted.mean
#'@importFrom MatrixGenerics colWeightedMeans

setMethod("colWeightedMeans", "Zi", function(x, w, rows = NULL, cols = NULL, na.rm = FALSE, useNames = TRUE) {
  if(is.null(rows)){
    rows <-  c(1:nrow(x@inputcounts))}
  if(is.null(cols)){
    cols <-  c(1:ncol(x@inputcounts))}
  mapply(
    weighted.mean,
    as.data.frame(x@inputcounts)[rows,cols],
    as.data.frame(x@weights*w)[rows,cols], USE.NAMES = useNames, na.rm = na.rm)
})

#'@export
#'@name weightedSd
#'@aliases weightedSd,Zi-method
#'@rdname weightedVar
#'@importFrom matrixStats weightedSd

setGeneric("weightedSd", function(x, w = NULL, idxs = NULL, na.rm = FALSE, center = NULL,
                                   ...) standardGeneric("weightedSd"))

setMethod("weightedSd", "Zi", function(x, w, idxs = NULL, na.rm = FALSE, center = NULL, ...) {
  sqrt(weightedVar(x = x@inputcounts,
              w = w * x@weights, idxs = idxs, na.rm = na.rm, center = NULL, ...))
})


#'@export
#'@name rowWeightedSds
#'@aliases rowWeightedSds,Zi-method
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
#'rowWeightedSds(Zi, w = runif(ncol(Zi@inputcounts), 0.1,1))
#'colWeightedSds(Zi, w = runif(nrow(Zi@inputcounts), 0.1,1))
#'rowWeightedVars(Zi, w = runif(ncol(Zi@inputcounts), 0.1,1))
#'colWeightedVars(Zi, w = runif(nrow(Zi@inputcounts), 0.1,1))
#'

setMethod("rowWeightedSds", "Zi", function(x,
                                             w,
                                             rows = NULL,
                                             cols = NULL,
                                             na.rm = FALSE,
                                             useNames = TRUE) {
  if (is.null(rows)) {
    rows <- c(1:nrow(t(x@inputcounts)))
  }
  if (is.null(cols)) {
    cols <- c(1:ncol(t(x@inputcounts)))
  }
  mapply(
    weightedSd,
    as.data.frame(t(x@inputcounts))[rows, cols],
    as.data.frame(t(x@weights * w))[rows, cols],
    USE.NAMES = useNames,
    na.rm = na.rm
  )
})

#'@name colWeightedSds
#'@aliases colWeightedSds,Zi-method
#'@export
#'@rdname rowWeightedSds
#'@importFrom matrixStats weightedSd
#'@importFrom MatrixGenerics colWeightedSds

setMethod("colWeightedSds", "Zi", function(x, w, rows = NULL, cols = NULL, na.rm = FALSE, useNames = TRUE) {
  if(is.null(rows)){
    rows <-  c(1:nrow(x@inputcounts))}
  if(is.null(cols)){
    cols <-  c(1:ncol(x@inputcounts))}
  mapply(
    weightedSd,
    as.data.frame(x@inputcounts)[rows,cols],
    as.data.frame(x@weights*w)[rows,cols], USE.NAMES = useNames, na.rm = na.rm)
})

#'@export
#'@name weightedVar
#'@aliases weightedVar,Zi-method
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
#' @param ... \code{\link[matrixStats]{weightedVar}}
#'@description  Calculate a weighted variance and standard deviation of zero inflated count data,
#'additionally taking weights for structural zeros into account
#'@returns a \code{\link{numeric}} scalar
#'@importFrom matrixStats weightedVar
#'@seealso \link[matrixStats]{weightedVar}, \link[zitools]{rowWeightedVars}, \link[zitools]{colWeightedVars}
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'weight <- runif(length(Zi@inputcounts), 0.1, 1)
#'weightedVar(Zi, w= weight)
#'weightedSd(Zi, w = weight)

setGeneric("weightedVar", function(x, w = NULL, idxs = NULL, na.rm = FALSE, center = NULL,
                                   ...) standardGeneric("weightedVar"))

setMethod("weightedVar", "Zi", function(x, w, idxs = NULL, na.rm = FALSE, center = NULL, ...) {
  weightedVar(x = x@inputcounts,
              w = w * x@weights, idxs = idxs, na.rm = na.rm, center = NULL, ...)
})

#'@export
#'@name rowWeightedVars
#'@aliases rowWeightedVars,Zi-method
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
    rows <- c(1:nrow(t(x@inputcounts)))
  }
  if (is.null(cols)) {
    cols <- c(1:ncol(t(x@inputcounts)))
  }
  mapply(
    weightedVar,
    as.data.frame(t(x@inputcounts))[rows, cols],
    as.data.frame(t(x@weights * w))[rows, cols],
    USE.NAMES = useNames,
    na.rm = na.rm
  )
})

#'@name colWeightedVars
#'@aliases colWeightedVars,Zi-method
#'@export
#'@rdname rowWeightedSds
#'@importFrom matrixStats weightedVar
#'@importFrom MatrixGenerics colWeightedVars


setMethod("colWeightedVars", "Zi", function(x, w, rows = NULL, cols = NULL, na.rm = FALSE, useNames = TRUE) {
  if(is.null(rows)){
    rows <-  c(1:nrow(x@inputcounts))}
  if(is.null(cols)){
    cols <-  c(1:ncol(x@inputcounts))}
  mapply(
    weightedVar,
    as.data.frame(x@inputcounts)[rows,cols],
    as.data.frame(x@weights*w)[rows,cols], USE.NAMES = useNames, na.rm = na.rm)
})


#'@name log1p
#'@aliases log1p,Zi-method
#'@title log(1+x)
#'@description Calculate log(1+x) of all 'matrix' objects of a 'Zi'-class
#'object, log calculates by default natural logarithms
#'@param x \code{\linkS4class{Zi}}-class object
#'@export
#'@seealso \link[base]{log1p}, \link[zitools]{log2p}
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'log1p(Zi)

setMethod("log1p", signature ="Zi", definition = function(x){
  inputcounts <- log1p(x@inputcounts)
  deinflatedcounts <- log1p(x@deinflatedcounts)
  weights <- log1p(x@weights)
  result <- new(
    Class = "Zi",
    inputdata = x@inputdata,
    inputcounts = inputcounts,
    model = x@model,
    deinflatedcounts = deinflatedcounts,
    weights = weights)
})

#'@name log2p
#'@aliases log2,Zi-method
#'@title log2(x+1)
#'@description Calculate log2(x+1) of all 'matrix' objects of a 'Zi'-class
#'object
#'@param x \code{\linkS4class{Zi}}-class object
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'log2p(Zi)
#'
setGeneric("log2p", function(x){
  log2(x+1)
})

setMethod("log2p", signature ="Zi",definition = function(x){
inputcounts <- log2p(x@inputcounts)
deinflatedcounts <- log2p(x@deinflatedcounts)
weights <- log2p(x@weights)
result <- new(
  Class = "Zi",
  inputdata = x@inputdata,
  inputcounts = inputcounts,
  model = x@model,
  deinflatedcounts = deinflatedcounts,
  weights = weights)
return(result)
})


#'@name +
#'@aliases +,Zi-method
#'@title Arithmetic Operators
#'@description Arithmetic operators for a Zi-class object
#'@param e1 \code{\linkS4class{Zi}}-class object, matrix or number
#'@param e2 \code{\linkS4class{Zi}}-class object, matrix or number
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'Zi+Zi
#'Zi+2
#'
setMethod("+", signature = "Zi", definition = function(e1,e2){
  if(class(e1)=="Zi")
    values1 <- e1@inputcounts
  else
    values1 <- e1
  if(class(e2)=="Zi")
    values2 <- e2@inputcounts
  else
    values2 <- e2
  adinputcounts <- values1 + values2

  if(class(e1)=="Zi")
    values1 <- e1@deinflatedcounts
  else
    values1 <- e1
  if(class(e2)=="Zi")
    values2 <- e2@deinflatedcounts
  else
    values2 <- e2
  addeinflatedcounts <- values1+values2

  if(class(e1)=="Zi")
    values1 <- e1@weights
  else
    values1 <- e1
  if(class(e2)=="Zi")
    values2 <- e2@weights
  else
    values2 <- e2
  adweights <- pmin(values1,values2) # currently the point-wise min is used, no error propagation

  result <- new(
    Class = "Zi",
    inputdata = e1@inputdata,
    inputcounts = adinputcounts,
    model = e1@model,
    deinflatedcounts = addeinflatedcounts,
    weights = adweights)
  return(result)
})

#'@name -
#'@aliases -,Zi-method
#'@title Arithmetic Operators
#'@description Arithmetic operators for a Zi-class object
#'@param e1 \code{\linkS4class{Zi}}-class object, matrix or number
#'@param e2 \code{\linkS4class{Zi}}-class object, matrix or number
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'Zi-Zi
#'Zi-2
#'
#'
setMethod("-", signature = "Zi", definition = function(e1,e2){
  if(class(e1)=="Zi")
    values1 <- e1@inputcounts
  else
    values1 <- e1
  if(class(e2)=="Zi")
    values2 <- e2@inputcounts
  else
    values2 <- e2
  adinputcounts <- values1-values2
  if(sum(adinputcounts<0,na.rm=T)>0)
    warning("Negative counts produced.")

  if(class(e1)=="Zi")
    values1 <- e1@deinflatedcounts
  else
    values1 <- e1
  if(class(e2)=="Zi")
    values2 <- e2@deinflatedcounts
  else
    values2 <- e2
  addeinflatedcounts <- values1-values2

  if(class(e1)=="Zi")
    values1 <- e1@weights
  else
    values1 <- e1
  if(class(e2)=="Zi")
    values2 <- e2@weights
  else
    values2 <- e2
  adweights <- pmin(values1,values2) # currently the point-wise min is used, no error propagation

  result <- new(
    Class = "Zi",
    inputdata = e1@inputdata,
    inputcounts = adinputcounts,
    model = e1@model,
    deinflatedcounts = addeinflatedcounts,
    weights = adweights)
  return(result)
})

#'@name *
#'@aliases *,Zi-method
#'@title Arithmetic Operators
#'@description Arithmetic operators for a Zi-class object
#'@param e1 \code{\linkS4class{Zi}}-class object, matrix or number
#'@param e2 \code{\linkS4class{Zi}}-class object, matrix or number
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'Zi*Zi
#'Zi*2
#'
setMethod("*", signature = "Zi", definition = function(e1,e2){
  if(class(e1)=="Zi")
    values1 <- e1@inputcounts
  else
    values1 <- e1
  if(class(e2)=="Zi")
    values2 <- e2@inputcounts
  else
    values2 <- e2
  adinputcounts <- values1*values2

  if(class(e1)=="Zi")
    values1 <- e1@deinflatedcounts
  else
    values1 <- e1
  if(class(e2)=="Zi")
    values2 <- e2@deinflatedcounts
  else
    values2 <- e2
  addeinflatedcounts <- values1*values2

  if(class(e1)=="Zi")
    values1 <- e1@weights
  else
    values1 <- e1
  if(class(e2)=="Zi")
    values2 <- e2@weights
  else
    values2 <- e2
  adweights <- pmin(values1,values2) # currently the point-wise min is used, no error propagation

  result <- new(
    Class = "Zi",
    inputdata = e1@inputdata,
    inputcounts = adinputcounts,
    model = e1@model,
    deinflatedcounts = addeinflatedcounts,
    weights = adweights)
  return(result)
})

#'@name /
#'@aliases /,Zi-method
#'@title Arithmetic Operators
#'@description Arithmetic operators for a Zi-class object
#'@param e1 \code{\linkS4class{Zi}}-class object, matrix or number
#'@param e2 \code{\linkS4class{Zi}}-class object, matrix or number
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'Zi/Zi
#'Zi/2
#'
#'
setMethod("/", signature = "Zi", definition = function(e1,e2){
  if(class(e1)=="Zi")
    values1 <- e1@inputcounts
  else
    values1 <- e1
  if(class(e2)=="Zi")
    values2 <- e2@inputcounts
  else
    values2 <- e2
  adinputcounts <- values1/values2

  if(class(e1)=="Zi")
    values1 <- e1@deinflatedcounts
  else
    values1 <- e1
  if(class(e2)=="Zi")
    values2 <- e2@deinflatedcounts
  else
    values2 <- e2
  addeinflatedcounts <- values1/values2

  if(class(e1)=="Zi")
    values1 <- e1@weights
  else
    values1 <- e1
  if(class(e2)=="Zi")
    values2 <- e2@weights
  else
    values2 <- e2
  adweights <- pmin(values1,values2) # currently the point-wise min is used, no error propagation

  result <- new(
    Class = "Zi",
    inputdata = e1@inputdata,
    inputcounts = adinputcounts,
    model = e1@model,
    deinflatedcounts = addeinflatedcounts,
    weights = adweights)
  return(result)
})




#'@export
#'@name show
#'@aliases show,Zi-method
#'@title Show summary of Zi object
#'@description Message printed at command line
#'@param  object        \code{\linkS4class{Zi}}-class object
#'
#'@returns returns a numeric vector of row/column length
#'@importFrom methods show
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'Zi
#'show(Zi)
#'
setMethod("show", "Zi" , function(object) {
  dims <- dim(object@inputdata)
  formel <- as.character(Zi@model[[1]]$formula)
  anz <- length(object@inputcounts) # data points
  anz0 <- sum(object@inputcounts==0,na.rm=T) # zeros
  anz00 <- sum(is.na(object@deinflatedcounts)) # structural zeros

  cat("Formal class 'Zi' [package \"zitools\"]\n")
  cat(" ",dims[1],"features (rows),",dims[2],"samples (columns)\n")
  cat(paste0(" ",anz," data points, ",
      anz0, "(",round(anz0/anz*100,3),"%) zeros, ",
      anz00,"(",round(anz00/anz*100,3),"%) structual zeros estimated with ",
      formel[2]," ",formel[1]," ",formel[3],"\n"))
  str(object, list.len = 10,max.level = 2)
  cat("Use str(object) to inspect the whole object structure.")
})

#'@name tax_table
#'@title Access the taxonomy table
#'@aliases tax_table,Zi-method
#'@param object \code{\linkS4class{Zi}}-class object
#'@description access the taxonomy table (\link[phyloseq]{tax_table}) of an \code{\linkS4class{Zi}}-class object if the
#'inputdata slot is a phyloseq object
#'@importFrom phyloseq tax_table
#'@returns tax_table
#'@export


setMethod("tax_table", signature = "Zi", function(object){
  if("phyloseq" %in% class(object@inputdata)) {
    tax_table <- tax_table(object@inputdata)}
  if("matrix" %in% class(object@inputdata)){
    tax_table <- NULL}
  if("SummarizedExperiment" %in% class(object@inputdata)){
    tax_table <- NULL}
  return(tax_table)
})

#'@name sample_data
#'@aliases sample_data,Zi-method
#'@title Access the sample data
#'@param object \code{\linkS4class{Zi}}-class object
#'@description access the \link[phyloseq]{sample_data} of an \code{\linkS4class{Zi}}-class object if the
#'inputdata slot is a phyloseq object
#'@importFrom phyloseq sample_data
#'@returns sample_data
#'
#'@export
#'

setMethod("sample_data", signature = "Zi", function(object){
  if("phyloseq" %in% class(object@inputdata)) {
    sample_data <- sample_data(object@inputdata)}
  if("matrix" %in% class(object@inputdata)){
    sample_data <- NULL}
  if("SummarizedExperiment" %in% class(object@inputdata)){
    sample_data <- NULL}
  return(sample_data)
})

#'@name otu_table
#'@aliases otu_table,Zi-method
#'@title Access the otu table
#'@param object \code{\linkS4class{Zi}}-class object
#'@description access the \link[phyloseq]{otu_table} of an \code{\linkS4class{Zi}}-class object if the
#'inputdata slot is a phyloseq object
#'@returns otu_table
#'@importFrom phyloseq otu_table
#'
#'@export
#'

setMethod("otu_table", signature = "Zi", function(object){
  if("phyloseq" %in% class(object@inputdata)) {
    otu_table <- otu_table(object@inputdata)}
  if("matrix" %in% class(object@inputdata)){
    otu_table <- NULL}
  if("SummarizedExperiment" %in% class(object@inputdata)){
    otu_table <- NULL}
  return(otu_table)
})

#'@name phy_tree
#'@aliases phy_tree,Zi-method
#'@title Access the phylogenetic tree
#'@param physeq \code{\linkS4class{Zi}}-class object
#'@description access the phylogenetic tree (\link[phyloseq]{phy_tree}) of an object of the class
#'"Zi" if the inputdata slot is a phyloseq object
#'@returns phy_tree
#'@importFrom phyloseq phy_tree
#'
#'@export
#'

setMethod("phy_tree", signature = "Zi", function(physeq){
  if("phyloseq" %in% class(physeq@inputdata)) {
    phy_tree <- phy_tree(physeq@inputdata)}
  if("matrix" %in% class(physeq@inputdata)){
    phy_tree <- NULL}
  if("SummarizedExperiment" %in% class(physeq@inputdata)){
    phy_tree <- NULL}
  return(phy_tree)
})

#'@name rowData
#'@aliases rowData,Zi-method
#'@title Access the row data
#'@param x \code{\linkS4class{Zi}}-class object
#'@param useNames returns a rowData dataframe with rownames
#'@param ... \code{\link[SummarizedExperiment]{rowData}}
#'@description access the \link[SummarizedExperiment]{rowData} of an \code{\linkS4class{Zi}}-class object if the inputdata
#'is an object of the class SummarizedExperiment
#'@returns DFrame
#'@importFrom SummarizedExperiment rowData
#'
#'@export
#'

setMethod("rowData", signature = "Zi", function(x, useNames = TRUE, ...){
  if("SummarizedExperiment" %in% class(x@inputdata)) {
    rowData <- rowData(x@inputdata, useNames = useNames, ...)}
  if("matrix" %in% class(x@inputdata)){
    rowData <- NULL}
  if("phyloseq" %in% class(x@inputdata)){
    rowData <- NULL}
  return(rowData)
})

#'@name assays
#'@aliases assays,Zi-method
#'@title Access assays
#'@param x \code{\linkS4class{Zi}}-class object
#'@param withDimnames A \code{logical}, indicating whether the dimnames of the
#'SummarizedExperiment object should be applied (i.e. copied) to the extracted
#'assays. see \code{\link[SummarizedExperiment]{assays}}
#'@param ... see \code{\link[SummarizedExperiment]{assays}}
#'@description access  \link[SummarizedExperiment]{assays} of an \code{\linkS4class{Zi}}-class object if the inputdata is an
#'object of the class SummarizedExperiment
#'@importFrom SummarizedExperiment assays
#'@returns list
#'@export


setMethod("assays", signature = "Zi", function(x, withDimnames=TRUE,  ...){
  if("SummarizedExperiment" %in% class(x@inputdata)) {
    assays <- assays(x@inputdata, withDimnames=withDimnames, ...)}
  if("matrix" %in% class(x@inputdata)){
    assays <- NULL}
  if("phyloseq" %in% class(x@inputdata)){
    assays <- NULL}
  return(assays)
})


#'@name colData
#'@aliases colData,Zi-method
#'@title Access the col Data
#'@param x \code{\linkS4class{Zi}}-class object
#'@param ... \code{\link[SummarizedExperiment]{colData}}
#'@description access the \link[SummarizedExperiment]{colData} of an \code{\linkS4class{Zi}}-class object if the inputdata is an
#'object of the class SummarizedExperiment
#'@importFrom SummarizedExperiment colData
#'@returns DFrame
#'@seealso \code{\link[SummarizedExperiment]{colData}}
#'
#'@export

setMethod("colData", signature = "Zi", function(x, ...){
  if("SummarizedExperiment" %in% class(x@inputdata)) {
    colData <- colData(x@inputdata, ...)}
  if("matrix" %in% class(x@inputdata)){
    colData <- NULL}
  if("phyloseq" %in% class(x@inputdata)){
    colData <- NULL}
  return(colData)
})

#'@name t
#'@aliases t,Zi-method
#'@title Transpose a \code{\linkS4class{Zi}}-class object
#'@param x \code{\linkS4class{Zi}}-class object
#'@description transpose all matrizes of a \code{\linkS4class{Zi}}-class object
#'@returns \code{\linkS4class{Zi}}-class object
#'@export
#'

setMethod("t", signature = "Zi", definition = function(x){
  inputcounts <- t(x@inputcounts)
  deinflatedcounts <- t(x@deinflatedcounts)
  weights <- t(x@weights)
  result <- new(
    Class = "Zi",
    inputdata = x@inputdata,
    inputcounts = inputcounts,
    model = x@model,
    deinflatedcounts = deinflatedcounts,
    weights = weights
  )
  return(result)
})
