#median - total, row, column
median.Zi <- function(result_zi, na.rm = TRUE, ...)
{
  median(result_zi@output, na.rm = na.rm,  ...)
}

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

#quantiles - total, row, column
quantile.Zi <- function(result_zi, na.rm = TRUE, ...) {
  quantile(result_zi@output, na.rm = na.rm, ...)
}

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


setMethod("colQuantiles", "zi", function(x,
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

#mean - matrix, row, column
mean.Zi <- function(zi_result, ...) {
  x <- zi_result@inputmatrix
  w <- zi_result@weights
  mean <- weighted.mean(x, w, ...)
  return(mean)
}

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

#standard deviation and variance - matrix, row, column
setMethod("sd", "Zi", function(x) {
  sd <- matrixStats::weightedSd(x = x@inputmatrix,
                                w = x@weights)
  return(sd)
})
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

setMethod("var", "Zi", function(x) {
  var <- matrixStats::weightedVar(x = x@inputmatrix,
                                  w = x@weights)
  return(var)
})
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


