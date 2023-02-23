#median - total, row, column
median.zi <- function(result_zi, na.rm =TRUE, ...)
{
  median(result_zi$ziOutput, na.rm = na.rm,  ...)
}

colMedians <- function(input, ...)
{
  UseMethod("colMedians")
}

colMedians.zi <- function(result_zi, na.rm = TRUE, ...)
{
  MatrixGenerics::colMedians(result_zi$ziOutput, na.rm = na.rm, ...)
}

rowMedians <- function(input, ...)
{
  UseMethod("rowMedians")
}
rowMedians.zi <- function(result_zi, na.rm = TRUE,  ...)
{
  MatrixGenerics::rowMedians(result_zi$ziOutput, na.rm = na.rm,  ...)
}

#quantiles - total, row, column
quantile.zi <- function(result_zi, na.rm = TRUE, ...)
{
  quantile(result_zi$ziOutput, na.rm = na.rm, ...)
}


rowQuantiles <- function(input, ...)
{
  UseMethod("rowQuantiles")
}
rowQuantiles.zi <- function(result_zi, na.rm = TRUE, ...)
{
  MatrixGenerics::rowQuantiles(result_zi$ziOutput, na.rm = na.rm, ...)
}

colQuantiles <- function(input, ...)
{
  UseMethod("colQuantiles")
}
colQuantiles.zi <- function(result_zi, na.rm = TRUE, ...)
{
  MatrixGenerics::colQuantiles(result_zi$ziOutput, na.rm = na.rm, ...)
}

weighted.mean.zi <- function(zi_result, ...){
  x <- zi2inputMatrix(zi_result)
  w <- zi_result[["weights"]]
  mean <- weighted.mean(x, w, ...)
  return(mean)
}

weightedMedian <- function(input, ...) {
  UseMethod("weightedMedian")
}

weightedMedian.zi <- function(zi_result, ...) {
  x <- zi2inputMatrix(zi_result)
  w <- zi_result[["weights"]]
  median <- weightedMedian(x, w,  ...)
  return(median)
}

weightedMean <- function(input, ...){
  UseMethod("weightedMean")
}

weightedMean.zi <- function(zi_result, ...){
  x <- zi2inputMatrix(zi_result)
  w <- zi_result[["weights"]]
  mean <- weighted.mean(x, w, ...)
  return(mean)
}


# function only for NA matrix
mean.zi <- function(zi_result, na.rm=TRUE, ...) {
  x <- zi_result[["ziOutput"]]
  m <- mean(x, na.rm =na.rm, ...)
  return(m)
}
