#'@name boxplot
#'
#'@title Box Plots
#'
#'@param ZiObject ZiObject, result of the ziMain function
#'@param log1p logical, default = FALSE, if TRUE log(1+p) transformation takes
#'place
#'@param ... see boxplot documentation
#'
#'@description produce box-and-whisker plot(s) of the given (grouped) values.
#'boxplot.Zi uses the output matrix (drawn structural zeros) to produce box-and
#'whisker plots
#'
#'@returns boxplot
#'@example

boxplot.Zi <- function(result_zi, log1p =FALSE, ...)
{
  if (log1p == TRUE) {
    boxplot(log1p(result_zi@output), ...)
  }
  if (log1p == FALSE) {
    boxplot(result_zi@output)
  }
}

#'@name heatmap
#'@title Draw a Heat Map
#'
#'@param ziObject
#'@description draw a heatmap of a given ziObject, heatmap.Zi uses the output
#'matrix (drawn structural zeros) to produce a heatmap. NA values are white
#'
#'@returns heatmap
#'
#'@example
#'
#'
#'

heatmap <- function(input, ...) {
  UseMethod("heatmap")
}

heatmap.Zi <- function(result_zi, ...) {
  df <- as.data.frame(result_zi@output)
  df <- df %>%
    filter_all(any_vars(!is.na(.)))%>%
    filter_all(any_vars(. != 0))
  mtx <- as.matrix(df)
  stats::heatmap(mtx, ...)
}
#'@name weightedCor
#'
#'@title weighted Correlation
#'
#'@description Calculate the weighted pearson correlation coefficients of a
#'Zi object. The inputmatrix is used to calculate column wise correlations taking
#'the weights matrix of the Zi object into account.
#'
#'@param x Zi Object
#'@param y another Zi Object (same dimension as x)
#'@param na.rm If TRUE (default), missing values are excluded
#'@param transpose If TRUE, rowwise correlations are calculated, default = FALSE
#'
#'@returns a matrix of weighted pearson correlation coefficients
#'@export
#'@example
#'
#'

weightedCor <- function(x, y = NULL, na.rm=TRUE, transpose = FALSE, ...) {
  my_vector <- numeric()
  wx <- x@weights
  x <- x@inputmatrix
  if (is.null(y)) {
    y <- x
    wy <- wx
  }
  if(transpose == TRUE) {
    x <- t(x)
    wx <- t(wx)
    y <- t(y)
    wy <- t(wy)
  }
  colnames <- colnames(x)
  rownames <- colnames(y)
  for (a in 1:ncol(x)) {
    for (b in 1:ncol(x)) {
      col_a <- x[,a]
      col_b <- y[,b]
      weights_a <- wx[,a]
      weights_b <- wy[,b]
      mean_a <- sum(weights_a*col_a, na.rm = na.rm)/(sum(weights_a, na.rm=na.rm))
      mean_b <- sum(weights_b*col_b, na.rm = na.rm)/(sum(weights_b, na.rm = na.rm))
      var_a <- sum(weights_a*(col_a-mean_a)^2, na.rm = na.rm)/(sum(weights_a, na.rm = na.rm)-1)
      var_b <- sum(weights_b*(col_b-mean_b)^2, na.rm = na.rm)/(sum(weights_b, na.rm = na.rm)-1)
      cov <-
        sum(sqrt(weights_a)*(col_a - mean_a) * sqrt(weights_b)*(col_b - mean_b), na.rm = na.rm) / sqrt( (sum(weights_a, na.rm = na.rm)-1) * (sum(weights_b, na.rm = na.rm)-1))
      cor <- cov / sqrt(var_a * var_b)
      my_vector <- c(my_vector, cor)
    }
  }
  mtx <- matrix(my_vector, ncol(x))
  colnames(mtx) <- colnames
  rownames(mtx) <- rownames
  return(mtx)
}
