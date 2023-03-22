#'@name boxplot
#'
#'@title Box Plots
#'
#'@param x ZiObject, result of the ziMain function
#'@param log1p logical, default = FALSE, if TRUE log(1+p) transformation takes
#'place
#'@param ... see boxplot documentation
#'@export
#'@example

boxplot.Zi <- function(x, log1p =FALSE, ...)
{
  if (log1p == TRUE) {
    boxplot(log1p(x@output), ...)
  }
  if (log1p == FALSE) {
    boxplot(x@output)
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
?heatmap
heatmap <- function(x, ...) {
  UseMethod("heatmap")
}

heatmap.Zi <- function(x, ...) {
  df <- as.data.frame(x@output)
  df <- df %>%
    filter_all(any_vars(!is.na(.)))%>%
    filter_all(any_vars(. != 0))
  mtx <- as.matrix(df)
  stats::heatmap(mtx, ...)
}

#'@name MissingValueHeatmap
#'@title Missing Value Heatmap
#'@description
#'
#'@param ZiObject ZiObject, result of the ziMain function
#'
#'@returns heatmap
#'
#'
#'
MissingValueHeatmap <- function(ZiObject) {
  mtx <- ZiObject@output
  # Sort matrix according to number of missing values
  mtx.heatmap <- mtx[order(-rowSums(is.na(mtx))), ]
  stats::heatmap(mtx.heatmap,Rowv=NA,Colv=NA,labRow=FALSE,
                 na.rm=T,col=RColorBrewer::brewer.pal(n = 9, name = "Blues"),
                 scale = "none", margins = c(11,0), cexCol=1)
}
#'@name weightedCor
#'
#'@title weighted Correlation
#'
#'@description calculate the weighted pearson correlation coefficients of a
#'matrix x (and y) taking a weight matrix for x (and y) into account
#'
#'@param x matrix
#'@param wx weight matrix (same dimension as x)
#'@param y matrix
#'@param wy weight matrix (same dimension as y)
#'@param na.rm If TRUE (default), missing values are excluded
#'
#'@returns a matrix of weighted pearson correlation coefficients
#'
#'@example
#'mtx <- matrix(runif(1000,0,1000), 10, 100)
#'w <- matrix(runif(1000, 0.01,1), 10, 100)
#'msb.WeightedCor(x = mtx, wx = w)
#'

setGeneric("weightedCor", function(x, wx, y = NULL, wy = NULL, na.rm=TRUE, ...) {
  my_vector <- numeric()
  if (is.null(y)) {
    y <- x
  }
  if(is.null(wy)) {
    wy <- wx
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
})


setMethod("weightedCor", signature = "Zi", definition = function(x, y = NULL, na.rm=TRUE, transpose = FALSE, ...) {
  my_vector <- numeric()
  wx <- x@weights
  cx <- x@countmatrix
  if (is.null(y)) {
    y <- cx
    wy <- wx
  }
  mtx <- weightedCor(x = cx, wx = wx, y=y, wy =wy, na.rm=na.rm, ...)
  return(mtx)
})
