#'@name boxplot
#'@title Create boxplots of a 'Zi'-class object
#'@description Create boxplots of a 'Zi'-class object.
#'@param x 'Zi'-class object
#'@param log1p logical, default = FALSE, if TRUE log(1+p) transformation takes
#'place
#'@param ... see boxplot documentation
#'@seealso \link[graphics]{boxplot}
#'@importFrom graphics boxplot
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

heatmap <- function(x, ...) {
  UseMethod("heatmap")
}

#'@name heatmap
#'@title Draw a Heat Map
#'@param x 'Zi'-class object
#'@description draw a heatmap of a given 'Zi'-class object, heatmap.Zi uses the output
#'matrix (drawn structural zeros) to produce a heatmap. NA values are white
#'@returns heatmap
#'@export
#'@example
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




setGeneric("cor", function(x, y = NULL, use = "everything",
                           method = c("pearson", "kendall", "spearman")) standardGeneric("cor"))

#'@name cor
#'@title Calculate weighted Pearson Correlation coeffiecients
#'@description calculate the weighted pearson correlation coefficients of a count matrix
#'of an Zi object taking weights for structural zeros into account
#'@param x 'Zi'-class object
#'@importFrom stats cor
#'@export
setMethod("cor", signature = "Zi", function(x, y = NULL, use = "everything", method = "pearson"){
  my_vector <- numeric()
  wx <- x@weights
  cx <- x@countmatrix
  if (is.null(y)) {
    y <- cx
    wy <- wx
  }
  colnames <- colnames(cx)
  rownames <- colnames(y)
  for (a in 1:ncol(cx)) {
    for (b in 1:ncol(cx)) {
      col_a <- cx[,a]
      col_b <- y[,b]
      weights_a <- wx[,a]
      weights_b <- wy[,b]
      mean_a <- sum(weights_a*col_a)/(sum(weights_a))
      mean_b <- sum(weights_b*col_b)/(sum(weights_b))
      var_a <- sum(weights_a*(col_a-mean_a)^2)/(sum(weights_a)-1)
      var_b <- sum(weights_b*(col_b-mean_b)^2)/(sum(weights_b)-1)
      cov <-
        sum(sqrt(weights_a)*(col_a - mean_a) * sqrt(weights_b)*(col_b - mean_b)) / sqrt( (sum(weights_a)-1) * (sum(weights_b)-1))
      cor <- cov / sqrt(var_a * var_b)
      my_vector <- c(my_vector, cor)
    }
  }
  mtx <- matrix(my_vector, ncol(cx))
  colnames(mtx) <- colnames
  rownames(mtx) <- rownames
  print(message("Currently implemented only for Pearson Correlation"))
  return(mtx)
})

