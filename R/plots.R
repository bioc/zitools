#'@name boxplot
#'@title Create boxplots of a 'Zi'-class object
#'@description Create boxplots of a 'Zi'-class object.
#'@param x 'Zi'-class object
#'@param ... see \link[graphics]{boxplot.default}
#'@seealso \link[graphics]{boxplot.default}
#'@importFrom graphics boxplot
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'boxplot(Zi)
#'boxplot(log1p(Zi))
#'
#'

boxplot.Zi <- function(x, ...)
{
  boxplot(x@output, ...)
}

heatmap <- function(x, ...) {
  UseMethod("heatmap")
}

#'@name heatmap
#'@title Draw a Heat Map
#'@param x 'Zi'-class object
#'@param ... see \link[stats]{heatmap}
#'@description draw a heatmap of a given 'Zi'-class object, heatmap.Zi uses the output
#'matrix (drawn structural zeros) to produce a heatmap. NA values are white
#'@returns heatmap
#'@importFrom dplyr filter_all
#'@importFrom stats heatmap
#'@importFrom dplyr any_vars
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'heatmap(Zi) # error because too many NA
#'heatmap(Zi, Rowv=NA) # no clustering of rows
#'heatmap(Zi, Rowv=NA, Colv=NA) # no clustering of rows and cols
#'
heatmap.Zi <- function(x, ...) {
  df <- as.data.frame(x@output)
  mtx <- as.matrix(df)
  stats::heatmap(mtx, ...)
}

#'@name MissingValueHeatmap
#'@title Missing Value Heatmap
#'@description Missing Value Heatmap
#'
#'@param ZiObject ZiObject, result of the ziMain function
#'
#'@returns heatmap
#'
#'@import ggplot2
#'@importFrom reshape2 melt

MissingValueHeatmap <- function(ZiObject,title = "", xlab = "", ylab = "") {

  mtx <- ZiObject@output
  mtx.heatmap.sorted <- data.frame(mtx[order(-rowSums(is.na(mtx)), rowMeans(mtx, na.rm = TRUE)),])
  mtx.heatmap.sorted$Feature <- row.names(mtx.heatmap.sorted)
  mtx.heatmap.sorted$FeatureIdx <- seq(1, nrow(mtx.heatmap.sorted))

  mtx.long <- reshape2::melt(mtx)
  colnames(mtx.long) <- c("Feature", "Sample", "value")
  heatmap.df <-merge(mtx.long, mtx.heatmap.sorted[, c("Feature", "FeatureIdx")], by="Feature")
  gg.heatmap <-   ggplot(data = heatmap.df, aes(x = Sample, y = FeatureIdx, fill = value))  +
    geom_tile() +
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 9, name = "Blues"), na.value = 'red') +
    labs(x = xlab, y = ylab) +
    ggtitle(title)
  return(gg.heatmap)}

setGeneric("cor", function(x, y = NULL, use = "everything",
                           method = c("pearson", "kendall", "spearman")) standardGeneric("cor"))

#'@name cor
#'@title Calculate weighted Pearson Correlation coeffiecients
#'@description calculate the weighted pearson correlation coefficients of a
#'count matrix of an Zi object taking weights for zero counts into account
#'@param x    'Zi'-class object
#'@param y    'Zi'-class object
#'@param use "everything" see \link[stats]{cor}
#'@param method default = "pearson", weighted correlation only implemented for
#'person correlation
#'@importFrom stats cor
#'@export
#'@examples
#'cor(Zi)
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

