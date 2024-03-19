#'@include ziMain.R
NULL


#'@name boxplot
#'@title Create boxplots of a 'Zi'-class object
#'@description Create boxplots of a 'Zi'-class object.
#'@param x 'Zi'-class object
#'@param ... see \link[graphics]{boxplot.default}
#'@seealso \link[graphics]{boxplot.default}
#'@importFrom graphics boxplot
#'@returns  A List with all information to create a boxplot see
#'\link[graphics]{boxplot.default}
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
  boxplot(x@deinflatedcounts, ...)
}
#'@export
#'@importFrom stats heatmap
heatmap <- function(x, ...) {
  UseMethod("heatmap")
}

#'@name heatmap
#'@title Draw a Heat Map
#'@param x 'Zi'-class object
#'@param ... see \link[stats]{heatmap}
#'@description draw a heatmap of a given 'Zi'-class object, heatmap.Zi uses the deinflatedcounts
#'matrix (drawn structural zeros) to produce a heatmap. NA values are white
#'@returns heatmap
#'@importFrom dplyr filter_all
#'@importFrom stats heatmap
#'@importFrom dplyr any_vars
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'#heatmap(Zi) # Error, clustering not possible
#'heatmap(Zi, Rowv=NA) # no clustering of rows
#'heatmap(Zi, Rowv=NA, Colv=NA) # no clustering of rows and cols
#'
heatmap.Zi <- function(x, ...) {
  df <- as.data.frame(x@deinflatedcounts)
  mtx <- as.matrix(df)
  stats::heatmap(mtx, ...)
}

#'@name MissingValueHeatmap
#'@title Missing Value Heatmap
#'@description Missing Value Heatmap
#'
#'@param ZiObject ZiObject, result of the ziMain function
#'@param title Title of the plot .
#'@param ylab Title of the y axis.
#'@param xlab Title of the x axis.
#'
#'@returns heatmap
#'
#'@import ggplot2
#'@importFrom RColorBrewer brewer.pal
#'@importFrom reshape2 melt
#'@examples
#'data(mtx)
#'
#'@export

MissingValueHeatmap <- function(ZiObject,title = "", xlab = "", ylab = "") {

  mtx <- ZiObject@deinflatedcounts
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
#'@aliases cor,Zi,ANY-method
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
#'@returns correlation matrix
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'cor(Zi)
setMethod("cor", c("Zi", "ANY"), function(x, y = NULL, use = "everything", method = "pearson"){
  if(use!="everything")
    stop("zitools::cor only implemented so far for use=\"everything\"")
  if(method!="pearson")
    stop("zitools::cor only implemented so far for Pearson Correlation.")

  my_vector <- numeric()
  wx <- x@weights
  cx <- x@inputcounts
  if (is.null(y)) {
    y <- cx
    wy <- wx
  }
  colnames <- colnames(cx)
  rownames <- colnames(y)
  for (a in seq_len(ncol(cx))) {
    for (b in seq_len(ncol(cx))) {
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
  return(mtx)
})

#'@name cov
#'@aliases cov,Zi,ANY-method
#'@title Calculate weighted Covariance
#'@description calculate the weighted covariance of the columns of the
#'count matrix of an Zi object taking weights for possible structural zero counts into account
#'@param x    'Zi'-class object
#'@param y    'Zi'-class object
#'@param use "everything"
#'@importFrom stats cor
#'@returns covariance matrix
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'cov(Zi)
setMethod("cov", c("Zi","ANY"), function(x, y = NULL, use = "everything"){
  if(use!="everything")
    stop("zitools::cor only implemented so far for use=\"everything\"")
  my_vector <- numeric()
  wx <- x@weights
  cx <- x@inputcounts
  nc <- ncol(cx)
  C <- matrix(nrow = nc,ncol = nc) # empty matrix, correct size

  for (a in (1:(nc-1))) {
    C[a,a] <- 1
    for (b in ((a+1):nc)) {
      col_a <- cx[,a]
      col_b <- cx[,b]
      weights_a <- wx[,a]
      weights_b <- wx[,b]
      mean_a <- sum(weights_a*col_a)/(sum(weights_a))
      mean_b <- sum(weights_b*col_b)/(sum(weights_b))
      var_a <- sum(weights_a*(col_a-mean_a)^2)/(sum(weights_a)-1)
      var_b <- sum(weights_b*(col_b-mean_b)^2)/(sum(weights_b)-1)
      C[a,b] <-
        sum(sqrt(weights_a)*(col_a - mean_a) * sqrt(weights_b)*(col_b - mean_b)) / sqrt( (sum(weights_a)-1) * (sum(weights_b)-1))
      C[b,a] <- C[a,b]
    }
  }
  C[nc,nc] <- 1
  colnames(C) <- colnames(cx)
  rownames(C) <- colnames(cx)
  return(C)
})


#'@export
#'@name plot
#'@aliases plot,Zi,ANY-method
#'@title Plotting
#'@description plot
#'@param  x      \code{\linkS4class{Zi}}-class object
#'@param y the y coordinates of points in the plot, optional if x is an appropriate
#'structure
#'@param ... Arguments to be passed to plot
#'
#'@returns returns plot object
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'plot(Zi)
#'
setMethod("plot", c("Zi","ANY"), function(x, y, ...) {
  plot(x@deinflatedcounts, y = NULL, ...)
})
