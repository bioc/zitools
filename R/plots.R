#'@name boxplot
#'@title Create boxplots of a 'Zi'-class object
#'@description Create boxplots of a 'Zi'-class object.
#'@param x 'Zi'-class object
#'@param log1p logical, default = FALSE, if TRUE log(1+p) transformation takes
#'place
#'@param ... see \link[graphics]{boxplot.default}
#'@seealso \link[graphics]{boxplot.default}
#'@importFrom graphics boxplot
#'@export
#'@examples
#'boxplot(Zi, log1p = TRUE)
#'
#'

boxplot.Zi <- function(x, log1p =FALSE, ...)
{
  if (log1p == TRUE) {
    boxplot(log1p(x@output), ...)
  }
  if (log1p == FALSE) {
    boxplot(x@output, ...)
  }
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
#'@description
#'This function illustrates replacement of structural zeros with NA (zero-deinflation) as heamap.
#'NAs are highlighted in black color.
#'The function also enables reordering of the rows according to proportions of NA.
#'
#'
#'@param ZiObject ZiObject, result of the ziMain function
#'@param reorderRows If TRUE, rows with least NAs are at the top.
#'
#'@returns heatmap
#'
#'
MissingValueHeatmap <- function(ZiObject,reorderRows=FALSE) {
  mtx <- log1p(ZiObject@output)
  # Sort matrix according to number of missing values
  if(reorderRows)
    mtx <- mtx[order(rowSums(mtx,na.rm=T)), ]
  stats::heatmap(mtx,Rowv=NA,Colv=NA,labRow=FALSE,
                 na.rm=T,col=c(RColorBrewer::brewer.pal(n = 9, name = "Blues")),
                 scale = "none", margins = c(11,0), cexCol=1, na.value="red")
}




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

