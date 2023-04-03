#'@include zi_function.R
NULL

#'@name cor
#'@title Calculate weighted Pearson Correlation coeffiecients
#'
#'@description calculate the weighted pearson correlation coefficients of a count matrix
#'of an Zi object taking weights for structural zeros into account
#'
#'@param x Zi-class object


setMethod("cor", signature = "Zi", definition = function(x, y = NULL, use = "everything", method = "pearson"){
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
