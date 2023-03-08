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
