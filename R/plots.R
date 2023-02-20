#plots
# new boxplot method for class(zi)
boxplot.zi <- function(result_zi, log1p, ...)
{
  if (log1p == TRUE) {
    boxplot(log1p(result_zi$ziOutput), ...)
  }
  if (log1p == FALSE) {
    boxplot(result_zi$ziOutput)
  }
}

#heatmap function
heatmap <- function(input, ...) {
  UseMethod("heatmap")
}

heatmap.zi <- function(result_zi, ...) {
  df <- as.data.frame(result_zi$ziOutput)
  df <- df %>%
    filter_all(any_vars(!is.na(.)))
  mtx <- as.matrix(df)
  stats::heatmap(mtx, ...)
}
