#plots
# new boxplot method for class(zi)
boxplot.Zi <- function(result_zi, log1p =FALSE, ...)
{
  if (log1p == TRUE) {
    boxplot(log1p(result_zi@output), ...)
  }
  if (log1p == FALSE) {
    boxplot(result_zi@output)
  }
}

#heatmap function
heatmap <- function(input, ...) {
  UseMethod("heatmap")
}

heatmap.Zi <- function(result_zi, ...) {
  df <- as.data.frame(result_zi@output)
  df <- df %>%
    filter_all(any_vars(!is.na(.)))
  mtx <- as.matrix(df)
  stats::heatmap(mtx, ...)
}
