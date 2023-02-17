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
