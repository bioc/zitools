


# constructor function to create a zi object
zi <- function(ziInput, ziModel, ziOutput) {
  structure(list("ziInput" = ziInput, "ziModel" = ziModel, "ziOutput" = ziOutput), class = "zi")
}

#new method for generic boxplot()
boxplot.zi <- function(zi, plot = c("before", "after"), log1p = TRUE, ...) {
  if (log1p == TRUE ){
    if (plot == "before") {
      boxplot.default(log1p(zi$ziInput), ...)
    }
    if (plot == "after") {
      boxplot.default(log1p(zi$ziOutput), ...)
    }
  }
  if(log1p == FALSE) {
    if (plot == "before") {
      boxplot.default(zi$ziInput, ...)
    }
    if (plot == "after") {
      boxplot.default(zi$ziOutput, ...)
    }
  }
}

