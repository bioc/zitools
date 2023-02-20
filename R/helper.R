
# constructor function to create a zi object
zi <- function(ziInput, ziModel, ziOutput, weights) {
  structure(list("ziInput" = ziInput, "ziModel" = ziModel, "ziOutput" = ziOutput, "weights" = weights), class = c("zi","matrix"))
}

# functions to replace otu table of ps object with new otu table(str0 = NA) and get ps object
ps_replace <- function(result_zi)
{
  ps <- result_zi$ziInput
  new_otu <- otu_table(result_zi$ziOutput, taxa_are_rows(result_zi$ziInput))
  otu_table(ps)<-new_otu
  return(ps)
}

# get Output Matrix
zi2OutputMatrix <- function(result_zi) {
  result_zi$ziOutput
}

#get input matrix
zi2inputMatrix <- function(result_zi) {
  if ("matrix" %in% class(result_zi$ziInput)){
    mtx <- result_zi$ziInput
  }
  if ("phyloseq" %in% class(result_zi$ziInput)) {
    mtx <- otu_table(result_zi$ziInput) %>%
      as.matrix()
  }
  if ("SummarizedExperiment" %in% class(result_zi$ziInput) | "RangedSummarizedExperiment" %in% class(result_zi$ziInput)) {
    mtx <- assays(result_zi$ziInput)$counts
  }
  class(mtx) <- c("matrix", "array")
  return(mtx)
}


