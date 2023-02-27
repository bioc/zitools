
# setClass for Zi object
setClass(
  Class = "Zi",
  slots = list(
    design = "ANY",
    inputmatrix = "matrix",
    model = "list",
    output = "matrix",
    weights = "matrix"
  )
)
# functions to replace otu table of ps object with new otu table(str0 = NA) and get ps object
ps_replace <- function(result_zi)
{
  ps <- result_zi@design
  new_otu <- otu_table(result_zi@output, taxa_are_rows(result_zi@design))
  otu_table(ps)<-new_otu
  return(ps)
}

zi2OutputMatrix <- function(result_zi) {
  result_zi@output
}

zi2inputMatrix <- function(result_zi) {
  mtx <- result_zi@inputmatrix
  return(mtx)
}


