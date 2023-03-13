#'Class Zi
#'
#'Objects of this class store all the results of the ZiMain function to continue
#'zero inflated data analysis
#'@slot design a matrix, phyloseq or summarized experiment object.
#'@slot inputmatrix matrix. The design matrix, features as rows, samples as columns
#'@slot model list. The result of fitting a zero inflated model using zeroinfl of
#'the pscl package
#'@slot output matrix. The matrix where predicted structural zeros are ommitted
#'and stored as NA values
#'@slot weights matrix. A matrix containing weights for zero counts
#'
#'@exportClass Zi
#'
#'
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

#'@name ps_replace
#'
#'@title
#'
#'@param ZiObject Zi Object with a phyloseq object as input
#'
#'@description Replace the OTU table of a phyloseq object with an OTU table where
#'predicted structural zeros are replaced with NA
#'
#'@returns a phyloseq object
#'
#'@example
#'
#'
#'
ps_replace <- function(result_zi)
{
  ps <- result_zi@design
  new_otu <- otu_table(result_zi@output, taxa_are_rows(result_zi@design))
  otu_table(ps)<-new_otu
  return(ps)
}

#'@name zi2OutputMatrix
#'@title
#'
#'@param  ZiObject ZiObject, result of the ZiMain function
#'@description extract the NA, structural zeros replaced with NA matrix from a
#'given ZiObject
#'
#'@returns matrix
#'
#'
zi2OutputMatrix <- function(result_zi) {
  result_zi@output
}

#'@name zi2inputMatrix
#'@title
#'
#'@param ZiObject result of the ZiMain function
#'@description extract the count matrix, that has been used to fit a zero inflation
#'model
#'
#'@returns matrix
#'
#'
#'
#'


zi2inputMatrix <- function(result_zi) {
  mtx <- result_zi@inputmatrix
  return(mtx)
}


