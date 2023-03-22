#'Class Zi
#'
#'Objects of this class store all the results of the ZiMain function to continue
#'zero inflated data analysis
#'@slot datafile a matrix, phyloseq or summarized experiment object.
#'@slot countmatrix matrix. The design matrix, features as rows, samples as columns
#'@slot ZiModel list. The result of fitting a zero inflated model using zeroinfl of
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
    datafile = "ANY",
    countmatrix = "matrix",
    ZINBModel = "list",
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
ps_replace <- function(ZiObject)
{
  ps <- ZiObject@datafile
  new_otu <- otu_table(ZiObject@output, taxa_are_rows(ZiObject@datafile))
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
#'@export
#'
zi2OutputMatrix <- function(ZiObject) {
  ZiObject@output
}

#'@name zi2inputMatrix
#'@title
#'
#'@param ZiObject result of the ZiMain function
#'@description extract the count matrix, that has been used to fit a zero inflation
#'model
#'
#'@returns matrix
#'@export
#'
#'
#'
zi2inputMatrix <- function(ZiObject) {
  ZiObject@countmatrix
}

#'@name tax_table
#'@title
#'
#'@param ZiObject result of the ZiMain function
#'@description acces the tax table of an Zi object if the datafile slot is a phyloseq
#'object
#'@returns tax_table
#'
#'@export
#'
#'
#'

setMethod("tax_table", signature = "Zi", function(object){
  tax_table <- phyloseq::tax_table(object@datafile)
  return(tax_table)
})

#'@name sample_data
#'@title
#'@param ZiObject result of the ZiMain function
#'@description acces the sample data of an Zi object if the datafile slot is a phyloseq
#'object
#'@returns sample_data
#'
#'@export
#'
#'
#'
#'
setMethod("sample_data", signature = "Zi", function(object){
  sample_data <- sample_data(object@datafile)
  return(sample_data)
})

#'@name otu_table
#'@title
#'@param ZiObject result of the ZiMain function
#'@description acces the otu table  of an Zi object if the datafile slot is a phyloseq
#'object
#'@returns otu_table
#'
#'@export
#'
#'
setMethod("otu_table", signature = "Zi", function(object){
  otu_table <- otu_table(object@datafile)
  return(otu_table)
})

#'@name phy_tree
#'@title
#'@param ZiObject result of the ZiMain function
#'@description acces the phy tree  of an Zi object if the datafile slot is a phyloseq
#'object
#'@returns phy_tree
#'
#'@export
#'
#'
setMethod("phy_tree", signature = "Zi", function(physeq){
  phy_tree <- phy_tree(physeq@datafile)
  return(phy_tree)
})

#'@name rowData
#'@title
#'@param ZiObject result of the ZiMain function
#'@description access the rowData of an Zi object if the datafile is an SummarizedExperiment
#'object
#'@returns DFrame
#'
#'@export
#'
#'
#'
#'
setMethod("rowData", signature = "Zi", function(x, ...){
  rowData <- rowData(x@datafile, ...)
  return(rowData)
})

#'@name assays
#'@title
#'@param ZiObject result of the ZiMain function
#'@description access the assays of an Zi object if the datafile is an SummarizedExperiment
#'object
#'@returns list
#'@export
#'
setMethod("assays", signature = "Zi", function(x, ...){
  assays <- assays(x@datafile, ...)
  return(assays)
})

#'@name metadata
#'@title
#'@param ZiObject result of the ZiMain function
#'@description access the metadata of an Zi object if the datafile is an SummarizedExperiment
#'object
#'@returns list
#'@export
#'
setMethod("metadata", signature = "Zi", function(x, ...){
  metadata <- metadata(x@datafile, ...)
  return(metadata)
})

#'@name colData
#'@title
#'@param ZiObject result of the ZiMain function
#'@description access the colData of an Zi object if the datafile is an SummarizedExperiment
#'object
#'@returns DFrame
#'
#'@export
#'
setMethod("colData", signature = "Zi", function(x, ...){
  colData <- colData(x@datafile, ...)
  return(colData)
})

#'@name t
#'@title
#'@param ZiObject result of the ZiMain function
#'@description transpose of matrix of a Zi object
#'@returns Zi object
#'@export
#'
#'
#'

setMethod("t", signature = "Zi", definition = function(x){
  countmatrix <- t(x@countmatrix)
  output <- t(x@output)
  weights <- t(x@weights)
  result <- new(
    Class = "Zi",
    datafile = x@datafile,
    countmatrix = countmatrix,
    ZiModel = x@ZiModel,
    output = output,
    weights = weights
  )
  return(result)
})

#'@name zi_to_deseq2
#'@title Convert zi data to DESeq2 dds object
#'@description Zi data is converted to a DESeqDataSet object, which can be used for
#'DESeq2 Analysis. Both, weight and count matrices will be stored in assays of
#'the DESeqDataSet
#'
#'@param zi Zi-class object
#'@param design  A formula which specifies the design of the experiment, taking
#' the form formula(~ x + y + z). That is, a formula with right-hand side only.
#' By default, the functions in this package and DESeq2 will use the last variable
#' in the formula (e.g. z) for presenting results (fold changes, etc.) and plotting.
#' When considering your specification of experimental design, you will want to re-order the
#' levels so that the NULL set is first.
#'@param colData for matrix input: a DataFrame or data.frame with at least a single
#'column. Rows of colData correspond to columns of countData
#'
#'
#'
?phyloseq_to_deseq2

zi_to_deseq2 <- function(zi, design, colData, ... ){
  if (is(zi@datafile, "phyloseq") == TRUE) {
    dds <- phyloseq_to_deseq2(zi@datafile, design = design, ...)
  }
  if (is(zi@datafile, "SummarizedExperiment") == TRUE) {
    dds <- DESeqDataSet(zi@datafile, design = design)
  }
  if (is(zi@datafile, "matrix") == TRUE){
    dds <- DESeqDataSetFromMatrix(zi@countmatrix, colData = colData, design = design, ...)
  }
  assays(dds)[["weights"]] <- zi@weights
  return(dds)
}

