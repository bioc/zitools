#'@name replace_phyloseq
#'@title Replace the otu table of a phyloseq object
#'@param ZiObject 'Zi'-class object with a phyloseq object as input
#'@description Replace the OTU table of a phyloseq object with an OTU table
#'where predicted structural zeros are replaced with NA
#'@returns a "phyloseq"-class object
#'@importFrom phyloseq otu_table
#'@importFrom phyloseq otu_table<-
#'@importFrom phyloseq taxa_are_rows
#'@export

replace_phyloseq <- function(ZiObject)
{
  ps <- ZiObject@datafile
  new_otu <- otu_table(ZiObject@output, taxa_are_rows(ZiObject@datafile))
  otu_table(ps)<-new_otu
  return(ps)
}

#'@name zi2outputMatrix
#'@title Access the output matrix of an 'Zi'-class object
#'
#'@param  ZiObject 'Zi'-class object
#'@description extract the output matrix of an 'Zi'-class object. The
#'output matrix is a count matrix (column = sample, row = feature) where
#'predicted structural zeros are replaced with NA
#'
#'@returns matrix
#'@export
#'
zi2outputMatrix <- function(ZiObject) {
  ZiObject@output
}

#'@name zi2inputMatrix
#'@title Access the input matrix of an 'Zi'-class object
#'
#'@param ZiObject 'Zi'-class object
#'@description extract the count matrix (column = sample, row = feature), that
#'has been used to fit a zero inflation model
#'
#'@returns matrix
#'@export
#'

zi2countmatrix <- function(ZiObject) {
  ZiObject@countmatrix
}

#'@name tax_table
#'@title Access the taxonomy table
#'
#'@param ZiObject 'Zi'-class object
#'@description access the taxonomy table of an 'Zi'-class object if the
#'datafile slot is a phyloseq object
#'@importFrom phyloseq tax_table
#'@returns tax_table
#'@export

setMethod("tax_table", signature = "Zi", function(object){
  tax_table <- tax_table(object@datafile)
  return(tax_table)
})

#'@name sample_data
#'@title Access the sample data
#'@param ZiObject 'Zi'-class object
#'@description access the sample_data of an 'Zi'-class object if the
#'datafile slot is a phyloseq object
#'@importFrom phyloseq sample_data
#'@returns sample_data
#'
#'@export
#'

setMethod("sample_data", signature = "Zi", function(object){
  sample_data <- sample_data(object@datafile)
  return(sample_data)
})

#'@name otu_table
#'@title Access the otu table
#'@param ZiObject 'Zi'-class object
#'@description access the otu_table of an 'Zi'-class object if the
#'datafile slot is a phyloseq object
#'@returns otu_table
#'@importFrom phyloseq otu_table
#'
#'@export
#'

setMethod("otu_table", signature = "Zi", function(object){
  otu_table <- otu_table(object@datafile)
  return(otu_table)
})

#'@name phy_tree
#'@title Access the phylogenetic tree
#'@param ZiObject 'Zi'-class object
#'@description access the phylogenetic tree (phy_tree) of an object of the class
#'"Zi" if the datafile slot is a phyloseq object
#'@returns phy_tree
#'@importFrom phyloseq phy_tree
#'
#'@export
#'

setMethod("phy_tree", signature = "Zi", function(physeq){
  phy_tree <- phy_tree(physeq@datafile)
  return(phy_tree)
})

#'@name rowData
#'@title Access the row data
#'@param ZiObject 'Zi'-class object
#'@description access the rowData of an 'Zi'-class object if the datafile
#'is an object of the class SummarizedExperiment
#'@returns DFrame
#'@importFrom SummarizedExperiment rowData
#'
#'@export
#'

setMethod("rowData", signature = "Zi", function(x, ...){
  rowData <- rowData(x@datafile, ...)
  return(rowData)
})

#'@name assays
#'@title Access assays
#'@param ZiObject 'Zi'-class object
#'@description access  assays of an 'Zi'-class object if the datafile is an
#'object of the class SummarizedExperiment
#'@importFrom SummarizedExperiment assays
#'@returns list
#'@export

setMethod("assays", signature = "Zi", function(x, ...){
  assays <- assays(x@datafile, ...)
  return(assays)
})


#'@name colData
#'@title Access the col Data
#'@param ZiObject 'Zi'-class object
#'@description access the colData of an 'Zi'-class object if the datafile is an
#'object of the class SummarizedExperiment
#'@importFrom SummarizedExperiment colData
#'@returns DFrame
#'
#'@export

setMethod("colData", signature = "Zi", function(x, ...){
  colData <- colData(x@datafile, ...)
  return(colData)
})

#'@name t
#'@title Transpose a 'Zi'-class object
#'@param ZiObject 'Zi'-class object
#'@description transpose all matrizes of a 'Zi'-class object
#'@returns 'Zi'-class object
#'@export
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
#'@title Convert a 'Zi'-class object to a DESeq2 dds object
#'@description A 'Zi'-class object is converted to a DESeqDataSet object, which
#'can be used for DESeq2 analysis. Both, weight and count matrices will be
#'stored in assays of the DESeqDataSet.
#'
#'@param ZiObject 'Zi'-class object
#'@param design  A formula which specifies the design of the experiment, taking
#' the form formula(~ x + y + z). That is, a formula with right-hand side only.
#' By default, the functions in this package and DESeq2 will use the last variable
#' in the formula (e.g. z) for presenting results (fold changes, etc.) and plotting.
#' When considering your specification of experimental design, you will want to re-order the
#' levels so that the NULL set is first.
#'@param colData if the datafile of the 'Zi'-class object is a matrix: a
#'DataFrame or data.frame with at least a single column. Rows of colData
#'correspond to columns of countData
#'@importFrom phyloseq phyloseq_to_deseq2
#'@importFrom DESeq2 DESeqDataSet
#'@importFrom DESeq2 DESeqDataSetFromMatrix
#'@importFrom SummarizedExperiment assays<-
#'

zi_to_deseq2 <- function(ZiObject, design, colData, ... ){
  if (is(ZiObject@datafile, "phyloseq") == TRUE) {
    dds <- phyloseq_to_deseq2(ZiObject@datafile, design = design, ...)
  }
  if (is(ZiObject@datafile, "SummarizedExperiment") == TRUE) {
    dds <- DESeqDataSet(ZiObject@datafile, design = design)
  }
  if (is(ZiObject@datafile, "matrix") == TRUE){
    dds <- DESeqDataSetFromMatrix(ZiObject@countmatrix, colData = colData, design = design, ...)
  }
  assays(dds)[["weights"]] <- ZiObject@weights
  return(dds)
}

