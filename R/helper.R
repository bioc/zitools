#'@include zi_function.R
NULL

#'@name replace_phyloseq
#'@title Replace the otu table of a phyloseq object
#'@param ZiObject \code{\linkS4class{Zi}}-class object with a phyloseq object as input
#'@description Replace the OTU table of a phyloseq object with the OTU table
#'of zero de-inflated count data
#'@returns a "phyloseq"-class object
#'@importFrom phyloseq otu_table
#'@importFrom phyloseq otu_table<-
#'@importFrom phyloseq taxa_are_rows
#'@export

replace_phyloseq <- function(ZiObject)
{
  ps <- ZiObject@inputdata
  new_otu <- otu_table(ZiObject@output, taxa_are_rows(ZiObject@inputdata))
  otu_table(ps)<-new_otu
  return(ps)
}

#'@name zi2outputMatrix
#'@title Access the zero de-inflated matrix of an \code{\linkS4class{Zi}}-class object
#'
#'@param  ZiObject \code{\linkS4class{Zi}}-class object
#'@description extract the zero de-inflated matrix of an \code{\linkS4class{Zi}}-class object. The
#'output matrix is a count matrix (column = sample, row = feature) where
#'predicted structural zeros are replaced with NA
#'
#'@returns matrix
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'zi2outputMatrix(Zi)
#'
zi2outputMatrix <- function(ZiObject) {
  ZiObject@output
}

#'@name zi2countmatrix
#'@title Access the count matrix of an \code{\linkS4class{Zi}}-class object
#'
#'@param ZiObject \code{\linkS4class{Zi}}-class object
#'@description extract the count matrix (column = sample, row = feature), that
#'has been used to fit a zero inflation model
#'
#'@returns matrix
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'zi2countmatrix(Zi)

zi2countmatrix <- function(ZiObject) {
  ZiObject@countmatrix
}

#'@name tax_table
#'@title Access the taxonomy table
#'@aliases tax_table,Zi-method
#'@param object \code{\linkS4class{Zi}}-class object
#'@description access the taxonomy table of an \code{\linkS4class{Zi}}-class object if the
#'inputdata slot is a phyloseq object
#'@importFrom phyloseq tax_table
#'@returns tax_table
#'@export


setMethod("tax_table", signature = "Zi", function(object){
  tax_table <- tax_table(object@inputdata)
  return(tax_table)
})

#'@name sample_data
#'@aliases sample_data,Zi-method
#'@title Access the sample data
#'@param object \code{\linkS4class{Zi}}-class object
#'@description access the sample_data of an \code{\linkS4class{Zi}}-class object if the
#'inputdata slot is a phyloseq object
#'@importFrom phyloseq sample_data
#'@returns sample_data
#'
#'@export
#'

setMethod("sample_data", signature = "Zi", function(object){
  sample_data <- sample_data(object@inputdata)
  return(sample_data)
})

#'@name otu_table
#'@aliases otu_table,Zi-method
#'@title Access the otu table
#'@param object \code{\linkS4class{Zi}}-class object
#'@description access the otu_table of an \code{\linkS4class{Zi}}-class object if the
#'inputdata slot is a phyloseq object
#'@returns otu_table
#'@importFrom phyloseq otu_table
#'
#'@export
#'

setMethod("otu_table", signature = "Zi", function(object){
  otu_table <- otu_table(object@inputdata)
  return(otu_table)
})

#'@name phy_tree
#'@aliases phy_tree,Zi-method
#'@title Access the phylogenetic tree
#'@param physeq \code{\linkS4class{Zi}}-class object
#'@description access the phylogenetic tree (phy_tree) of an object of the class
#'"Zi" if the inputdata slot is a phyloseq object
#'@returns phy_tree
#'@importFrom phyloseq phy_tree
#'
#'@export
#'

setMethod("phy_tree", signature = "Zi", function(physeq){
  phy_tree <- phy_tree(physeq@inputdata)
  return(phy_tree)
})

#'@name rowData
#'@aliases rowData,Zi-method
#'@title Access the row data
#'@param x \code{\linkS4class{Zi}}-class object
#'@param useNames returns a rowData dataframe with rownames
#'@param ... \code{\link[SummarizedExperiment]{rowData}}
#'@description access the rowData of an \code{\linkS4class{Zi}}-class object if the inputdata
#'is an object of the class SummarizedExperiment
#'@returns DFrame
#'@importFrom SummarizedExperiment rowData
#'
#'@export
#'

setMethod("rowData", signature = "Zi", function(x, useNames = TRUE, ...){
  rowData <- rowData(x@inputdata, useNames = useNames, ...)
  return(rowData)
})

#'@name assays
#'@aliases assays,Zi-method
#'@title Access assays
#'@param x \code{\linkS4class{Zi}}-class object
#'@param withDimnames A \code{logical}, indicating whether the dimnames of the
#'SummarizedExperiment object should be applied (i.e. copied) to the extracted
#'assays. see \code{\link[SummarizedExperiment]{assays}}
#'@param ... see \code{\link[SummarizedExperiment]{assays}}
#'@description access  assays of an \code{\linkS4class{Zi}}-class object if the inputdata is an
#'object of the class SummarizedExperiment
#'@importFrom SummarizedExperiment assays
#'@returns list
#'@export


setMethod("assays", signature = "Zi", function(x, withDimnames=TRUE,  ...){
  assays <- assays(x@inputdata, withDimnames=withDimnames, ...)
  return(assays)
})


#'@name colData
#'@aliases colData,Zi-method
#'@title Access the col Data
#'@param x \code{\linkS4class{Zi}}-class object
#'@param ... \code{\link[SummarizedExperiment]{colData}}
#'@description access the colData of an \code{\linkS4class{Zi}}-class object if the inputdata is an
#'object of the class SummarizedExperiment
#'@importFrom SummarizedExperiment colData
#'@returns DFrame
#'@seealso \code{\link[SummarizedExperiment]{colData}}
#'
#'@export

setMethod("colData", signature = "Zi", function(x, ...){
  colData <- colData(x@inputdata, ...)
  return(colData)
})

#'@name t
#'@aliases t,Zi-method
#'@title Transpose a \code{\linkS4class{Zi}}-class object
#'@param x \code{\linkS4class{Zi}}-class object
#'@description transpose all matrizes of a \code{\linkS4class{Zi}}-class object
#'@returns \code{\linkS4class{Zi}}-class object
#'@export
#'

setMethod("t", signature = "Zi", definition = function(x){
  countmatrix <- t(x@countmatrix)
  output <- t(x@output)
  weights <- t(x@weights)
  result <- new(
    Class = "Zi",
    inputdata = x@inputdata,
    countmatrix = countmatrix,
    ZiModel = x@ZiModel,
    output = output,
    weights = weights
  )
  return(result)
})

#'@name zi_to_deseq2
#'@title Convert a \code{\linkS4class{Zi}}-class object to a DESeq2 dds object
#'@description A \code{\linkS4class{Zi}}-class object is converted to a DESeqDataSet object, which
#'can be used for DESeq2 analysis. Both, weight and count matrices will be
#'stored in assays of the DESeqDataSet.
#'
#'@param ZiObject \code{\linkS4class{Zi}}-class object
#'@param design  A formula which specifies the design of the experiment, taking
#' the form formula(~ x + y + z). That is, a formula with right-hand side only.
#' By default, the functions in this package and DESeq2 will use the last variable
#' in the formula (e.g. z) for presenting results (fold changes, etc.) and plotting.
#' When considering your specification of experimental design, you will want to re-order the
#' levels so that the NULL set is first.
#'@param colData if the inputdata of the \code{\linkS4class{Zi}}-class object is a matrix: a
#'DataFrame or data.frame with at least a single column. Rows of colData
#'correspond to columns of countData
#'@param ...  [phyloseq::phyloseq_to_deseq2] if the inputdata of the 'Zi'-object
#'is a phyloseq object       [DESeq2::DESeqDataSet] if the inputdata the '
#'Zi'-object is a SummarizedExperiment object
#'@importFrom phyloseq phyloseq_to_deseq2
#'@importFrom DESeq2 DESeqDataSet
#'@importFrom DESeq2 DESeqDataSetFromMatrix
#'@importFrom SummarizedExperiment assays<-
#'

zi_to_deseq2 <- function(ZiObject, design, colData, ... ){
  if (is(ZiObject@inputdata, "phyloseq") == TRUE) {
    dds <- phyloseq_to_deseq2(ZiObject@inputdata, design = design, ...)
  }
  if (is(ZiObject@inputdata, "SummarizedExperiment") == TRUE) {
    dds <- DESeqDataSet(ZiObject@inputdata, design = design, ...)
  }
  if (is(ZiObject@inputdata, "matrix") == TRUE){
    dds <- DESeqDataSetFromMatrix(ZiObject@countmatrix, colData = colData, design = design, ...)
  }
  assays(dds)[["weights"]] <- ZiObject@weights
  return(dds)
}

#'@name subset_sample
#'@title Subset a \code{\linkS4class{Zi}}-class object based on sample data
#'
#'@description Subset a \code{\linkS4class{Zi}}-class object based on sample_data of an phyloseq
#'object or on colData based on a SummarizedExperiment object
#'
#'@param Zi \code{\linkS4class{Zi}}-class object
#'@param ... The subsetting expression that should be applied, see \link[base]{subset}
#'for more details
#'
#'@export
#'@importFrom phyloseq sample_data
#'@importFrom phyloseq sample_data<-
#'@importFrom SummarizedExperiment colData
#'

subset_sample <- function(Zi, ...){
  if (is(Zi@inputdata, "phyloseq") == TRUE){
    newDF <- subset(as(sample_data(Zi@inputdata), "data.frame"), ...)
    colnames <- rownames(newDF)
    sample_data(Zi@inputdata) <- sample_data(newDF)
  }
  if (is(Zi@inputdata, "SummarizedExperiment") == TRUE){
    newDF <- subset(as(colData(Zi@inputdata), "DataFrame"), ...)
    colnames <- rownames(newDF)
    Zi@inputdata <- Zi@inputdata[,colnames]
  }
  countmatrix <- Zi@countmatrix[,colnames]
  output <- Zi@output[,colnames]
  weights <- Zi@weights[,colnames]
  result <- new(
    Class = "Zi",
    inputdata = Zi@inputdata,
    countmatrix = countmatrix,
    ZiModel = Zi@ZiModel,
    output = output,
    weights = weights
  )
  return(result)
}

#'@name subset_sample
#'@title Subset a \code{\linkS4class{Zi}}-class object based on sample data
#'
#'@description Subset a \code{\linkS4class{Zi}}-class object based on sample_data of an phyloseq
#'object or on colData of a SummarizedExperiment object
#'
#'@param Zi \code{\linkS4class{Zi}}-class object
#'@param ... The subsetting expression that should be applied, see \link[base]{subset}
#'for more details
#'
#'@export
#'@importFrom phyloseq sample_data
#'@importFrom phyloseq sample_data<-
#'@importFrom SummarizedExperiment colData
#'

subset_sample <- function(Zi, ...){
  if (is(Zi@inputdata, "phyloseq") == TRUE){
    newDF <- subset(as(tax_table(Zi@inputdata), "data.frame"), ...)
    rownames <- rownames(newDF)
    tax_table(Zi@inputdata) <- sample_data(newDF)
  }
  if (is(Zi@inputdata, "SummarizedExperiment") == TRUE){
    newDF <- subset(as(colData(Zi@inputdata), "DataFrame"), ...)
    colnames <- rownames(newDF)
    Zi@inputdata <- Zi@inputdata[,colnames]
  }
  countmatrix <- Zi@countmatrix[,colnames]
  output <- Zi@output[,colnames]
  weights <- Zi@weights[,colnames]
  result <- new(
    Class = "Zi",
    inputdata = Zi@inputdata,
    countmatrix = countmatrix,
    ZiModel = Zi@ZiModel,
    output = output,
    weights = weights
  )
  return(result)
}
#'@name subset_feature
#'@title Subset a \code{\linkS4class{Zi}}-class object based on feature data
#'
#'@description Subset a \code{\linkS4class{Zi}}-class object based on tax_table of a phyloseq
#'object or on rowData of a SummarizedExperiment object
#'
#'@param Zi \code{\linkS4class{Zi}}-class object
#'@param ... The subsetting expression that should be applied, see \link[base]{subset}
#'for more details
#'
#'@export
#'@importFrom phyloseq tax_table
#'@importFrom phyloseq tax_table<-
#'@importFrom SummarizedExperiment rowData
#'
subset_feature <- function(Zi, ...){
  if (is(Zi@inputdata, "phyloseq") == TRUE){
    mtx <- as(tax_table(Zi@inputdata), "matrix")
    newdf <- subset(data.frame(mtx), ...)
    newmtx <- as(newdf, "matrix")
    rownames <- rownames(newdf)
    tax_table(Zi@inputdata) <- tax_table(newmtx)
  }
  if (is(Zi@inputdata, "SummarizedExperiment") == TRUE){
    newDF <- subset(as(rowData(Zi@inputdata), "DataFrame"), ...)
    rownames <- rownames(newDF)
    Zi@inputdata <- Zi@inputdata[rownames,]
  }
  countmatrix <- Zi@countmatrix[rownames,]
  output <- Zi@output[rownames,]
  weights <- Zi@weights[rownames,]
  result <- new(
    Class = "Zi",
    inputdata = Zi@inputdata,
    countmatrix = countmatrix,
    ZiModel = Zi@ZiModel,
    output = output,
    weights = weights
  )
  return(result)
}
