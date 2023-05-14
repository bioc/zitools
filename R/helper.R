#'@include ziMain.R
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
  new_otu <- otu_table(ZiObject@deinflatedcounts, taxa_are_rows(ZiObject@inputdata))
  otu_table(ps)<-new_otu
  return(ps)
}

#'@name inputdata
#'@title Access the inputdata
#'@param x \code{\linkS4class{Zi}}-class object
#'@param value replacement of inputdata
#'@description access the inputdata of an \code{\linkS4class{Zi}}-class object
#'@returns inputdata
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'inputdata(Zi)

setGeneric("inputdata", function(x) standardGeneric("inputdata"))
setGeneric("inputdata<-", function(x, value) standardGeneric("inputdata<-"))

#'@name inputdata
#'@aliases inputdata,Zi-method
#'@export
#'@rdname inputdata
#'
#'
setMethod("inputdata", "Zi", function(x) x@inputdata)
setMethod("inputdata<-", "Zi", function(x, value) {
  x@inputdata <- value
  x
})
#'@name inputcounts
#'@title Access the inputcounts
#'@param x \code{\linkS4class{Zi}}-class object
#'@param value replacement of inputcounts
#'@description access the inputcounts of an \code{\linkS4class{Zi}}-class object
#'@returns inputcounts
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'inputcounts(Zi)

setGeneric("inputcounts", function(x) standardGeneric("inputcounts"))
setGeneric("inputcounts<-", function(x, value) standardGeneric("inputcounts<-"))

#'@name inputcounts
#'@aliases inputcounts,Zi-method
#'@export
#'@rdname inputcounts

setMethod("inputcounts", "Zi", function(x) x@inputcounts)
setMethod("inputcounts<-", "Zi", function(x, value) {
  x@inputcounts <- value
  x
})

#'@name model
#'@title Access the model
#'@param x \code{\linkS4class{Zi}}-class object
#'@param value replacement of model
#'@description access the model of an \code{\linkS4class{Zi}}-class object
#'@returns model
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'model(Zi)

setGeneric("model", function(x) standardGeneric("model"))
setGeneric("model<-", function(x, value) standardGeneric("model<-"))

#'@name model
#'@aliases model,Zi-method
#'@export
#'@rdname model

setMethod("model", "Zi", function(x) x@model)
setMethod("model<-", "Zi", function(x, value) {
  x@model <- value
  x
})

#'@name deinflatedcounts
#'@title Access the model
#'@param x \code{\linkS4class{Zi}}-class object
#'@param value replacement of deinflatedcounts
#'@description access the deinflatedcounts of an \code{\linkS4class{Zi}}-class object
#'@returns deinflatedcounts
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'deinflatedcounts(Zi)

setGeneric("deinflatedcounts", function(x) standardGeneric("deinflatedcounts"))
setGeneric("deinflatedcounts<-", function(x, value) standardGeneric("deinflatedcounts<-"))

#'@name deinflatedcounts
#'@aliases deinflatedcounts,Zi-method
#'@export
#'@rdname deinflatedcounts

setMethod("deinflatedcounts", "Zi", function(x) x@deinflatedcounts)
setMethod("deinflatedcounts<-", "Zi", function(x, value) {
  x@deinflatedcounts <- value
  x
})

#'@name weights
#'@title Access the weights
#'@param x \code{\linkS4class{Zi}}-class object
#'@param value replacement of deinflatedcounts
#'@description access the weights of an \code{\linkS4class{Zi}}-class object
#'@returns weights
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'weights(Zi)

setGeneric("weights", function(x) standardGeneric("weights"))
setGeneric("weights<-", function(x, value) standardGeneric("weights<-"))

#'@name weights
#'@aliases weights,Zi-method
#'@export
#'@rdname weights

setMethod("weights", "Zi", function(x) x@weights)
setMethod("weights<-", "Zi", function(x, value) {
  x@weights <- value
  x
})

#'@name tax_table
#'@title Access the taxonomy table
#'@aliases tax_table,Zi-method
#'@param object \code{\linkS4class{Zi}}-class object
#'@description access the taxonomy table (\link[phyloseq]{tax_table}) of an \code{\linkS4class{Zi}}-class object if the
#'inputdata slot is a phyloseq object
#'@importFrom phyloseq tax_table
#'@returns tax_table
#'@export


setMethod("tax_table", signature = "Zi", function(object){
  if("phyloseq" %in% class(object@inputdata)) {
  tax_table <- tax_table(object@inputdata)}
  if("matrix" %in% class(object@inputdata)){
    tax_table <- NULL}
  if("SummarizedExperiment" %in% class(object@inputdata)){
    tax_table <- NULL}
  return(tax_table)
})

#'@name sample_data
#'@aliases sample_data,Zi-method
#'@title Access the sample data
#'@param object \code{\linkS4class{Zi}}-class object
#'@description access the \link[phyloseq]{sample_data} of an \code{\linkS4class{Zi}}-class object if the
#'inputdata slot is a phyloseq object
#'@importFrom phyloseq sample_data
#'@returns sample_data
#'
#'@export
#'

setMethod("sample_data", signature = "Zi", function(object){
  if("phyloseq" %in% class(object@inputdata)) {
  sample_data <- sample_data(object@inputdata)}
  if("matrix" %in% class(object@inputdata)){
    sample_data <- NULL}
  if("SummarizedExperiment" %in% class(object@inputdata)){
    sample_data <- NULL}
  return(sample_data)
})

#'@name otu_table
#'@aliases otu_table,Zi-method
#'@title Access the otu table
#'@param object \code{\linkS4class{Zi}}-class object
#'@description access the \link[phyloseq]{otu_table} of an \code{\linkS4class{Zi}}-class object if the
#'inputdata slot is a phyloseq object
#'@returns otu_table
#'@importFrom phyloseq otu_table
#'
#'@export
#'

setMethod("otu_table", signature = "Zi", function(object){
  if("phyloseq" %in% class(object@inputdata)) {
  otu_table <- otu_table(object@inputdata)}
  if("matrix" %in% class(object@inputdata)){
    otu_table <- NULL}
  if("SummarizedExperiment" %in% class(object@inputdata)){
    otu_table <- NULL}
  return(otu_table)
})

#'@name phy_tree
#'@aliases phy_tree,Zi-method
#'@title Access the phylogenetic tree
#'@param physeq \code{\linkS4class{Zi}}-class object
#'@description access the phylogenetic tree (\link[phyloseq]{phy_tree}) of an object of the class
#'"Zi" if the inputdata slot is a phyloseq object
#'@returns phy_tree
#'@importFrom phyloseq phy_tree
#'
#'@export
#'

setMethod("phy_tree", signature = "Zi", function(physeq){
  if("phyloseq" %in% class(physeq@inputdata)) {
  phy_tree <- phy_tree(physeq@inputdata)}
  if("matrix" %in% class(physeq@inputdata)){
    phy_tree <- NULL}
  if("SummarizedExperiment" %in% class(physeq@inputdata)){
    phy_tree <- NULL}
  return(phy_tree)
})

#'@name rowData
#'@aliases rowData,Zi-method
#'@title Access the row data
#'@param x \code{\linkS4class{Zi}}-class object
#'@param useNames returns a rowData dataframe with rownames
#'@param ... \code{\link[SummarizedExperiment]{rowData}}
#'@description access the \link[SummarizedExperiment]{rowData} of an \code{\linkS4class{Zi}}-class object if the inputdata
#'is an object of the class SummarizedExperiment
#'@returns DFrame
#'@importFrom SummarizedExperiment rowData
#'
#'@export
#'

setMethod("rowData", signature = "Zi", function(x, useNames = TRUE, ...){
  if("SummarizedExperiment" %in% class(x@inputdata)) {
  rowData <- rowData(x@inputdata, useNames = useNames, ...)}
  if("matrix" %in% class(x@inputdata)){
  rowData <- NULL}
  if("phyloseq" %in% class(x@inputdata)){
  rowData <- NULL}
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
#'@description access  \link[SummarizedExperiment]{assays} of an \code{\linkS4class{Zi}}-class object if the inputdata is an
#'object of the class SummarizedExperiment
#'@importFrom SummarizedExperiment assays
#'@returns list
#'@export


setMethod("assays", signature = "Zi", function(x, withDimnames=TRUE,  ...){
  if("SummarizedExperiment" %in% class(x@inputdata)) {
  assays <- assays(x@inputdata, withDimnames=withDimnames, ...)}
  if("matrix" %in% class(x@inputdata)){
    assays <- NULL}
  if("phyloseq" %in% class(x@inputdata)){
    assays <- NULL}
  return(assays)
})


#'@name colData
#'@aliases colData,Zi-method
#'@title Access the col Data
#'@param x \code{\linkS4class{Zi}}-class object
#'@param ... \code{\link[SummarizedExperiment]{colData}}
#'@description access the \link[SummarizedExperiment]{colData} of an \code{\linkS4class{Zi}}-class object if the inputdata is an
#'object of the class SummarizedExperiment
#'@importFrom SummarizedExperiment colData
#'@returns DFrame
#'@seealso \code{\link[SummarizedExperiment]{colData}}
#'
#'@export

setMethod("colData", signature = "Zi", function(x, ...){
  if("SummarizedExperiment" %in% class(x@inputdata)) {
  colData <- colData(x@inputdata, ...)}
  if("matrix" %in% class(x@inputdata)){
    colData <- NULL}
  if("phyloseq" %in% class(x@inputdata)){
    colData <- NULL}
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
  inputcounts <- t(x@inputcounts)
  deinflatedcounts <- t(x@deinflatedcounts)
  weights <- t(x@weights)
  result <- new(
    Class = "Zi",
    inputdata = x@inputdata,
    inputcounts = inputcounts,
    model = x@model,
    deinflatedcounts = deinflatedcounts,
    weights = weights
  )
  return(result)
})

#'@name zi2deseq2
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
#'@export
#'@examples
#'data(mtx)
#'Zi <- ziMain(mtx)
#'colData <- data.frame(group = factor(x = c(1,1,1,1,1,2,2,2,2,2)))
#'zi2deseq2(Zi, ~group, colData)

zi2deseq2 <- function(ZiObject, design, colData, ... ){
  if (is(ZiObject@inputdata, "phyloseq") == TRUE) {
    dds <- phyloseq_to_deseq2(ZiObject@inputdata, design = design, ...)
  }
  if (is(ZiObject@inputdata, "SummarizedExperiment") == TRUE) {
    dds <- DESeqDataSet(ZiObject@inputdata, design = design, ...)
  }
  if (is(ZiObject@inputdata, "matrix") == TRUE){
    dds <- DESeqDataSetFromMatrix(ZiObject@inputcounts, colData = colData, design = design, ...)
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
  inputcounts <- Zi@inputcounts[,colnames]
  deinflatedcounts <- Zi@deinflatedcounts[,colnames]
  weights <- Zi@weights[,colnames]
  result <- new(
    Class = "Zi",
    inputdata = Zi@inputdata,
    inputcounts = inputcounts,
    model = Zi@model,
    deinflatedcounts = deinflatedcounts,
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
  inputcounts <- Zi@inputcounts[,colnames]
  deinflatedcounts <- Zi@deinflatedcounts[,colnames]
  weights <- Zi@weights[,colnames]
  result <- new(
    Class = "Zi",
    inputdata = Zi@inputdata,
    inputcounts = inputcounts,
    model = Zi@model,
    deinflatedcounts = deinflatedcounts,
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
  inputcounts <- Zi@inputcounts[rownames,]
  deinflatedcounts <- Zi@deinflatedcounts[rownames,]
  weights <- Zi@weights[rownames,]
  result <- new(
    Class = "Zi",
    inputdata = Zi@inputdata,
    inputcounts = inputcounts,
    model = Zi@model,
    deinflatedcounts = deinflatedcounts,
    weights = weights
  )
  return(result)
}


#'@name resample_deinflatedcounts
#'@title Resample a \code{\linkS4class{Zi}}-class object
#'
#'@description Resample the deinflatedcounts matrix of an \code{\linkS4class{Zi}}-class object.
#'Resampling is done by drawing from a binomial distribution with a given probability
#'that a count value (zero and non-zero) is a structural zero.
#'
#'
#'@param x \code{\linkS4class{Zi}}-class object
#'
#'@export

resample_deinflatedcounts <- function(x) {
  mtx <- x@inputcounts
  mtx_new <- mtx[rowSums(mtx[]) > 0,]
  feature <- colnames(x@model[[1]][["model"]])[3]
  rownames <- rownames(mtx)
  colnames <- colnames(mtx)
  list_deinflatedcounts <- list()
  for(i in 1:length(x@model)){
    vec <- x@model[[i]][["model"]][, 3]
    vec <- vec[1:(length(vec) / ncol(mtx_new))]
    count_sub <- mtx_new[vec, ]
    count_long <- reshape_zi(count_sub, feature = feature)
    new_deinflatedcounts <-
      omit_str_zero(x@model[[i]], count_long, feature = feature)
    new_deinflatedcounts <- new_deinflatedcounts %>%
      spread(key = "sample", value = "count") %>%
      column_to_rownames(var = feature) %>%
      as.matrix()
    list_deinflatedcounts[[i]] <- new_deinflatedcounts
  }
  deinflatedcounts <- do.call(rbind, list_deinflatedcounts)
  deinflatedcounts <- rbind(deinflatedcounts, mtx[rowSums(mtx[])==0,])
  deinflatedcounts <- deinflatedcounts[rownames,colnames]
  result <- new(Class = "Zi",
                inputdata = x@inputdata,
                inputcounts = x@inputcounts,
                model = x@model,
                deinflatedcounts = deinflatedcounts,
                weights = x@weights)
  return(result)
}
