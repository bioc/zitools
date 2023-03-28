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
    ZiModel = "list",
    output = "matrix",
    weights = "matrix")
)

#'@name ziMain
#'
#'@title  ziMain
#'
#'
#'
#'@param datafile matrix (rows = features, columns = samples), phyloseq object or
#'SummarizedExperiment object
#'@param feature "gene", "OTU", "phylum", etc.
#'@param formula  formula to fit the model response ~ predictor1 + predictor2 + ..., default = count ~ sample+OTU
#'@param dist = distribution, either poisson, negative binomial = negbin or geometric
#'@param link = link function, either logit, probit, cloglog, cauchit
#'@param zeroRows.rm = logical, TRUE if rows that only contain zeros should be removed
#'or not (they are removed to fit a zero inflated model and will be added afterwards
#'count matrix per default = 0 and weights = 1)
#'
#'
#'@description The ziMain function uses a matrix, phyloseq, or SummarizedExperiment
#'object, extracts the count matrix to fit a zero inflation model to the data.
#'The matrix is divided into blocks of around 5000 count values to improve run time.
#'Predicted probabilities given that a zero in the count matrix  is a
#'structural zero are used to draw structural zeros and replace them with NA. Further,
#'weights for all zero counts are calculated given the following formula:
#'w=...
#'
#'
#'@returns S4 list with slots for the input object, the extracted count matrix if
#'the input is not a matrix, the results of the fitted zero inflation model, a
#'matrix where the predicted structural zeros are replaced with NA, a matrix
#'containing the calculated weights
#'
#'@export
#'@example
#'
#'

setGeneric("ziMain", function(datafile,
                              feature = "",
                              formula,
                              dist = c("poisson", "negbin", "geometric"),
                              link = c("logit", "probit", "cloglog", "cauchit",                                  "log"),
                              zeroRows.rm = FALSE,
                              ...)
{
  mtx <- as.matrix(datafile)
  if(is.null(rownames(mtx))) {
    rownames(mtx) <- c(1:nrow(mtx))
  }
  if(is.null(colnames(mtx))) {
    colnames(mtx) <- c(1:ncol(mtx))
  }
  rownames <- rownames(mtx)
  colnames <- colnames(mtx)
  mtx_new <- mtx[rowSums(mtx[]) > 0,] #remove rows that contain only 0
  mtx_random <- preprocess_mtx(mtx)
  list_subset <- subset_mtx(mtx_random)
  list_core <- list()
  for (i in 1:length(list_subset)) {
    list_core[[i]] <-
      zi_core(
        list_subset[[i]],
        feature = feature,
        formula = formula,
        dist = dist,
        link = link,
        ...
      )
  }
  ziInput <- do.call(rbind, lapply(list_core, '[[', "ziInput"))
  ziModel <- lapply(list_core, '[[', "ziModel")
  ziOutput <- do.call(rbind, lapply(list_core, '[[', "ziOutput"))
  weights <- do.call(rbind, lapply(list_core, '[[', "weights"))
  if (zeroRows.rm == FALSE) {
    mtx_new <- mtx
    ziOutput <- rbind(ziOutput, mtx[rowSums(mtx[]) == 0,])
    zero_weights <- mtx[rowSums(mtx[]) == 0,]
    zero_weights[] <- 1
    weights <- rbind(weights, zero_weights)
  }
  if(zeroRows.rm == TRUE) {
    rownames <- rownames(mtx_new)
    colnames <- colnames(mtx_new)
  }
  mtx_new <- mtx_new[rownames,colnames]
  mode(mtx_new) <- "integer"
  ziOutput <- ziOutput[rownames,colnames]
  mode(ziOutput) <- "integer"
  weights <- weights[rownames,colnames]
  result <- new(
    Class = "Zi",
    datafile = datafile,
    countmatrix = mtx_new,
    ZiModel = ziModel,
    output = ziOutput,
    weights = weights
  )
  return(result)
})

setMethod(
  "ziMain",
  signature = c("phyloseq"),
  definition = function(datafile,
                        feature = "",
                        formula,
                        dist = c("poisson", "negbin", "geometric"),
                        link = c("logit", "probit", "cloglog", "cauchit", "log"),
                        zeroRows.rm = FALSE,
                        ...) {
    matrix <- as.matrix(otu_table(datafile))
    suppressWarnings(class(matrix) <- "matrix")
    zi_result <-
      ziMain(
        matrix,
        feature = feature,
        formula = formula,
        dist = dist,
        link = link,
        zeroRows.rm = zeroRows.rm,
        ...
      )
    result <- new(
      Class = "Zi",
      datafile = datafile,
      countmatrix = zi_result@countmatrix,
      ZiModel = zi_result@ZiModel,
      output = zi_result@output,
      weights = zi_result@weights
    )
    return(result)
  }
)

setMethod(
  "ziMain",
  signature = c("SummarizedExperiment"),
  definition = function(datafile,
                        feature = "",
                        formula,
                        dist = c("poisson", "negbin", "geometric"),
                        link = c("logit", "probit", "cloglog", "cauchit", "log"),
                        zeroRows.rm = FALSE,
                        ...) {
    matrix <- as.matrix(assays(datafile)$counts)
    zi_result <-
      ziMain(
        matrix,
        feature = feature,
        formula = formula,
        dist = dist,
        link = link,
        zeroRows.rm = zeroRows.rm,
        ...
      )
    result <- new(
      Class = "Zi",
      datafile = datafile,
      countmatrix = zi_result@countmatrix,
      ZiModel = zi_result@ZiModel,
      output = zi_result@output,
      weights = zi_result@weights
    )
    return(result)
  }
)



