#'Class Zi
#'
#'Objects of this class store all the results of the ZiMain function to continue
#'zero inflated data analysis
#'@slot inputdata a matrix, phyloseq or SummarizedExperiment object.
#'@slot countmatrix matrix. The count matrix, features as rows, samples as columns
#'@slot ZiModel list. The result of fitting a zero inflated model using
#'\link[pscl]{zeroinfl}
#'@slot output matrix. The matrix where predicted structural zeros are omitted
#'and stored as NA values
#'@slot weights matrix. A matrix containing weights for zero counts
#'@exportClass Zi
#'
#'
setClass(
  Class = "Zi",
  slots = list(
    inputdata = "ANY",
    countmatrix = "matrix",
    ZiModel = "list",
    output = "matrix",
    weights = "matrix")
)

#'@name ziMain
#'@aliases ziMain,SummarizedExperiment-method
#'@aliases ziMain,phyloseq-method
#'@title  ziMain
#'@param inputdata phyloseq object, SummarizedExperiment object, or matrix (rows
#'=features, columns=samples)
#'@param feature "feature", "gene", "OTU", "phylum", etc.
#'@param formula  formula to fit the zero inflated model
#'y ~ x1 + x2 + ..., default = count ~ sample + feature.
#'A different set of regressors can be specified using y ~ x1 +x2 + ...|z1 + z2 + ...
#'where the first part describes the count data model and the second part
#'describes the zero inflation model
#'@param dist = distribution, either poisson, negative binomial ("negbin") or
#'geometric
#'@param link = link function, either "logit", "probit", "cloglog", "cauchit"
#'@param zeroRows.rm = logical, TRUE if rows that only contain zeros should be removed
#'or not (they are removed to fit a zero inflated model and will be added afterwards
#'count matrix per default = 0 and weights = 1)
#'@param ... additional parameters to describe the model, see \link[pscl]{zeroinfl}
#'
#'@description
#'This function fits a zero-inflated mixture model (either
#'Poisson or negative binomial distribution) to count data and calculates
#'weights for all zeros indicating whether a zero is a real count (weight close to 1)
#'or whether it is a structural zeros (weight close to 0). The default model is
#'a zero inflated negative binomial model.
#'
#'The input inputdata of the ziMain function is either a phyloseq
#'object, SummarizedExperiment object or count matrix.
#'
#'In order to reduce calculation times, the count matrix is divided into
#'blocks of around 5000 count values. Then, a zero inflation model (either
#'Poisson or negative binomial distribution) is fitted to the data. Using the
#'fitted zero inflated model, probabilities given that a zero in the count matrix
#'is a structural zero are predicted. Those probabilities are used in two ways:
#' 1) A zero-deinflated count matrix is generated where a appropriate proportion
#' of zeros are randomly replaced by NA. This count matrix can be used for analysis
#' methods which cannot deal with weights.
#' 2) Weights
#' #'\deqn{w = \frac{\left(1 - \pi\right) f_{\text{NB}}\left(y; \mu, \theta \right) }{f_{\text{ZINB}}\left(y;\mu, \theta, \pi\right)}.}
#' (see Van den Berge, K., Perraudeau, F., Soneson, C. et al.)
#' are calculated in order to down-weight structrual zeros in analyses
#' which can account for weighting of individual data points.
#'
#'all zero counts are calculated given the following formula:
#'
#'The result of the ziMain function can be used to analyze zero inflated count data.
#'@slot inputdata a matrix, phyloseq or SummarizedExperiment object.
#'@slot countmatrix matrix. The count matrix, features as rows, samples as columns
#'@slot ZiModel list. The result of fitting a zero inflated model using
#'\code{\link[pscl]{zeroinfl}}
#'@slot output matrix. The matrix where predicted structural zeros are omitted
#'and stored as NA values
#'@slot weights matrix. A matrix containing weights for zero counts
#'@references
#'Van den Berge, K., Perraudeau, F., Soneson, C. et al. Observation weights
#'unlock bulk RNA-seq tools for zero inflation and single-cell applications.
#'Genome Biol 19, 24 (2018). https://doi.org/10.1186/s13059-018-1406-4
#'
#'@seealso \code{\link{pscl::zeroinfl}}
#'@returns \code{\linkS4class{Zi}}-class object
#'@export
#'@examples
#'# zero-inflated count matrix
#'data(mtx)
#'# calling ziMain function:
#'Zi <- ziMain(mtx)
#'
#'

setGeneric("ziMain", function(inputdata,
                              feature = "feature",
                              formula = count ~ sample + feature,
                              dist = "negbin",
                              link = "logit",
                              zeroRows.rm = FALSE,
                              ...)
{

  mtx <- as.matrix(inputdata)

  if(is.null(rownames(mtx))) {
    rownames(mtx) <- c(1:nrow(mtx))
  }
  if(is.null(colnames(mtx))) {
    colnames(mtx) <- c(1:ncol(mtx))
  }
  if (any(mtx == 0) == FALSE) {
    result <- new(
      Class = "Zi",
      inputdata = inputdata,
      countmatrix = mtx,
      ZiModel = vector(mode = "list"),
      output = mtx,
      weights = matrix(1, nrow = nrow(mtx), ncol = ncol(mtx)))
  }
  else{
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
    inputdata = inputdata,
    countmatrix = mtx_new,
    ZiModel = ziModel,
    output = ziOutput,
    weights = weights
  )}
  return(result)
})



#'@importFrom phyloseq otu_table
setMethod(
  "ziMain",
  signature = c("phyloseq"),
  definition = function(inputdata,
                        feature = "feature",
                        formula = count ~ sample + feature,
                        dist = "negbin",
                        link = "logit",
                        zeroRows.rm = FALSE,
                        ...) {
    matrix <- as.matrix(otu_table(inputdata))
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
      inputdata = inputdata,
      countmatrix = zi_result@countmatrix,
      ZiModel = zi_result@ZiModel,
      output = zi_result@output,
      weights = zi_result@weights
    )
    return(result)
  }
)

#'@importFrom SummarizedExperiment assays
setMethod(
  "ziMain",
  signature = c("SummarizedExperiment"),
  definition = function(inputdata,
                        feature = "feature",
                        formula = count ~ sample + feature,
                        dist = "negbin",
                        link = "logit",
                        zeroRows.rm = FALSE,
                        ...) {
    matrix <- as.matrix(assays(inputdata)$counts)
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
      inputdata = inputdata,
      countmatrix = zi_result@countmatrix,
      ZiModel = zi_result@ZiModel,
      output = zi_result@output,
      weights = zi_result@weights
    )
    return(result)
  }
)



