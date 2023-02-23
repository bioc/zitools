library("usethis")
library("devtools")
library("roxygen2")


#given matrix - reshape it into long format as input for zeroinflation model, df_wide = given matrix, feature = rows, sample = columns, feature = "OTU", "gene", "species", etc
reshape_zi <- function(mtx, feature = "") {
  zi_long <- data.frame(mtx) %>%
    rownames_to_column(var = feature)%>%
    gather(key = "sample", value = "count", 1:ncol(mtx)+1)
  return(zi_long)
}

#draw structural zeros with a given probability (p_str_zero) compare if
#predicted zero = zero count –> set to NA, repeat until NA >= sum(p_str_zero) –> TRUE replace excess of NA with zero (random sampling)
#zi_prediction_wide = new dataframe as matrix, zi_prediction_long = new dataframe long format
#zi = zero inflation model (output of pscl package)
#zi_input = long format of the input matrix
#zi_prediction = dataframe with calculated count and p_str_zero
omit_str_zero <- function(zi, zi_input, feature = "") {
  zi_prediction <- data.frame(predicted_count = predict(zi, type = "count"))%>%
    bind_cols(data.frame(p_str_zero = predict(zi, type = "zero")))%>%
    bind_cols(zi_input, .)
  str_zero <- bind_cols(zi_prediction, predicted_zero= c(rbinom(nrow(zi_prediction), 1, prob = zi_prediction$p_str_zero)))
  repeat {
    str_zero <- str_zero %>%
      mutate(predicted_zero = replace(predicted_zero, 1:nrow(str_zero), c(rbinom(nrow(str_zero),1, prob= str_zero$p_str_zero))))
    str_zero$count <- ifelse(str_zero$count == 0 & str_zero$predicted_zero == 1, NA, str_zero$count)
    str_zero$na <- ifelse(is.na(str_zero$count),1,0)
    if(sum(str_zero$na) >= round(sum(str_zero$p_str_zero), digit = 0)) {
      break
    }
  }
  replaced_zero <- str_zero %>%
    filter(is.na(count)) %>%
    sample_n(sum(str_zero$na)-round(sum(str_zero$p_str_zero)), replace = FALSE) %>%
    replace_na(list(count = 0))
  zi_replaced <- merge(str_zero, replaced_zero, by.x=1:2, by.y=1:2, all = TRUE)
  zi_replaced <- cbind(zi_replaced[c(1,2)], "count" = with(zi_replaced, ifelse(is.na(count.y), count.x, count.y)))

  return(zi_replaced)
}

#function to calculate weights - formula: w=(1-pi)*fnb/fzinb
calcWeights <- function(zi_input, zi, feature) {
  df<-data.frame(zi_input, pstr0 = predict(zi, type = "zero"))
  df <- mutate(df, pstr0 = case_when( count != 0~0, TRUE ~ pstr0))
  f_nb<- dnbinom(x = zi_input$count, size =zi[["theta"]], mu = predict(zi, type = "count"), log = FALSE )
  f_zinb <- dzinegbin(x = zi_input$count, size = zi[["theta"]], munb = predict(zi, type = "count"), pstr0 = df$pstr0, log = FALSE)
  df$f_nb <- f_nb
  df$f_zinb <- f_zinb
  df$weights <- c(((1-df$pstr0)*df$f_nb)/df$f_zinb)
  df$"1-pi" <- c(1-df$pstr0)
  df_wide <- df %>%
    select(c(feature, "sample", "weights"))%>%
    spread(key="sample", value = "weights")%>%
    column_to_rownames(var=feature)%>%
    as.matrix(.)
  return(df_wide)
}

#subset matrix function - subset the matrix into groups of roughly 5000 datapoints per block
subset_mtx <- function(mtx)
{
  a <- round(nrow(mtx)*ncol(mtx)/5000) # number of blocks
  if (a == 0) {a <- 1}
  b <- ceiling(nrow(mtx)/a) #number of rows for each block
  index <- rep(c(1:a),each=b) #index for subsetted groups
  index <- index[1:nrow(mtx)]
  df <- as.data.frame(mtx)
  df$index <- index
  list_subset <- split(df, index) # split df based on index
  list_subset <- lapply(list_subset, function(x) x[!(names(x) %in% c("index"))]) # remove index
  return(list_subset)}

#zi core function
zi_core <- function(input, feature = "",  formula, dist = c("poisson", "negbin", "geometric"), link = c("logit", "probit", "cloglog", "cauchit", "log"), ...)
{
  matrix <- as.matrix(input)
  zi_input <- reshape_zi(input, feature)
  zi <- zeroinfl(formula, data = zi_input, dist = dist,
                 link = link, ...)
  zi_prediction_long <- omit_str_zero(zi, zi_input, feature)
  zi_prediction_wide <- zi_prediction_long %>%
    spread(key = "sample", value = "count") %>%
    column_to_rownames(var = feature)%>%
    as.matrix()
  weights <- calcWeights(zi_input, zi, feature)
  result <- list(ziInput = matrix, ziModel = zi, ziOutput = zi_prediction_wide, weights = weights)
  return(result)
}

preprocess_mtx <- function(mtx){
  mtx <- mtx[rowSums(mtx[]) >0, ]
  seed <- .Random.seed
  set.seed(seed)
  random_rows<-sample(nrow(mtx))
  mtx <- mtx[random_rows,]
  .Random.seed <- seed
  return(mtx)
}

#zi_main function - generic function: default method for class(object)=matrix, further methods for phyloseq and summarized experiment objects defined
zi_main <- function(object, ...) UseMethod("zi_main")

zi_main.default <- function(input, feature = "",  formula, dist = c("poisson", "negbin", "geometric"), link = c("logit", "probit", "cloglog", "cauchit", "log"), zeroRows.rm = FALSE,  ...)
{
  mtx <- as.matrix(input)
  mtx_new <- preprocess_mtx(mtx)
  list_subset <- subset_mtx(mtx_new)
  list_core <- list()
  for(i in 1:length(list_subset)){
    list_core[[i]] <- zi_core(list_subset[[i]], feature = feature, formula = formula, dist = dist, link = link, ...)
  }
  ziInput <- do.call(rbind, lapply(list_core, '[[', "ziInput"))
  ziModel <- lapply(list_core, '[[', "ziModel")
  ziOutput <- do.call(rbind, lapply(list_core, '[[', "ziOutput"))
  weights <- do.call(rbind, lapply(list_core, '[[', "weights"))
  if (zeroRows.rm == FALSE) {
    ziOutput <- rbind(ziOutput, mtx[rowSums(mtx[]) == 0, ])
    zero_weights <- mtx[rowSums(mtx[]) == 0, ]
    zero_weights[]<- 1
    weights <- rbind(weights, zero_weights) }
  ziOutput <- ziOutput[order(match(rownames(ziOutput), row.names(mtx))), , drop = FALSE]
  weights <- weights[order(match(rownames(weights), rownames(mtx))), , drop = FALSE]
  result <- zi(input, ziModel, ziOutput, weights)
  return(result)
}

zi_main.phyloseq <- function(input, feature = "",  formula, dist = c("poisson", "negbin", "geometric"), link = c("logit", "probit", "cloglog", "cauchit", "log"), ...)
{
  matrix <- as.matrix(otu_table(input))
  class(matrix) <- "matrix"
  zi_result <- zi_main(matrix, feature = feature, formula = formula, dist = dist, link = link, ...)
  result <- zi(input, zi_result$ziModel, zi_result$ziOutput, zi_result$weights)
  return(result)
}

zi_main.SummarizedExperiment <- function(input,  feature = "",  formula, dist = c("poisson", "negbin", "geometric"), link = c("logit", "probit", "cloglog", "cauchit", "log"), ...)
{
  matrix <- as.matrix(assays(input)$counts)
  zi_result <- zi_main(matrix, feature = feature, formula = formula, dist = dist, link = link, ...)
  #assays(input)$counts_str0 <- zi_result$ziOutput
  #assays(input)$weights <- zi_result$weights
  result <- zi(input, zi_result$ziModel, zi_result$ziOutput, zi_result$weights)
  return(result)
}

