library("usethis")
library("devtools")


#git repository
use_git(message = "Initial commit")
?use_git


#extract matrix from a given object - either phyloseq object or RangedSummarizedExperiment
getmtx <- function(x, y) {
  if(missing(y)) {
    if (class(x) == "phyloseq") {
      df <- data.frame(otu_table(x))
    }
    else {
      if (class(x) == "RangedSummarizedExperiment") {
        df <- data.frame(assays(x)$y)
      }
    }
  }
  return(df)
}
#given matrix - reshape it into long format as input for zeroinflation model, df_wide = given matrix, feature = rows, sample = columns, feature = "OTU", "Gene", "Species", etc
reshape_zi <- function(df_wide, feature = "") {
  zi_long <- data.frame(df_wide) %>%
    rownames_to_column(var = feature)%>%
    gather(key = "sample", value = "count", 1:ncol(df_wide)+1)
  return(zi_long)
}

#fit zero inflation model - zeroinfl() function of the pscl package
pscl::zeroinfl(formula, data, subset, na.action, weights, offset,
               dist = c("poisson", "negbin", "geometric"),
               link = c("logit", "probit", "cloglog", "cauchit", "log"),
               control = zeroinfl.control(...),
               model = TRUE, y = TRUE, x = FALSE)

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

#function to calculate weights (zinbwave), given phyloseq object
calcWeights <- function(x, y, ...) {
  ps_se <- phyloseq_to_deseq2(x, y)
  zinbwave <- zinbwave(Y = ps_se, observationalWeights = TRUE,
                       BPPARAM = BiocParallel::SerialParam(), ...)
  weigths <- assay(zinbwave, "weights")
  return(weights)
}

#zi_main <- function(df_wide, feature,  formula,  dist = c("poisson", "negbin", "geometric"), link = c("logit", "probit", "cloglog", "cauchit", "log"))
zi_main <- function(df_wide, feature = "",  formula, dist = c("poisson", "negbin", "geometric"), link = c("logit", "probit", "cloglog", "cauchit", "log"),ps, group)
{
  zi_input <- reshape_zi(df_wide, feature)
  zi <- zeroinfl(formula, data = zi_input, dist = dist,
                 link = link)
  zi_prediction_long <- omit_str_zero(zi, zi_input, feature)
  zi_prediction_wide <- zi_prediction_long %>%
    spread(key = "sample", value = "count") %>%
    column_to_rownames(var = feature)
  weights <- calcWeights(ps, group)
  result <- zi(df_wide, zi, zi_prediction_wide, weights)
  return(result)
}



