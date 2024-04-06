#'@include ziMain.R
NULL


#'@name reshape_zi
#'@title Reshape a given matrix into long format
#'
#'@description reshape a given matrix into a long format as input for for
#'fitting a zero inflated model
#'
#'@param mtx count matrix
#'@param feature character string characterizing the rows, e.g. gene, OTU
#'
#'@returns a dataframe (long format), 3 columns: count, sample, 'feature'
#'
#'@importFrom tibble rownames_to_column
#'@importFrom magrittr %>%
#'@importFrom tidyr gather
#'@noRd
reshape_zi <- function(mtx, feature = "") {
    zi_long <- as.data.frame(mtx) %>%
        rownames_to_column(var = feature) %>%
        gather(key = "sample", value = "count", 1:ncol(mtx) + 1)
    return(zi_long)
}

#'@name omit_str_zero
#'@title Draw Structural Zeros
#'
#'@description draw structural zeros with a given probability (p_str_zero),
#'compare if predicted zero = zero count, if YES set to NA, repeat until
#'NA >= sum(p_str_zero), replace excess of NA with zero (random sampling)
#'
#'@param zi result of fitting a zero inflation model using pscl::zeroinfl
#'@param zi_input result of reshape_zi, count data in a long format
#'@param feature character string characterizing the rows, e.g. gene, OTU, ...
#'
#'@returns a dataframe(long format), columns: count, sample, 'feature',
#'
#'@importFrom magrittr %>%
#'@importFrom dplyr bind_cols
#'@importFrom dplyr mutate
#'@importFrom stats rbinom
#'@importFrom dplyr filter
#'@importFrom dplyr sample_n
#'@importFrom tidyr replace_na
#'@import pscl
#'@noRd

omit_str_zero <- function(zi, zi_input, feature = "") {
    zi_prediction <- data.frame(predicted_count =
        predict(zi, type = "count")) %>%
        bind_cols(data.frame(p_str_zero = predict(zi, type = "zero"))) %>%
        bind_cols(zi_input, .)
    str_zero <- bind_cols(zi_prediction, predicted_zero =
        c(rbinom(nrow(zi_prediction), 1, prob = zi_prediction$p_str_zero)))
    repeat {
        str_zero <- str_zero %>%
        mutate(predicted_zero = replace(predicted_zero, 1:nrow(str_zero),
            c(rbinom(nrow(str_zero), 1, prob = str_zero$p_str_zero))))
        str_zero$count <- ifelse(str_zero$count == 0 &
            str_zero$predicted_zero == 1, NA, str_zero$count)
        str_zero$na <- ifelse(is.na(str_zero$count), 1, 0)
    if (sum(str_zero$na) >= round(sum(str_zero$p_str_zero), digits = 0)) {
        break
        }
    }
    replaced_zero <- str_zero %>%
        filter(is.na(count)) %>%
        sample_n(sum(str_zero$na) - round(sum(str_zero$p_str_zero)),
            replace = FALSE) %>%
        replace_na(list(count = 0))
    zi_replaced <- merge(str_zero, replaced_zero, by.x = 1:2, by.y = 1:2,
        all = TRUE)
    zi_replaced <- cbind(zi_replaced[c(1, 2)], count = with(zi_replaced,
        ifelse(is.na(count.y), count.x, count.y)))

    return(zi_replaced)
}

#'@name calcWeights
#'@title calculate weights
#'
#'@description the posterior probability that a given count y arises from the
#'NB count component can be used as weights and is calculated using the Bayes'
#'rule: w = ((1-pi)*fnb)/fzinb
#'
#'@param zi result of fitting a zero inflation model to the data using
#'pscl::zeroinfl
#'@param zi_input result of reshape_zi, count data in a long format
#'@param feature character string characterizing the rows, e.g. gene, OTU, ...
#'
#'@returns matrix of the calculated weights
#'
#'@import pscl
#'@importFrom dplyr mutate
#'@importFrom dplyr case_when
#'@importFrom stats dpois
#'@importFrom VGAM dzipois
#'@importFrom stats dnbinom
#'@importFrom VGAM dzinegbin
#'@importFrom dplyr select
#'@importFrom tidyr spread
#'@importFrom tibble column_to_rownames
#'@importFrom magrittr %>%
#'@noRd
#'
# function to calculate weights - formula: w=(1-pi)*fnb/fzinb
calcWeights <- function(zi_input, zi, feature, dist) {
    df <- data.frame(zi_input, pstr0 = predict(zi, type = "zero"))
    df <- mutate(df, pstr0 = case_when(count != 0 ~ 0, TRUE ~ pstr0))
    if (dist == "poisson") {
    f <- dpois(x = zi_input$count, lambda = predict(zi, type = "count"),
        log = FALSE)
    f_zi <- dzipois(x = zi_input$count, lambda = predict(zi, type = "count"),
        pstr0 = df$pstr0, log = FALSE)
    }
    if (dist == "negbin") {
        f <- dnbinom(x = zi_input$count, size = zi[["theta"]], mu = predict(zi,
            type = "count"), log = FALSE)
        f_zi <- dzinegbin(x = zi_input$count, size = zi[["theta"]],
            munb = predict(zi, type = "count"), pstr0 = df$pstr0, log = FALSE)
    }
    df$f <- f
    df$f_zi <- f_zi
    df$weights <- c(((1 - df$pstr0) * df$f)/df$f_zi)
    df$"1-pi" <- c(1 - df$pstr0)
    df_wide <- df %>%
        select(c(feature, "sample", "weights")) %>%
        spread(key = "sample", value = "weights") %>%
        column_to_rownames(var = feature) %>%
        as.matrix(.)
    return(df_wide)
}


#'@name preprocess_mtx
#'@title Preprocess Matrix
#'
#'@description Preprocess matrix before subsetting it into blocks. The matrix is
#'randomly sorted to make sure that no prior sorting interferes when fitting the
#'model
#'
#'@param mtx count matrix
#'
#'@returns randomly sorted matrix
#'@noRd
preprocess_mtx <- function(mtx) {
    set.seed(123)
    mtx <- mtx[rowSums(mtx[]) > 0, ]
    seed <- .Random.seed
    set.seed(seed)
    random_rows <- sample(nrow(mtx))
    mtx <- mtx[random_rows, ]
    .Random.seed <- seed
    return(mtx)
}

#'@name subset_mtx
#'@title subset matrix into groups of roughly 5000 datapoints per block
#'
#'@description subset a matrix into blocks of approximately 5000 datapoints per
#'block. First, the number of the blocks is calculated (n_blocks). Then, the
#'number of rows for each block is calculated and based on that an index for the
#'blocks is created. According to the index, the matrix is split into the blocks
#'and after splitting the matrix, the index is removed.
#'
#'@param mtx matrix to be subsetted
#'
#'@returns list containing the subsets
#'@noRd
subset_mtx <- function(mtx) {
    n_blocks <- round(nrow(mtx) * ncol(mtx)/5000)
    if (n_blocks == 0) {
    n_blocks <- 1
    }
    n_rowsperblock <- ceiling(nrow(mtx)/n_blocks)
    index <- rep(c(1:n_blocks), each = n_rowsperblock)
    index <- index[1:nrow(mtx)]
    df <- as.data.frame(mtx)
    df$index <- index
    list_subset <- split(df, index)
    list_subset <- lapply(list_subset, function(x) x[!(names(x) %in%
        c("index"))])
    return(list_subset)
}

#'@name zi_core
#'
#'
#'@description core function that is applied to the subsets of the count matrix
#'
#'@param input matrix
#'@param feature character string characterizing the rows, e.g. gene, OTU, ...
#'@param dist character specification of count model family, either poisson,
#'negbin, geometric (a log link is always used)
#'@param link character specification of link function in the binary zero-
#'inflation model, either logit, probit, cloglog, cauchit, log (a binomial
#'family is always used).
#'
#'@returns list including the input matrix, the zero inflation model, the matrix
#'where predicted structural zeros are replaced with NA, and the weight matrix
#'
#'@importFrom pscl zeroinfl
#'@importFrom tidyr spread
#'@importFrom tibble column_to_rownames
#'@importFrom magrittr %>%
#'@noRd
#'
zi_core <- function(input, feature = "", formula, dist, link, ...) {
    matrix <- as.matrix(input)
    zi_input <- reshape_zi(input, feature)
    zi <- pscl::zeroinfl(formula, data = zi_input, dist = dist, link = link,
        ...)
    zi_prediction_long <- omit_str_zero(zi, zi_input, feature)
    zi_prediction_wide <- zi_prediction_long %>%
        spread(key = "sample", value = "count") %>%
        column_to_rownames(var = feature) %>%
        as.matrix()
    weights <- calcWeights(zi_input, zi, feature, dist = dist)
    result <- list(ziInput = matrix, model = zi,
        zideinflatedcounts = zi_prediction_wide, weights = weights)
    return(result)
}
