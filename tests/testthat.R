# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(zitools)

# test_check("zitools")

# Unit tests for comparing the case where all weights are 1 to the case where zi@inputcounts is used

data(mtx)
Zi <- ziMain(mtx)

ZiW1 <- Zi # here all weights will be 1
ZiW1@weights <- Zi@weights * 0 + 1

# Testing whether the results form generic sd and zitools::sd are the same:
# test_that("Check results for weights equals to one", {
#     expect_identical(zitools::sd(ZiW1), sd(ZiW1@inputcounts))
#     expect_identical(zitools::mean(ZiW1), mean(ZiW1@inputcounts))
# })

##
sd(ZiW1)- sd(ZiW1@inputcounts)
mean(ZiW1)- mean(ZiW1@inputcounts)
sum(abs(rowMeans2(ZiW1)- rowMeans2(ZiW1@inputcounts)))
sum(abs(colMeans2(ZiW1)- colMeans2(ZiW1@inputcounts)))

## Testing whether weights =1/2 lead correct results
mtx2 <- rbind(mtx,mtx)
Zi2 <- ziMain(mtx2)
Zi2@weights <- rbind(Zi@weights,Zi@weights)/2

sum(abs(colMeans2(Zi2)- colMeans2(Zi)))
mean(Zi2) - mean(Zi)
sd(Zi2) - sd(Zi)


## Testing whether counts>0 lead correct results
mtxAll <- mtx
mtxAll[mtxAll==0] <- 1
ZiAll <- ziMain(mtxAll)
sum(ZiAll@weights<1)


