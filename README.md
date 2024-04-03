# zitools

<!-- badges: start -->

<!-- badges: end -->

One fundamental objective of microbiome analyses is to identify differentially abundant taxa among different experimental
groups or conditions. However, microbiome data are often overdispersed and zero inflated making data analysis extremely
challenging. Although there are several models considering zero inflation, none of them provides functionality for
subsequent analyses. Therefore, we propose zitools, an R package allowing for zero inflated count data analysis by
either using down-weighting of excess zeros or by replacing an appropriate proportion of excess zeros with NA. Through
overloading frequently used statistical functions (such as mean, median, standard deviation), plotting functions (such as
boxplots or heatmap) or differential abundance tests, it allows a wide range of downstream analyses for zero-inflated data
in a less biased manner.


## Installation

You can install the development version of zitools from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kreutz-lab/zitools")
```
``` r
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("zitools")
``` 

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(zitools)
## basic example code
```
