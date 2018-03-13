# Installation

In R, install the package "devtools" if you haven't already, and use `install_github` to get this package. 
```r
install.packages("devtools")
library("devtools")
install_github("jfukuyama/phyloseqGraphTest")
```

# Usage

This package is meant to be used with phyloseq, so suppose that you have a phyloseq object called ps. The sample data in the phyloseq object has a column corresponding to the groups you want to test, called testVar. This should be given as a character string. If applicable, you can also specify a stratification variable (like the strata argument in adonis) either by its name or as a vector with the same number of entries as samples in ps. So for example, to perform the test you would use
```r
gt = graph_perm_test(ps, sampletype = "testVar", grouping = "stratVar",
                     distance = "jaccard", type = "mst")
```
Then the p-value is stored in `gt$pval`, and you can plot the network constructed and the permutation distribution (with the observed test statistic marked in red) with the functions `plot_test_network` and `plot_permutations`.
```r
plot_test_network(gt)
plot_permutations(gt)
```

[![Travis-CI Build Status](https://travis-ci.org/jfukuyama/phyloseqGraphTest.svg?branch=master)](https://travis-ci.org/jfukuyama/phyloseqGraphTest)
