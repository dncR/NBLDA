## Introduction

[![CRAN_version](http://www.r-pkg.org/badges/version/NBLDA)](https://cran.r-project.org/package=NBLDA)
[![CRAN_download](http://cranlogs.r-pkg.org/badges/NBLDA)](https://cran.r-project.org/package=NBLDA)
[![CRAN_download_total](http://cranlogs.r-pkg.org/badges/grand-total/NBLDA)](https://cran.r-project.org/package=NBLDA)


<br>

We proposed a package for classification task which uses Negative Binomial distribution within Linear Discriminant Analysis. It is basically an extension of PoiClaClu package to Negative Binomial distribution. The classification algorithms are based on the papers Dong et al. (2016, ISSN: 1471-2105) and Witten, DM (2011, ISSN: 1932-6157) for NBLDA and PLDA respectively. Although PLDA is a sparse algorithm and can be used for variable selection, the algorithm proposed by Dong et. al. is not sparse, hence, it uses all variables in the classifier. Here, we extent Dong et. al.'s algorithm to sparse case by shrinking overdispersion towards 0 (Yu et. al., 2013, ISSN: 1367-4803) and offset parameter towards 1 (as proposed by Witten DM, 2011). We support only the classification task with this version. However, the clustering task will be included with the following versions.

To install the NBLDA package in R:

```{r, eval = FALSE, message=FALSE, warning=FALSE}
install.packages("NBLDA")
```

To install the developement version of NBLDA package from GitHub:

```{r, eval = FALSE, message=FALSE, warning=FALSE}
devtools::install_github("dncR/NBLDA")
```

If you use NBLDA package in your research, please cite it as below:

> Dincer Goksuluk, Gokmen Zararsiz, Selcuk Korkmaz and Ahmet Ergun Karaagaoglu (<YEAR_here>). NBLDA: Negative Binomial Linear Discriminant Analysis. R package
  version "use_package_version_here".


To get BibTeX entry for LaTeX users, type the following:

```{r, eval = FALSE}
citation("NBLDA")
```

<br>

Please contact us, if you have any questions or suggestions:

  dincergoksuluk@erciyes.edu.tr

## News:

#### First Release (1.0)
* Poisson (PLDA) and Negative Binomail Linear Discriminant (NBLDA) functions included.
* Sparse version of NBLDA is proposed.
* Vignette is not available. It will be released with next version.
