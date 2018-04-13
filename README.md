
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis-CI Build Status](https://travis-ci.org/SMAC-Group/classimu.svg?branch=master)](https://travis-ci.org/SMAC-Group/classimu) [![Project Status: Active](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![Licence](https://img.shields.io/badge/licence-CC%20BY--NC--SA%204.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html) [![minimal R version](https://img.shields.io/badge/R%3E%3D-3.4.0-6666ff.svg)](https://cran.r-project.org/) [![CRAN](http://www.r-pkg.org/badges/version/classimu)](https://cran.r-project.org/package=classimu) [![packageversion](https://img.shields.io/badge/Package%20version-0.1.0-orange.svg?style=flat-square)](commits/develop) [![Last-changedate](https://img.shields.io/badge/last%20change-2018--01--15-yellowgreen.svg)](/commits/master)

`mgmwm` Overview <a href="https://smac-group.com/"><img src="man/figures/logo.png" align="right" style="width: 20%; height: 20%"/></a>
=========================================================================================================================================

The Multivariate GMWM (`mgmwm`) R package estimates and select time series models:


To see what `mgmwm` is capable of, please refer to the "Vignettes" tabs above.

Install Instructions
--------------------

To install the `mgmwm` package, there is currently one option: [GitHub](https://github.com/SMAC-Group/classimu/).

### Installing the package through GitHub

For users who are interested in having the latest developments, this option is ideal. Though, more dependancies are required to run a stable version of the package. Most importantly, users **must** have a compiler installed on their machine that is compatible with R (e.g. Clang).

*The setup to obtain the development version of `classimu` is platform dependent.*

### Requirements and Dependencies

**OS X**

Some users report the need to use X11 to suppress shared library errors. To install X11, visit [xquartz.org](http://www.xquartz.org/).

**Linux**

Both curl and libxml are required.

For **Debian** systems, enter the following in terminal:

``` bash
sudo apt-get install curl libcurl3 libcurl3-dev libxml2 libxml2-dev
```

For **RHEL** systems, enter the following in terminal:

``` bash
sudo yum install curl curl-devel libxml2 libxml2-dev
```

**All Systems**

The following R packages are also required. If you have made it this far, run the following code in an R session and you will be ready to use the devlopment version of `mgmwm`.

``` r
# Install dependencies
install.packages(c("RcppArmadillo","devtools","knitr","rmarkdown"))

# Install the package from GitHub without Vignettes/User Guides
devtools::install_github("SMAC-Group/mgmwm")

# Install the package with Vignettes/User Guides 
devtools::install_github("SMAC-Group/mgmwm", build_vignettes = TRUE)
```
