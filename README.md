<!--
[![Build Status](https://travis-ci.org/dmlc/xgboost.svg?branch=master)](https://travis-ci.org/dmlc/xgboost)
[![Documentation Status](https://readthedocs.org/projects/xgboost/badge/?version=latest)](https://xgboost.readthedocs.org)
[![GitHub license](http://dmlc.github.io/img/apache2.svg)](./LICENSE)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/xgboost)](http://cran.r-project.org/web/packages/xgboost)
[![PyPI version](https://badge.fury.io/py/xgboost.svg)](https://pypi.python.org/pypi/xgboost/)
[![Gitter chat for developers at https://gitter.im/dmlc/xgboost](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/dmlc/xgboost?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[Documentation](https://xgboost.readthedocs.org) |
[Resources](demo/README.md) |
[Installation](https://xgboost.readthedocs.org/en/latest/build.html) |
[Release Notes](NEWS.md) |
[RoadMap](https://github.com/dmlc/xgboost/issues/873)
-->

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

This is an R packages to generate genomic and bioinformatic pipelines and submit jobs to a HPC running slurm system.

## Install

Install [devtools](https://github.com/hadley/devtools) first, and then use `devtools` to install `huskeR` from github.

```R
#install.packages(devtools)
devtools::install_github("yangjl/huskeR")

# OR specify a path to install
library(devtools)
withr::with_libpaths(new = "/home/jyanglab/jyang21/R/x86_64-pc-linux-gnu-library/3.5", install_github('yangjl/huskeR'))

library(huskeR)
```

List all the functions in the package and find help.

```R
## list all the functions
ls(getNamespace("huskeR"), all.names=TRUE)
## help info for a given function
?run_fastq_qc
```



## License

It is a free and open source software, licensed under [GPLv3](LICENSE).
This is an ongoing research project from **Yang Lab**. It was intended for internal lab usage. It has not been extensively tested. Use at your own risk.


