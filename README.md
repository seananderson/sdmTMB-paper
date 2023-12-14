## sdmTMB paper

Source code for

Anderson, S.C., E.J. Ward, P.A. English, L.A.K. Barnett, J.T. Thorson. sdmTMB: an R package for fast, flexible, and user-friendly generalized linear mixed effects models with spatial and spatiotemporal random fields.

A preprint of an earlier version is available at:

Anderson, S.C., E.J. Ward, P.A. English, L.A.K. Barnett. 2022. sdmTMB: an R package for fast, flexible, and user-friendly generalized linear mixed effects models with spatial and spatiotemporal random fields. bioRxiv 2022.03.24.485545; doi: <https://doi.org/10.1101/2022.03.24.485545>

sdmTMB:

<https://github.com/pbs-assess/sdmTMB>\
<https://CRAN.R-project.org/package=sdmTMB>

The paper source code is in `doc/paper-jss/sdmTMB-paper.Rnw`.

The paper can be built with:

```r
setwd("doc/paper-jss/")
knitr::knit("sdmTMB-paper.Rnw")
tinytex::latexmk("sdmTMB-paper.tex", engine = "xelatex")
```
