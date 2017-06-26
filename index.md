# Overview
`iviRA` is an R package that runs the [Innovation and Value Initiative's (IVIs)](http://www.thevalueinitiative.org/) cost-effectiveness model for rheumatoid arthritis. The IVI-RA model is an individual patient simulation (IPS) written mostly in C++ so that it can be run in a reasonable amount of time. The model is flexible and reflects a range of plausible structural assumptions and perspectives (e.g. health care sector, societal). Parameter uncertainty is quantified using probabilistic sensitivity analysis (PSA).

For a detailed explanation of the model see the [model description](model-description/model-description.pdf).

# Installation
`iviRA` can be installed from GitHub using `devtools`:

```r
# install.packages("devtools")
library(devtools)
devtools::install_github("InnovationValueInitiative/IVI-RA")
```

It can then be loaded into `R`:

```r
library(iviRA)
```