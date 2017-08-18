# Overview
`iviRA` is an R package that runs the [Innovation and Value Initiative's (IVIs)](http://www.thevalueinitiative.org/) individual patient simulation model for rheumatoid arthritis. The IVI-RA model is an iterative open-source model that is updated according to evidence-based suggestions for improvement. Information on the collaborative and consensus-based process can be found [here](articles/ways-to-contribute.html).

The model can be run with multiple perspectives (e.g., health care sector, societal, healthy or sick population) and accounts for both structural and parameter uncertainty. Structural uncertainty can be quantified by simulating treatment value across our 336 possible model structures and parameter uncertainty is quantified using probabilistic sensitivity analysis (PSA). The simulation is mostly written in C++ so that it can be run in a reasonable amount of time.

For a detailed explanation of the model see the [model description](model-description/model-description.pdf). To run the model with our interactive user-interface, use our [R Shiny web application](https://innovationandvalueinitiative.shinyapps.io/ivi-ra/).

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