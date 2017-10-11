# Overview
`iviRA` is an R package that runs the [Innovation and Value Initiative's (IVIs)](http://www.thevalueinitiative.org/) individual patient simulation model for rheumatoid arthritis. The IVI-RA model is an iterative open-source model that is updated according to evidence-based suggestions for improvement. Information on the collaborative and consensus-based process can be found [here](articles/ways-to-contribute.html).

The model simulates the costs, health outcomes, and risks associated with disease-modifying anti-rheumatic drugs (DMARDs) and biologic DMARDSs (bDMARDs) for patients with moderate to severe rheumatoid arthritis. The model contains 384 possible model structures, which can be used to quantify structural uncertainty or to evaluate the implications of different modeling assumptions. Parameter uncertainty is quantified using probabilistic sensitivity analysis (PSA). The simulation is mostly written in C++ so that it can be run in a reasonable amount of time. For a detailed explanation of the model see the [model description](model-description/model-description.pdf).

In addition to running the model with the R package, users can run the model online with extensive control over settings using the [IVI-RA model web interface](https://innovationandvalueinitiative.shinyapps.io/ivi-ra-expert/). Furthermore, stakeholders without expertise in RA modeling or value assessment can run the model with the [IVI-RA Value Tool](https://innovationandvalueinitiative.shinyapps.io/ivi-ra).

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