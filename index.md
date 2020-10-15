[![Travis-CI Build Status](https://travis-ci.org/InnovationValueInitiative/IVI-RA.svg?branch=master)](https://travis-ci.org/InnovationValueInitiative/IVI-RA)
[![Coverage Status](https://codecov.io/gh/InnovationValueInitiative/IVI-RA/branch/master/graph/badge.svg)](https://codecov.io/gh/InnovationValueInitiative/IVI-RA)

# Overview
`iviRA` is an R package that runs the [Innovation and Value Initiative's (IVI's)](http://www.thevalueinitiative.org/) individual patient simulation model for rheumatoid arthritis (RA) (the IVI-RA model). The model simulates the costs, health outcomes, and risks associated with disease-modifying anti-rheumatic drugs (DMARDs) including conventional DMARDs (cDMARDs), biologic DMARDs (bDMARDs), and Janus kinase/signal transducers and activators of transcription (JAK/STAT) inhibitors for patients with moderate to severe rheumatoid arthritis (RA). The model is intended to help decision-makers assess the value of treatments for a population of patients with RA. 

## Installation
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

## Documentation
* [Model description](model-description/model-description.pdf)
* [iviRA tutorial](articles/00-intro.html)
* [iviRA API](reference/index.html)


## Collaborate
The IVI-RA model is part of the [Open Source Value Project (OSVP)](http://www.thevalueinitiative.org/open-source-value-project/), a consensus-based process for the development of open-source cost-effectiveness models and other tools for value assessment of medical interventions. Learn more about how to collaborate [here](articles/how-to-contribute.html).

## Web applications
In addition to running the model with the R package, users can run the model online with our web interfaces:

* [IVI-RA Model Interface](https://ivi-ra-expert.clarityviz.com/): full control over treatment sequences, the patient population, model parameters, model structures, and the time horizon.
* [IVI-RA Value Tool](https://ivi-ra.clarityviz.com/): a more streamlined experience for users with less experience in decision-analytic modeling and rheumatoid arthritis.
