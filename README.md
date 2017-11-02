# The IVI-RA Individual Patient Simulation Model
`iviRA` is an R package that runs the [Innovation and Value Initiative's (IVI's)](http://www.thevalueinitiative.org/) individual patient simulation model for rheumatoid arthritis (RA) (the IVI-RA model). The model simulates the costs, health outcomes, and risks associated with disease-modifying anti-rheumatic drugs (DMARDs) including conventional DMARDs (cDMARDs), biologic DMARDs (bDMARDs), and Janus kinase/STAT pathway inhibitors for patients with moderate to severe rheumatoid arthritis (RA). The model is intended to help decision-makers assess the value of treatments for a population of patients with RA.

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
* [R package](https://innovationvalueinitiative.github.io/IVI-RA/)
* [Model description](https://innovationvalueinitiative.github.io/IVI-RA/model-description/model-description.pdf)

## Collaborate
* If you found a bug or would like to request a new feature, then submit an [issue](https://github.com/InnovationValueInitiative/IVI-RA/issues).
* If you would like to contribute to the code base, then create a pull request. We, for the most part, adhere to Hadley Wickham's [style guide](http://adv-r.had.co.nz/Style.html). The one exception is that we use dots to separate object names as in my.name.





