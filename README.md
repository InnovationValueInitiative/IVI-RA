# The IVI-RA Individual Patient Simulation Model
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
* [R package](https://innovationvalueinitiative.github.io/IVI-RA/)
* [Model description](https://innovationvalueinitiative.github.io/IVI-RA/model-description/model-description.pdf)

## Web interfaces
You can also run the model online using our web interfaces.

* [IVI-RA Model Interface](https://innovationandvalueinitiative.shinyapps.io/ivi-ra-expert/): full control over treatment sequences, the patient population, model parameters, model structures, and the time horizon. 
* [IVI-RA Value Tool](http://apps.thevalueinitiative.org/ivi-ra/): a more streamlined experience for users with less experience in decision-analytic modeling and rheumatoid arthritis. 

## Collaborate
* If you found a bug or would like to request a new feature, then submit an [issue](https://github.com/InnovationValueInitiative/IVI-RA/issues).
* If you would like to contribute to the code base, then create a pull request. For more about contributing to the code, see [here](https://innovationvalueinitiative.github.io/IVI-RA/articles/how-to-contribute.html). 






