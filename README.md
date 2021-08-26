# MSc Thesis

Thesis title: Bayesian inference with Expectation Maximisation for the characterisation of antibiotic treatment recovery in Cystic Fibrosis

Side project title: Estimation of the variability in Cystic Fibrosis patient’s FEV1 lung function measurements

Hosted at the University of Cambridge, supervised at EPFL

## Documentation
The code documentation is available [here](https://tristantreb.github.io/master_thesis_CF_ML/).

## Repository structure
```
    ├── FEV1variability               <- code base for the side project
    ├── docs                          <- code documentation
    ├── exploration                   <- code base for exploration of the data
    ├── helperfunctions               <- helper functions common to at least two other matlab scripts
    ├── msc-tristan                   <- important documents necessary to reproduce all results
        ├── report                      <- thesis report with LaTeX version, thesis defense presentation
    ├── recovery                      <- code base for the main project
        ├── updatedModel                <- updated version of the ML model
```

## Code base usage
1. Run `masterscriptRecovery` to load the data 
2. Explore data with functions in `exploration`
3. For the main project: run an alignment model with `runAlignmentModelEMMCRecoveryFcn`
For the side project: run `analyseFEV1Variability`

more details in the documentation

## Data
The data is not available as it is stored in a secured server due to privacy reasons.

## Contact
Contact me at [tristan.trebaol@alumni.epfl.ch](mailto:aiyu@di.ku.dk) if you would like to know more about the project!
