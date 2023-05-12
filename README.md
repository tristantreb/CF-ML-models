# Machine learning models for CF

Project 1: Bayesian inference with Expectation Maximisation for the characterisation of antibiotic treatment recovery in Cystic Fibrosis

Project 2: Estimation of the variability in Cystic Fibrosis patient’s FEV1 lung function measurements

Hosted at the University of Cambridge

## Documentation
The code documentation is available [here](https://tristantreb.github.io/CF-ML-models/).

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
1. Clone this repo in a folder called `Code/`
2. In that same folder, clone [this](https://github.com/damiansphd/smartcare/tree/9179127f15e96085db14e362f9d43b36f488472e) version of the repo - it shows the project version at commit 9179127f15e96085db14e362f9d43b36f488472e (26.07.2021)
3. You should obtain the following folder structure `Code/master_thesis_CF_ML` and `Code/smartcare`
4. Change the absolute pathes in the init.m files - see `FEV1variability`, `exploration`,  `recovery`
5. You are ready!

### Main project
1. Run `masterscriptRecovery` to load the data 
2. Explore data with functions in `exploration`
3. For the main project: run an alignment model with `runAlignmentModelEMMCRecoveryFcn`

### Side project
1. Run `analyseFEV1Variability`

## Data
The data is not available as it is stored in a secured server due to privacy reasons.

## Contact
Contact me at [tristan.trebaol@alumni.epfl.ch](mailto:aiyu@di.ku.dk) if you would like to know more about the project!
