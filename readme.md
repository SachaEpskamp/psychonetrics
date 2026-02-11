# psychonetrics

**Structural Equation Modeling and Confirmatory Network Analysis**

The `psychonetrics` R package provides a comprehensive framework for multi-group structural equation modeling (SEM) combined with confirmatory network analysis. It supports cross-sectional, time-series, and panel data through maximum likelihood estimation with core computations in C++ for fast performance.

## Modeling Frameworks

- **Variance-covariance modeling** (`varcov`) -- Gaussian graphical model (GGM) via `ggm()`
- **Latent variable modeling** (`lvm`) -- Confirmatory factor analysis via `lvm()`, latent network model via `lnm()`, residual network model via `rnm()`
- **Ising model** (`Ising`) -- Binary data network modeling via `Ising()`
- **Graphical VAR** -- Time-series modeling via `gvar()`

## Installation

The package can be installed from GitHub:

```r
library("devtools")
install_github("SachaEpskamp/psychonetrics")
```

## Documentation

Full documentation, tutorials, and examples are available at the package website:

**[https://psychonetrics.org](https://psychonetrics.org)**

## Key References

Epskamp, S. (2024). *psychonetrics: Structural Equation Modeling and Confirmatory Network Analysis.* R package. [doi:10.31234/osf.io/8ha93](https://doi.org/10.31234/osf.io/8ha93)

Epskamp, S., Rhemtulla, M.T., & Borsboom, D. (2017). Generalized Network Psychometrics: Combining Network and Latent Variable Models. *Psychometrika, 82*(4), 904-927. [doi:10.1007/s11336-017-9557-x](https://doi.org/10.1007/s11336-017-9557-x)

## Bug Reports

Please report bugs and feature requests at [https://github.com/SachaEpskamp/psychonetrics/issues](https://github.com/SachaEpskamp/psychonetrics/issues).
