Psychonetrics
================
Sacha Epskamp

The R package `psychonetrics` is designed to be a package for confirmatory multi-group network analysis. The goal of the package is to provide relatively fast maximum likelihood estimators for the following modeling frameworks:

-   Gaussian graphical model - `ggm()`
-   Residual network model- `rnm()`
-   Latent network model - `lnm()`
-   Graphical vector auto-regression for *n* = 1 time-series - `gvar()`
-   Graphical vector auto-regression for panel data - Not yet implemented
-   Ising model - Not yet implemented
-   Fused latent and graphical IRT (Chen et al. 2018) - Not yet implemented

Note that this version is highly unstable, and that all reported values such as standard errors, *p*-values and modification indices, have not yet been thoroughly validated (please let me know if you wish to help out with this).

Installation
============

The package can be installed from Github:

``` r
library("devtools")
install_github("sachaepskamp/psychonetrics")
```

and subsequently loaded as usual. I will also load the `dplyr` package to obtain the useful pipe operator `%>%`:

``` r
library("psychonetrics")
library("dplyr")
```

Gaussian graphical model
========================

Let's take the `bfi` dataset, but leave the first 1000 observations out to later test our model on:

``` r
library("psych")
# Load data:
data(bfi)

# Define variables:
vars <- names(bfi)[1:25]

# Let's use the first 1000 rows as test data:
testData <- bfi[1:1000, ]

# And the rest as training data:
trainingData <- bfi[-(1:1000), ]
```

This data does not have a lot of missing data:

``` r
mean(is.na(trainingData))
```

    ## [1] 0.009265873

So we can use the faster maximum likelihood estimation (especially given this large sample size), which uses listwise missing data deletion (set `estimator = "FIML"` for full information maximum likelihood instead).

The first step in using psychonetrics is to form a psychonetrics model:

``` r
model <- ggm(trainingData, vars = vars, omega = "full")
```

This creates an *unevaluated* model. The `omega` argument corresponds with the partial correlation matrix, using the following modeling framework for the GGM (Epskamp, Rhemtulla, and Borsboom 2017):

-   **μ** = **μ**
-   **Σ** = **Δ**(**I**−**Ω**)<sup>−1</sup>**Δ**

in which **Σ** represents the modeled variance--covariance matrix and **μ** the mean structure. The diagonal matrix **Δ** is a scaling matrix that is added by default (argument `delta`). This framework is especially useful in multi-group analysis, as it allows to fix partial correlations in **Ω** across groups while keeping the scaling unconstrained. We can use the `omega` argument in the `ggm` function to specify a network structure by assigning it an adjacency matrix (a matrix of zeroes and ones). The defaults `model = "full"` and `model = "empty"` instead simply create a fully connected network model or an empty network model.

The psychonetrics model is a rich model that contains a lot of information. For example, it stores the required summary statistics and fit functions (e.g,. gradient and Hessian), and even keeps a logbook of all changes made including time-stamps and `sessionInfo` information. Many other functions in psychonetrics will update the model in some way, or use the model to obtain some information to print.

The model, however, has not yet been evaluated. to do this, we need to apply the `runmodels` function:

``` r
model <- model %>% runmodel
```

    ## Estimating baseline model...

    ## Estimating model...

    ## Computing modification indices...

This will compute the baseline and saturated model as well, which are stored in the object so they don't have to be computed again. We can now print the model to obtain some information:

``` r
model
```

    ##                         _                      _        _          
    ##                        | |                    | |      (_)         
    ##    _ __  ___ _   _  ___| |__   ___  _ __   ___| |_ _ __ _  ___ ___ 
    ##   |  _ \/ __| | | |/ __|  _ \ / _ \|  _ \ / _ \ __|  __| |/ __/ __|
    ##   | |_) \__ \ |_| | (__| | | | (_) | | | |  __/ |_| |  | | (__\__ \
    ##   | .__/|___/\__, |\___|_| |_|\___/|_| |_|\___|\__|_|  |_|\___|___/
    ##   | |         __/ |                                                
    ##   |_|        |___/                                                 
    ##  
    ## 
    ## General: 
    ##  - psychonetrics version: 0.1.2 
    ##  - Model last edited at: 2019-03-05 21:46:17
    ## 
    ## Sample: 
    ##  - Number of cases: 1562 
    ##  - Number of groups: 1 
    ##  - Number of observed summary statistics: 350
    ## 
    ## Model: 
    ##  - model used: Gaussian graphical model (GGM) 
    ##  - Number of parameters: 350
    ## 
    ## Estimation: 
    ##  - Optimizer used: nlminb 
    ##  - Estimator used: Maximum likelihood estimation (ML) 
    ##  - Message: both X-convergence and relative convergence (5)
    ## 
    ## Fit: 
    ##  - Model Fit Test Statistic: < 0.0001 
    ##  - Degrees of freedom: 0 
    ##  - p-value (Chi-square): 1
    ## 
    ## Tips: 
    ##  - Use 'psychonetrics::compare' to compare psychonetrics models 
    ##  - Use 'psychonetrics::fit' to inspect model fit 
    ##  - Use 'psychonetrics::parameters' to inspect model parameters 
    ##  - Use 'psychonetrics::MIs' to inspect modification indices

Since the model is saturated, this is not that interesting. Let's check the estimated parameters:

``` r
model %>% parameters
```

    ## 
    ##  Parameters for group 1
    ##  -  mu  
    ##  var1 op var2  est    se        p row col par
    ##    A1 ~1      2.46 0.036 < 0.0001   1   1   1
    ##    A2 ~1      4.75 0.031 < 0.0001   2   1   2
    ##    A3 ~1      4.57 0.033 < 0.0001   3   1   3
    ##    A4 ~1      4.68 0.038 < 0.0001   4   1   4
    ##    A5 ~1      4.54 0.033 < 0.0001   5   1   5
    ##    C1 ~1      4.56 0.031 < 0.0001   6   1   6
    ##    C2 ~1      4.39 0.034 < 0.0001   7   1   7
    ##    C3 ~1      4.28 0.033 < 0.0001   8   1   8
    ##    C4 ~1      2.56 0.035 < 0.0001   9   1   9
    ##    C5 ~1      3.33 0.041 < 0.0001  10   1  10
    ##    E1 ~1      3.04 0.041 < 0.0001  11   1  11
    ##    E2 ~1      3.21 0.041 < 0.0001  12   1  12
    ##    E3 ~1      3.95 0.035 < 0.0001  13   1  13
    ##    E4 ~1      4.35 0.038 < 0.0001  14   1  14
    ##    E5 ~1      4.37 0.034 < 0.0001  15   1  15
    ##    N1 ~1      2.89 0.040 < 0.0001  16   1  16
    ##    N2 ~1      3.44 0.039 < 0.0001  17   1  17
    ##    N3 ~1      3.21 0.040 < 0.0001  18   1  18
    ##    N4 ~1      3.21 0.040 < 0.0001  19   1  19
    ##    N5 ~1      2.98 0.041 < 0.0001  20   1  20
    ##    O1 ~1      4.85 0.028 < 0.0001  21   1  21
    ##    O2 ~1      2.64 0.038 < 0.0001  22   1  22
    ##    O3 ~1      4.45 0.030 < 0.0001  23   1  23
    ##    O4 ~1      4.93 0.030 < 0.0001  24   1  24
    ##    O5 ~1      2.46 0.033 < 0.0001  25   1  25
    ## 
    ##  -  omega (symmetric) 
    ##  var1 op var2      est    se        p row col par
    ##    A2 --   A1    -0.23 0.024 < 0.0001   2   1  26
    ##    A3 --   A1    -0.16 0.025 < 0.0001   3   1  27
    ##    A4 --   A1   -0.017 0.025     0.51   4   1  28
    ##    A5 --   A1   -0.019 0.025     0.46   5   1  29
    ##    C1 --   A1    0.044 0.025    0.081   6   1  30
    ##    C2 --   A1    0.052 0.025    0.038   7   1  31
    ##    C3 --   A1    0.023 0.025     0.37   8   1  32
    ##    C4 --   A1    0.097 0.025  0.00011   9   1  33
    ##    C5 --   A1   -0.021 0.025     0.40  10   1  34
    ##    E1 --   A1    0.067 0.025   0.0074  11   1  35
    ##    E2 --   A1    0.032 0.025     0.21  12   1  36
    ##    E3 --   A1    0.042 0.025    0.093  13   1  37
    ##    E4 --   A1    0.091 0.025  0.00028  14   1  38
    ##    E5 --   A1    0.020 0.025     0.43  15   1  39
    ##    N1 --   A1    0.071 0.025   0.0047  16   1  40
    ##    N2 --   A1    0.032 0.025     0.20  17   1  41
    ##    N3 --   A1    0.042 0.025    0.098  18   1  42
    ##    N4 --   A1   -0.057 0.025    0.023  19   1  43
    ##    N5 --   A1   -0.042 0.025    0.098  20   1  44
    ##    O1 --   A1    0.050 0.025    0.046  21   1  45
    ##    O2 --   A1    0.037 0.025     0.15  22   1  46
    ##    O3 --   A1  -0.0092 0.025     0.72  23   1  47
    ##    O4 --   A1   -0.088 0.025  0.00049  24   1  48
    ##    O5 --   A1    0.017 0.025     0.50  25   1  49
    ##    A3 --   A2     0.26 0.024 < 0.0001   3   2  50
    ##    A4 --   A2     0.19 0.024 < 0.0001   4   2  51
    ##    A5 --   A2     0.13 0.025 < 0.0001   5   2  52
    ##    C1 --   A2   -0.040 0.025     0.12   6   2  53
    ##    C2 --   A2   0.0078 0.025     0.76   7   2  54
    ##    C3 --   A2    0.095 0.025  0.00016   8   2  55
    ##    C4 --   A2   -0.029 0.025     0.25   9   2  56
    ##    C5 --   A2    0.052 0.025    0.041  10   2  57
    ##    E1 --   A2   -0.044 0.025    0.083  11   2  58
    ##    E2 --   A2   -0.035 0.025     0.17  12   2  59
    ##    E3 --   A2   -0.067 0.025   0.0075  13   2  60
    ##    E4 --   A2   0.0099 0.025     0.70  14   2  61
    ##    E5 --   A2     0.19 0.024 < 0.0001  15   2  62
    ##    N1 --   A2   -0.072 0.025   0.0044  16   2  63
    ##    N2 --   A2    0.067 0.025   0.0076  17   2  64
    ##    N3 --   A2   -0.017 0.025     0.50  18   2  65
    ##    N4 --   A2    0.052 0.025    0.040  19   2  66
    ##    N5 --   A2    0.038 0.025     0.13  20   2  67
    ##    O1 --   A2    0.011 0.025     0.66  21   2  68
    ##    O2 --   A2    0.042 0.025    0.095  22   2  69
    ##    O3 --   A2   -0.018 0.025     0.48  23   2  70
    ##    O4 --   A2    0.065 0.025    0.010  24   2  71
    ##    O5 --   A2   -0.033 0.025     0.20  25   2  72
    ##    A4 --   A3     0.15 0.025 < 0.0001   4   3  73
    ##    A5 --   A3     0.27 0.023 < 0.0001   5   3  74
    ##    C1 --   A3  -0.0012 0.025     0.96   6   3  75
    ##    C2 --   A3    0.021 0.025     0.41   7   3  76
    ##    C3 --   A3  -0.0057 0.025     0.82   8   3  77
    ##    C4 --   A3    0.023 0.025     0.37   9   3  78
    ##    C5 --   A3   -0.010 0.025     0.69  10   3  79
    ##    E1 --   A3    0.031 0.025     0.22  11   3  80
    ##    E2 --   A3   -0.017 0.025     0.50  12   3  81
    ##    E3 --   A3     0.16 0.025 < 0.0001  13   3  82
    ##    E4 --   A3    0.090 0.025  0.00033  14   3  83
    ##    E5 --   A3   -0.037 0.025     0.14  15   3  84
    ##    N1 --   A3    0.029 0.025     0.25  16   3  85
    ##    N2 --   A3   0.0083 0.025     0.74  17   3  86
    ##    N3 --   A3    0.012 0.025     0.63  18   3  87
    ##    N4 --   A3  -0.0030 0.025     0.91  19   3  88
    ##    N5 --   A3   -0.030 0.025     0.24  20   3  89
    ##    O1 --   A3    0.026 0.025     0.30  21   3  90
    ##    O2 --   A3    0.021 0.025     0.41  22   3  91
    ##    O3 --   A3    0.071 0.025   0.0049  23   3  92
    ##    O4 --   A3   -0.030 0.025     0.24  24   3  93
    ##    O5 --   A3    0.016 0.025     0.53  25   3  94
    ##    A5 --   A4    0.044 0.025    0.079   5   4  95
    ##    C1 --   A4   -0.043 0.025    0.086   6   4  96
    ##    C2 --   A4     0.13 0.025 < 0.0001   7   4  97
    ##    C3 --   A4   -0.026 0.025     0.31   8   4  98
    ##    C4 --   A4   -0.037 0.025     0.14   9   4  99
    ##    C5 --   A4    -0.12 0.025 < 0.0001  10   4 100
    ##    E1 --   A4    0.021 0.025     0.40  11   4 101
    ##    E2 --   A4    0.021 0.025     0.41  12   4 102
    ##    E3 --   A4   0.0071 0.025     0.78  13   4 103
    ##    E4 --   A4     0.12 0.025 < 0.0001  14   4 104
    ##    E5 --   A4   -0.049 0.025    0.053  15   4 105
    ##    N1 --   A4    0.045 0.025    0.075  16   4 106
    ##    N2 --   A4    -0.10 0.025 < 0.0001  17   4 107
    ##    N3 --   A4    0.042 0.025     0.10  18   4 108
    ##    N4 --   A4   -0.035 0.025     0.17  19   4 109
    ##    N5 --   A4    0.031 0.025     0.22  20   4 110
    ##    O1 --   A4  -0.0083 0.025     0.74  21   4 111
    ##    O2 --   A4  -0.0081 0.025     0.75  22   4 112
    ##    O3 --   A4   -0.018 0.025     0.48  23   4 113
    ##    O4 --   A4   -0.024 0.025     0.34  24   4 114
    ##    O5 --   A4    0.032 0.025     0.21  25   4 115
    ##    C1 --   A5    0.028 0.025     0.27   6   5 116
    ##    C2 --   A5   -0.066 0.025   0.0091   7   5 117
    ##    C3 --   A5    0.041 0.025     0.11   8   5 118
    ##    C4 --   A5  -0.0014 0.025     0.96   9   5 119
    ##    C5 --   A5    0.025 0.025     0.32  10   5 120
    ##    E1 --   A5  -0.0087 0.025     0.73  11   5 121
    ##    E2 --   A5   -0.014 0.025     0.57  12   5 122
    ##    E3 --   A5     0.16 0.025 < 0.0001  13   5 123
    ##    E4 --   A5     0.20 0.024 < 0.0001  14   5 124
    ##    E5 --   A5   0.0042 0.025     0.87  15   5 125
    ##    N1 --   A5   -0.045 0.025    0.078  16   5 126
    ##    N2 --   A5   -0.075 0.025   0.0028  17   5 127
    ##    N3 --   A5   0.0018 0.025     0.94  18   5 128
    ##    N4 --   A5   -0.015 0.025     0.56  19   5 129
    ##    N5 --   A5    0.055 0.025    0.029  20   5 130
    ##    O1 --   A5 -0.00093 0.025     0.97  21   5 131
    ##    O2 --   A5    0.015 0.025     0.54  22   5 132
    ##    O3 --   A5    0.035 0.025     0.16  23   5 133
    ##    O4 --   A5    0.019 0.025     0.46  24   5 134
    ##    O5 --   A5   -0.014 0.025     0.57  25   5 135
    ##    C2 --   C1     0.28 0.023 < 0.0001   7   6 136
    ##    C3 --   C1     0.12 0.025 < 0.0001   8   6 137
    ##    C4 --   C1    -0.15 0.025 < 0.0001   9   6 138
    ##    C5 --   C1   -0.037 0.025     0.15  10   6 139
    ##    E1 --   C1    0.054 0.025    0.034  11   6 140
    ##    E2 --   C1    0.054 0.025    0.033  12   6 141
    ##    E3 --   C1  -0.0014 0.025     0.95  13   6 142
    ##    E4 --   C1    0.081 0.025   0.0012  14   6 143
    ##    E5 --   C1     0.13 0.025 < 0.0001  15   6 144
    ##    N1 --   C1   -0.039 0.025     0.13  16   6 145
    ##    N2 --   C1   -0.015 0.025     0.56  17   6 146
    ##    N3 --   C1    0.056 0.025    0.027  18   6 147
    ##    N4 --   C1   -0.030 0.025     0.23  19   6 148
    ##    N5 --   C1    0.018 0.025     0.47  20   6 149
    ##    O1 --   C1    0.029 0.025     0.25  21   6 150
    ##    O2 --   C1   -0.056 0.025    0.025  22   6 151
    ##    O3 --   C1    0.034 0.025     0.17  23   6 152
    ##    O4 --   C1    0.070 0.025   0.0053  24   6 153
    ##    O5 --   C1   -0.035 0.025     0.17  25   6 154
    ##    C3 --   C2     0.13 0.025 < 0.0001   8   7 155
    ##    C4 --   C2    -0.22 0.024 < 0.0001   9   7 156
    ##    C5 --   C2   -0.082 0.025   0.0012  10   7 157
    ##    E1 --   C2    0.087 0.025  0.00057  11   7 158
    ##    E2 --   C2    0.025 0.025     0.32  12   7 159
    ##    E3 --   C2    0.036 0.025     0.15  13   7 160
    ##    E4 --   C2    0.013 0.025     0.61  14   7 161
    ##    E5 --   C2    0.060 0.025    0.017  15   7 162
    ##    N1 --   C2   0.0049 0.025     0.85  16   7 163
    ##    N2 --   C2  -0.0047 0.025     0.85  17   7 164
    ##    N3 --   C2    0.027 0.025     0.28  18   7 165
    ##    N4 --   C2    0.058 0.025    0.021  19   7 166
    ##    N5 --   C2    0.092 0.025  0.00025  20   7 167
    ##    O1 --   C2    0.033 0.025     0.20  21   7 168
    ##    O2 --   C2    0.065 0.025    0.010  22   7 169
    ##    O3 --   C2     0.12 0.025 < 0.0001  23   7 170
    ##    O4 --   C2   0.0099 0.025     0.70  24   7 171
    ##    O5 --   C2    0.022 0.025     0.39  25   7 172
    ##    C4 --   C3    -0.17 0.025 < 0.0001   9   8 173
    ##    C5 --   C3    -0.18 0.025 < 0.0001  10   8 174
    ##    E1 --   C3    0.025 0.025     0.33  11   8 175
    ##    E2 --   C3    0.042 0.025     0.10  12   8 176
    ##    E3 --   C3   -0.023 0.025     0.36  13   8 177
    ##    E4 --   C3   -0.011 0.025     0.65  14   8 178
    ##    E5 --   C3    0.050 0.025    0.050  15   8 179
    ##    N1 --   C3   -0.011 0.025     0.67  16   8 180
    ##    N2 --   C3    0.021 0.025     0.41  17   8 181
    ##    N3 --   C3  -0.0031 0.025     0.90  18   8 182
    ##    N4 --   C3  -0.0062 0.025     0.81  19   8 183
    ##    N5 --   C3   0.0049 0.025     0.85  20   8 184
    ##    O1 --   C3    0.024 0.025     0.35  21   8 185
    ##    O2 --   C3    0.040 0.025     0.11  22   8 186
    ##    O3 --   C3   -0.033 0.025     0.19  23   8 187
    ##    O4 --   C3    0.047 0.025    0.065  24   8 188
    ##    O5 --   C3    0.071 0.025   0.0047  25   8 189
    ##    C5 --   C4     0.28 0.023 < 0.0001  10   9 190
    ##    E1 --   C4    0.058 0.025    0.020  11   9 191
    ##    E2 --   C4    0.066 0.025   0.0092  12   9 192
    ##    E3 --   C4    0.030 0.025     0.24  13   9 193
    ##    E4 --   C4    0.061 0.025    0.016  14   9 194
    ##    E5 --   C4  -0.0048 0.025     0.85  15   9 195
    ##    N1 --   C4    0.046 0.025    0.069  16   9 196
    ##    N2 --   C4   -0.057 0.025    0.023  17   9 197
    ##    N3 --   C4    0.046 0.025    0.069  18   9 198
    ##    N4 --   C4    0.055 0.025    0.029  19   9 199
    ##    N5 --   C4    0.065 0.025   0.0099  20   9 200
    ##    O1 --   C4    0.045 0.025    0.072  21   9 201
    ##    O2 --   C4    0.098 0.025 < 0.0001  22   9 202
    ##    O3 --   C4    0.089 0.025  0.00038  23   9 203
    ##    O4 --   C4    0.022 0.025     0.39  24   9 204
    ##    O5 --   C4     0.12 0.025 < 0.0001  25   9 205
    ##    E1 --   C5   -0.076 0.025   0.0026  11  10 206
    ##    E2 --   C5    0.086 0.025  0.00066  12  10 207
    ##    E3 --   C5   -0.033 0.025     0.20  13  10 208
    ##    E4 --   C5   -0.039 0.025     0.12  14  10 209
    ##    E5 --   C5   -0.054 0.025    0.032  15  10 210
    ##    N1 --   C5   -0.039 0.025     0.12  16  10 211
    ##    N2 --   C5     0.10 0.025 < 0.0001  17  10 212
    ##    N3 --   C5    0.019 0.025     0.45  18  10 213
    ##    N4 --   C5     0.17 0.025 < 0.0001  19  10 214
    ##    N5 --   C5   -0.025 0.025     0.32  20  10 215
    ##    O1 --   C5   -0.020 0.025     0.43  21  10 216
    ##    O2 --   C5    0.099 0.025 < 0.0001  22  10 217
    ##    O3 --   C5    0.034 0.025     0.18  23  10 218
    ##    O4 --   C5     0.10 0.025 < 0.0001  24  10 219
    ##    O5 --   C5   -0.027 0.025     0.28  25  10 220
    ##    E2 --   E1     0.22 0.024 < 0.0001  12  11 221
    ##    E3 --   E1   -0.074 0.025   0.0032  13  11 222
    ##    E4 --   E1    -0.18 0.024 < 0.0001  14  11 223
    ##    E5 --   E1    -0.10 0.025 < 0.0001  15  11 224
    ##    N1 --   E1   -0.031 0.025     0.22  16  11 225
    ##    N2 --   E1   -0.012 0.025     0.63  17  11 226
    ##    N3 --   E1   -0.022 0.025     0.38  18  11 227
    ##    N4 --   E1    0.095 0.025  0.00015  19  11 228
    ##    N5 --   E1   -0.084 0.025  0.00085  20  11 229
    ##    O1 --   E1   -0.026 0.025     0.31  21  11 230
    ##    O2 --   E1   -0.013 0.025     0.60  22  11 231
    ##    O3 --   E1   -0.075 0.025   0.0028  23  11 232
    ##    O4 --   E1    0.061 0.025    0.015  24  11 233
    ##    O5 --   E1    0.053 0.025    0.036  25  11 234
    ##    E3 --   E2    -0.10 0.025 < 0.0001  13  12 235
    ##    E4 --   E2    -0.31 0.023 < 0.0001  14  12 236
    ##    E5 --   E2    -0.12 0.025 < 0.0001  15  12 237
    ##    N1 --   E2   -0.031 0.025     0.22  16  12 238
    ##    N2 --   E2    0.050 0.025    0.048  17  12 239
    ##    N3 --   E2   -0.012 0.025     0.62  18  12 240
    ##    N4 --   E2    0.067 0.025   0.0080  19  12 241
    ##    N5 --   E2     0.14 0.025 < 0.0001  20  12 242
    ##    O1 --   E2  -0.0083 0.025     0.74  21  12 243
    ##    O2 --   E2    0.016 0.025     0.53  22  12 244
    ##    O3 --   E2   -0.021 0.025     0.41  23  12 245
    ##    O4 --   E2     0.13 0.025 < 0.0001  24  12 246
    ##    O5 --   E2    0.041 0.025     0.10  25  12 247
    ##    E4 --   E3     0.11 0.025 < 0.0001  14  13 248
    ##    E5 --   E3     0.14 0.025 < 0.0001  15  13 249
    ##    N1 --   E3   -0.013 0.025     0.60  16  13 250
    ##    N2 --   E3  -0.0011 0.025     0.97  17  13 251
    ##    N3 --   E3    0.085 0.025  0.00069  18  13 252
    ##    N4 --   E3   -0.012 0.025     0.64  19  13 253
    ##    N5 --   E3   -0.018 0.025     0.48  20  13 254
    ##    O1 --   E3     0.17 0.025 < 0.0001  21  13 255
    ##    O2 --   E3   0.0056 0.025     0.82  22  13 256
    ##    O3 --   E3     0.16 0.025 < 0.0001  23  13 257
    ##    O4 --   E3    0.033 0.025     0.20  24  13 258
    ##    O5 --   E3   0.0039 0.025     0.88  25  13 259
    ##    E5 --   E4    0.033 0.025     0.19  15  14 260
    ##    N1 --   E4   -0.032 0.025     0.20  16  14 261
    ##    N2 --   E4  -0.0052 0.025     0.84  17  14 262
    ##    N3 --   E4    0.017 0.025     0.51  18  14 263
    ##    N4 --   E4   -0.097 0.025  0.00010  19  14 264
    ##    N5 --   E4    0.043 0.025    0.086  20  14 265
    ##    O1 --   E4   -0.031 0.025     0.22  21  14 266
    ##    O2 --   E4     0.10 0.025 < 0.0001  22  14 267
    ##    O3 --   E4    0.064 0.025    0.011  23  14 268
    ##    O4 --   E4   0.0087 0.025     0.73  24  14 269
    ##    O5 --   E4     0.13 0.025 < 0.0001  25  14 270
    ##    N1 --   E5    0.098 0.025 < 0.0001  16  15 271
    ##    N2 --   E5    0.091 0.025  0.00027  17  15 272
    ##    N3 --   E5   -0.057 0.025    0.024  18  15 273
    ##    N4 --   E5   -0.063 0.025    0.013  19  15 274
    ##    N5 --   E5   -0.091 0.025  0.00029  20  15 275
    ##    O1 --   E5     0.12 0.025 < 0.0001  21  15 276
    ##    O2 --   E5    0.023 0.025     0.36  22  15 277
    ##    O3 --   E5    0.088 0.025  0.00046  23  15 278
    ##    O4 --   E5  -0.0073 0.025     0.77  24  15 279
    ##    O5 --   E5   0.0044 0.025     0.86  25  15 280
    ##    N2 --   N1     0.55 0.018 < 0.0001  17  16 281
    ##    N3 --   N1     0.23 0.024 < 0.0001  18  16 282
    ##    N4 --   N1     0.11 0.025 < 0.0001  19  16 283
    ##    N5 --   N1     0.11 0.025 < 0.0001  20  16 284
    ##    O1 --   N1    0.021 0.025     0.41  21  16 285
    ##    O2 --   N1    0.012 0.025     0.62  22  16 286
    ##    O3 --   N1    0.011 0.025     0.68  23  16 287
    ##    O4 --   N1   -0.040 0.025     0.11  24  16 288
    ##    O5 --   N1    0.056 0.025    0.028  25  16 289
    ##    N3 --   N2     0.17 0.025 < 0.0001  18  17 290
    ##    N4 --   N2    0.030 0.025     0.23  19  17 291
    ##    N5 --   N2    0.045 0.025    0.078  20  17 292
    ##    O1 --   N2   -0.021 0.025     0.41  21  17 293
    ##    O2 --   N2    0.034 0.025     0.18  22  17 294
    ##    O3 --   N2    0.017 0.025     0.51  23  17 295
    ##    O4 --   N2    0.034 0.025     0.18  24  17 296
    ##    O5 --   N2   -0.042 0.025    0.094  25  17 297
    ##    N4 --   N3     0.27 0.023 < 0.0001  19  18 298
    ##    N5 --   N3     0.17 0.025 < 0.0001  20  18 299
    ##    O1 --   N3  -0.0021 0.025     0.93  21  18 300
    ##    O2 --   N3   0.0029 0.025     0.91  22  18 301
    ##    O3 --   N3   -0.039 0.025     0.12  23  18 302
    ##    O4 --   N3    0.074 0.025   0.0032  24  18 303
    ##    O5 --   N3  -0.0034 0.025     0.89  25  18 304
    ##    N5 --   N4     0.14 0.025 < 0.0001  20  19 305
    ##    O1 --   N4    0.043 0.025    0.088  21  19 306
    ##    O2 --   N4   -0.047 0.025    0.063  22  19 307
    ##    O3 --   N4    0.042 0.025     0.10  23  19 308
    ##    O4 --   N4    0.059 0.025    0.019  24  19 309
    ##    O5 --   N4    0.016 0.025     0.54  25  19 310
    ##    O1 --   N5   -0.055 0.025    0.029  21  20 311
    ##    O2 --   N5    0.096 0.025  0.00013  22  20 312
    ##    O3 --   N5  -0.0096 0.025     0.71  23  20 313
    ##    O4 --   N5    0.012 0.025     0.63  24  20 314
    ##    O5 --   N5    0.045 0.025    0.076  25  20 315
    ##    O2 --   O1    -0.13 0.025 < 0.0001  22  21 316
    ##    O3 --   O1     0.16 0.025 < 0.0001  23  21 317
    ##    O4 --   O1     0.12 0.025 < 0.0001  24  21 318
    ##    O5 --   O1   -0.087 0.025  0.00050  25  21 319
    ##    O3 --   O2    -0.18 0.024 < 0.0001  23  22 320
    ##    O4 --   O2   -0.011 0.025     0.67  24  22 321
    ##    O5 --   O2     0.20 0.024 < 0.0001  25  22 322
    ##    O4 --   O3    0.084 0.025  0.00087  24  23 323
    ##    O5 --   O3    -0.17 0.025 < 0.0001  25  23 324
    ##    O5 --   O4    -0.12 0.025 < 0.0001  25  24 325
    ## 
    ##  -  delta (diagonal) 
    ##  var1  op var2  est    se        p row col par
    ##    A1 ~/~   A1 1.27 0.023 < 0.0001   1   1 326
    ##    A2 ~/~   A2 0.93 0.017 < 0.0001   2   2 327
    ##    A3 ~/~   A3 0.96 0.017 < 0.0001   3   3 328
    ##    A4 ~/~   A4 1.26 0.023 < 0.0001   4   4 329
    ##    A5 ~/~   A5 0.96 0.017 < 0.0001   5   5 330
    ##    C1 ~/~   C1 1.02 0.018 < 0.0001   6   6 331
    ##    C2 ~/~   C2 1.07 0.019 < 0.0001   7   7 332
    ##    C3 ~/~   C3 1.14 0.020 < 0.0001   8   8 333
    ##    C4 ~/~   C4 1.06 0.019 < 0.0001   9   9 334
    ##    C5 ~/~   C5 1.26 0.023 < 0.0001  10  10 335
    ##    E1 ~/~   E1 1.34 0.024 < 0.0001  11  11 336
    ##    E2 ~/~   E2 1.18 0.021 < 0.0001  12  12 337
    ##    E3 ~/~   E3 1.04 0.019 < 0.0001  13  13 338
    ##    E4 ~/~   E4 1.06 0.019 < 0.0001  14  14 339
    ##    E5 ~/~   E5 1.09 0.019 < 0.0001  15  15 340
    ##    N1 ~/~   N1  1.0 0.018 < 0.0001  16  16 341
    ##    N2 ~/~   N2 1.00 0.018 < 0.0001  17  17 342
    ##    N3 ~/~   N3 1.15 0.021 < 0.0001  18  18 343
    ##    N4 ~/~   N4 1.17 0.021 < 0.0001  19  19 344
    ##    N5 ~/~   N5 1.36 0.024 < 0.0001  20  20 345
    ##    O1 ~/~   O1 0.94 0.017 < 0.0001  21  21 346
    ##    O2 ~/~   O2 1.31 0.023 < 0.0001  22  22 347
    ##    O3 ~/~   O3 0.97 0.017 < 0.0001  23  23 348
    ##    O4 ~/~   O4 1.08 0.019 < 0.0001  24  24 349
    ##    O5 ~/~   O5 1.16 0.021 < 0.0001  25  25 350

The edge *O1 -- A5* is estimated to be near zero, so let's try removing this edge:

``` r
# fix the parameter to zero:
model2 <- model %>% fixpar("omega", "O1", "A5", value=0)
```

    ## Fixed 1 parameters!

``` r
# Evaluate the new model:
model2 <- model2 %>% runmodel
```

    ## Estimating model...

    ## Computing modification indices...

We can compare the two models:

``` r
compare(model, model2)
```

    ##    model DF       AIC       BIC    Chisq Chisq_diff DF_diff p_value
    ##  Model 1  0 125975.19 127848.99 < 0.0001                           
    ##  Model 2  1 125973.20 127841.64   0.0033     0.0033       1    0.95
    ## 
    ## Note: Chi-square difference test assumes models are nested.

Removing the edge improved our model. The *prune* function can be used to automatically and recursively remove any parameter that is not significant at some *α* level

``` r
prunedmodel <- model2 %>% prune(alpha = 0.01, recursive = TRUE)
```

    ## Clearing 199 parameters!

    ## Estimating model...

    ## Computing modification indices...

    ## Clearing 6 parameters!

    ## Estimating model...

    ## Computing modification indices...

``` r
compare(saturated = model, nearsaturated = model2, pruned = prunedmodel)
```

    ##          model  DF       AIC       BIC    Chisq Chisq_diff DF_diff
    ##      saturated   0 125975.19 127848.99 < 0.0001                   
    ##  nearsaturated   1 125973.20 127841.64   0.0033     0.0033       1
    ##         pruned 206 126082.77 126853.71   519.58     519.57     205
    ##   p_value
    ##          
    ##      0.95
    ##  < 0.0001
    ## 
    ## Note: Chi-square difference test assumes models are nested.

The pruned model has a much better BIC, although the comparison with regard to AIC and *χ*<sup>2</sup> is less strong in favor of this pruned model. We can also look at the modification indices:

``` r
prunedmodel %>% MIs
```

    ## 
    ## Top 10 modification indices:
    ## 
    ##  var1 op var2 est    mi      pmi epc matrix row col group
    ##    C1 --   A1   0 19.53 < 0.0001      omega   6   1     1
    ##    C2 --   A1   0 14.80  0.00012      omega   7   1     1
    ##    E5 --   C2   0 14.54  0.00014      omega  15   7     1
    ##    N5 --   A4   0 13.62  0.00022      omega  20   4     1
    ##    O5 --   N5   0 11.79  0.00060      omega  25  20     1
    ##    O4 --   N4   0 11.67  0.00064      omega  24  19     1
    ##    E5 --   C5   0 11.56  0.00067      omega  15  10     1
    ##    C5 --   A1   0 11.25  0.00079      omega  10   1     1
    ##    O5 --   N1   0 10.82   0.0010      omega  25  16     1
    ##    O5 --   E2   0 10.73   0.0011      omega  25  12     1

By default, only the top 10 modification indices are shown, but this can be changed with `MIs(all = TRUE)`. We can try to add one edge to the model:

``` r
# Add an edge and evaluate:
model3 <- prunedmodel %>% freepar("omega","A1","C1") %>% runmodel
```

    ## No parameters need to be freed

    ## Estimating model...

    ## Computing modification indices...

``` r
# compare:
compare(saturated = model, nearsaturated = model2,
        pruned = prunedmodel, lastmodel = model3)
```

    ##          model  DF       AIC       BIC    Chisq Chisq_diff  DF_diff
    ##      saturated   0 125975.19 127848.99 < 0.0001                    
    ##  nearsaturated   1 125973.20 127841.64   0.0033     0.0033        1
    ##         pruned 206 126082.77 126853.71   519.58     519.57      205
    ##      lastmodel 206 126082.69 126853.63   519.50      0.076 < 0.0001
    ##   p_value
    ##          
    ##      0.95
    ##  < 0.0001
    ##  < 0.0001
    ## 
    ## Note: Chi-square difference test assumes models are nested.

We can see that this model fits better. There is an automated function in *psychonetrics* that adds edges according to modification indices at a given level of *α*, which will stop searching by default if BIC is no longer increased:

``` r
# Stepup search:
model_stepup <- model3 %>% stepup

# compare:
compare(saturated = model, nearsaturated = model2, pruned = prunedmodel, oneMIadded = model3, stepup = model_stepup)
```

    ##          model  DF       AIC       BIC    Chisq Chisq_diff  DF_diff
    ##      saturated   0 125975.19 127848.99 < 0.0001                    
    ##  nearsaturated   1 125973.20 127841.64   0.0033     0.0033        1
    ##         stepup 192 125954.23 126800.11   363.03     363.03      191
    ##         pruned 206 126082.77 126853.71   519.58     156.54       14
    ##     oneMIadded 206 126082.69 126853.63   519.50      0.076 < 0.0001
    ##   p_value
    ##          
    ##      0.95
    ##  < 0.0001
    ##  < 0.0001
    ##  < 0.0001
    ## 
    ## Note: Chi-square difference test assumes models are nested.

Which leads to a best fitting model according to all indices.. This leads to an efficient estimation algorithm by combining the two. For example, `model %>% runmodel %>% prune %>% stepup` will estimate a model structure relatively fast, while `model %>% runmodel %>% stepup %>% prune` will be slower but generally more conservative. all these model modifications we made are actually stored in the *psychonetrics* object. For example:

``` r
data.frame(
  timestamp = sapply(model_stepup@log,function(x)as.character(x@time)),
  event = sapply(model_stepup@log,function(x)as.character(x@event))
)
```

             timestamp

1 2019-03-05 14:25:49 2 2019-03-05 21:46:17 3 2019-03-05 21:46:17 4 2019-03-05 21:47:22 5 2019-03-05 21:47:22 6 2019-03-05 21:47:36 7 2019-03-05 21:47:36 8 2019-03-05 21:47:46 9 2019-03-05 21:47:46 10 2019-03-05 21:47:54 11 2019-03-05 21:50:43 event 1 Model created 2 Evaluated model 3 Fixed element(s) of omega: 1 parameters! 4 Evaluated model 5 Pruned all parameters in matrices omega at alpha = 0.01 6 Evaluated model 7 Pruned all parameters in matrices omega at alpha = 0.01 8 Evaluated model 9 Pruned all parameters in matrices omega at alpha = 0.01 (none were removed) 10 Evaluated model 11 Performed step-up model search

Obtaining a network from bootnet
--------------------------------

The *psychonetrics* is not intended, however, to replace exploratory model search functions. Instead, the bootnet package provides more algorithms that may be faster. For example, we could also use *qgraph*'s `ggmModSelect` function:

``` r
library("bootnet")
net <- estimateNetwork(trainingData[,vars], default = "ggmModSelect", verbose = FALSE)
```

We can transform this into a psychonetrics object and refit the model:

``` r
model_frombootnet <- frombootnet(net) %>% runmodel
```

    ## Estimating baseline model...

    ## Estimating model...

    ## Computing modification indices...

And compare it to our model:

``` r
compare(ggmModSelect = model_frombootnet, psychonetrics = model_stepup)
```

    ##          model  DF       AIC       BIC  Chisq Chisq_diff DF_diff  p_value
    ##   ggmModSelect 173 125905.33 126852.94 276.14                            
    ##  psychonetrics 192 125954.23 126800.11 363.03      86.90      19 < 0.0001
    ## 
    ## Note: Chi-square difference test assumes models are nested.

Which results in a comparable fit.

Plotting the network
--------------------

The package does not contain plotting methods yet, but does return the estimated network structure, which we can plot:

``` r
library("qgraph")
stepup_net <- model_stepup@modelmatrices$`1`$omega
bootnet_net <- model_frombootnet@modelmatrices$`1`$omega
L <- averageLayout(as.matrix(stepup_net), as.matrix(bootnet_net))
layout(t(1:2))
qgraph(stepup_net, labels = vars, theme = "colorblind", 
       title = "Psychonetrics estimation", layout = L)
qgraph(bootnet_net, labels = vars, theme = "colorblind", 
       title = "ggmModSelect estimation", layout = L)
```

![](readme_files/figure-markdown_github/unnamed-chunk-19-1.png)

Confirmatory fit
----------------

Now let's take our model (let's use only the psychonetrics estimated model) and fit it to the test data:

``` r
adjacency <- 1*(model_stepup@modelmatrices$`1`$omega!=0)
confirmatory <- ggm(testData, vars = vars, omega = adjacency)
confirmatory <- confirmatory %>% runmodel
```

    ## Estimating baseline model...

    ## Estimating model...

    ## Computing modification indices...

The model shows very good fit to the test data:

``` r
confirmatory %>% fit
```

    ##            Measure     Value
    ##               logl -35112.68
    ##  unrestricted.logl -34874.93
    ##      baseline.logl -38155.76
    ##               nvar        25
    ##               nobs       350
    ##               npar       158
    ##                 df       192
    ##          objective     34.40
    ##              chisq    475.49
    ##             pvalue        ~0
    ##     baseline.chisq   6561.66
    ##        baseline.df       300
    ##    baseline.pvalue        ~0
    ##                nfi      0.93
    ##               pnfi      0.59
    ##                tli      0.93
    ##               nnfi      0.93
    ##                rfi      0.89
    ##                ifi      0.96
    ##                rni      0.95
    ##                cfi      0.95
    ##              rmsea     0.041
    ##     rmsea.ci.lower     0.036
    ##     rmsea.ci.upper     0.046
    ##       rmsea.pvalue       1.0
    ##             aic.ll  70541.36
    ##            aic.ll2  70611.63
    ##              aic.x     91.49
    ##             aic.x2    791.49
    ##                bic  71295.50
    ##               bic2  70793.73
    ##               ebic  71804.09
    ##         ebicTuning      0.25

Multi-group analysis
--------------------

We can create a multi-group model using the `groups` argument (see also Kan, Maas, and Levine (2019)). Let's fit our model to both groups separately:

``` r
groupmodel <- ggm(trainingData, vars = vars, omega = adjacency, groups = "gender") %>% runmodel
```

This model fits very well:

``` r
groupmodel %>% fit
```

    ##            Measure     Value
    ##               logl -62572.95
    ##  unrestricted.logl -62243.20
    ##      baseline.logl -68447.73
    ##               nvar        25
    ##               nobs       700
    ##               npar       316
    ##                 df       384
    ##          objective     34.17
    ##              chisq    659.50
    ##             pvalue        ~0
    ##     baseline.chisq  12409.06
    ##        baseline.df       600
    ##    baseline.pvalue        ~0
    ##                nfi      0.95
    ##               pnfi      0.61
    ##                tli      0.96
    ##               nnfi      0.96
    ##                rfi      0.92
    ##                ifi      0.98
    ##                rni      0.98
    ##                cfi      0.98
    ##              rmsea     0.030
    ##     rmsea.ci.lower     0.026
    ##     rmsea.ci.upper     0.034
    ##       rmsea.pvalue         1
    ##             aic.ll 125777.90
    ##            aic.ll2 125938.82
    ##              aic.x   -108.50
    ##             aic.x2   1291.50
    ##                bic 127469.68
    ##               bic2 126465.82
    ##               ebic 128486.84
    ##         ebicTuning      0.25

Next, I can constrain all edges to be equal:

``` r
# Run model:
groupmodel_2 <- groupmodel %>% groupequal("omega") %>% runmodel
```

    ## Constrained 108 parameters!

    ## Estimating model...

    ## Computing modification indices...

``` r
# Compare models:
compare(configural = groupmodel, metric = groupmodel_2)
```

    ##       model  DF       AIC       BIC  Chisq Chisq_diff DF_diff  p_value
    ##  configural 384 125777.90 127469.68 659.50                            
    ##      metric 492 125737.42 126851.00 835.02     175.52     108 < 0.0001
    ## 
    ## Note: Chi-square difference test assumes models are nested.

The model with constrained edges is confirmed. However, I can look at the equality-free modification indices (note: this is very experimental still and may not be the best method for this):

``` r
groupmodel_2 %>% MIs(type = "free")
```

    ## 
    ## Top 10 equality-free modification indices:
    ## 
    ##  var1 op var2        est mi_free pmi_free epc_free matrix row col group
    ##    E5 --   E4 0.00000000   25.04 < 0.0001           omega  15  14     2
    ##    E5 --   A4 0.00000000   12.93  0.00032           omega  15   4     1
    ##    E2 --   C5 0.08490288   12.74  0.00036           omega  12  10     2
    ##    C4 --   A3 0.00000000   12.36  0.00044           omega   9   3     2
    ##    N4 --   A3 0.00000000   10.86  0.00098           omega  19   3     2
    ##    E4 --   E3 0.12535518   10.81   0.0010           omega  14  13     2
    ##    C5 --   A3 0.00000000   10.63   0.0011           omega  10   3     2
    ##    E3 --   E1 0.00000000   10.23   0.0014           omega  13  11     2
    ##    E1 --   C5 0.00000000    9.98   0.0016           omega  11  10     1
    ##    O3 --   A4 0.00000000    8.77   0.0031           omega  23   4     1

And see that there are several edges that are included in the model but could improve fit when freed. Let's only look at the first two:

``` r
groupmodel_3 <- groupmodel_2 %>% 
  groupfree("omega","E2","C5") %>% 
  groupfree("omega","E5","A4") %>% 
  runmodel(addMIs = FALSE)
```

    ## Freed 1 parameters!

    ## No parameters need to be freed

    ## Estimating model...

``` r
compare(configural = groupmodel, metric = groupmodel_2, metric_adjusted = groupmodel_3)
```

    ##            model  DF       AIC       BIC  Chisq Chisq_diff DF_diff p_value
    ##       configural 384 125777.90 127469.68 659.50                           
    ##  metric_adjusted 491 125729.41 126848.34 825.01     165.50     107 0.00025
    ##           metric 492 125737.42 126851.00 835.02      10.02       1  0.0016
    ## 
    ## Note: Chi-square difference test assumes models are nested.

This improved the fit of the model, giving evidence that the edges E2 (find it difficult to approach others) - C5 (waste my time) and E5 (take charge) - A4 (love children) differ across genders. Let's try to confirm this in the test data (not computing modification indices to increase speed):

``` r
groupmodel_test_configural <-  ggm(testData, vars = vars, omega = adjacency, groups = "gender") %>% 
  runmodel(addMIs = FALSE)
```

    ## Estimating baseline model...

    ## Estimating model...

``` r
groupmodel_test_metric <- groupmodel_test_configural %>% groupequal("omega") %>% 
  runmodel(addMIs = FALSE)
```

    ## Constrained 108 parameters!
    ## Estimating model...

``` r
groupmodel_test_metric_adjusted <- groupmodel_test_metric  %>% 
  groupfree("omega","E2","C5") %>% 
  groupfree("omega","E5","A4") %>% 
  runmodel(addMIs = FALSE)
```

    ## Freed 1 parameters!

    ## No parameters need to be freed

    ## Estimating model...

``` r
compare(
  configural = groupmodel_test_configural, 
  metric = groupmodel_test_metric, 
  metric_adjusted = groupmodel_test_metric_adjusted)
```

    ##            model  DF      AIC      BIC  Chisq Chisq_diff DF_diff p_value
    ##       configural 384 70490.25 71998.54 713.95                           
    ##  metric_adjusted 491 70428.66 71426.24 866.37     152.42     107  0.0026
    ##           metric 492 70431.37 71424.17 871.08       4.71       1   0.030
    ## 
    ## Note: Chi-square difference test assumes models are nested.

Which provides mixed support, as BIC is lower for the metric model but AIC is lower for the adjusted model.

Latent network modeling (CFA)
=============================

The latent network model takes the following form:

-   **μ** = **τ** + **Λ****μ**<sub>**η**</sub>
-   **Σ** = **Λ****Σ**<sub>**η**</sub>**Λ**<sup>⊤</sup> + **Σ**<sub>**ε**</sub>
-   **Σ**<sub>**η**</sub> = **Δ**<sub>**η**</sub>(**I** − **Ω**<sub>**η**</sub>)<sup>−1</sup>**Δ**<sub>**η**</sub>

with **τ** indicating thresholds, **Λ** factor loadings, **μ**<sub>**η**</sub> latent means, **Ω**<sub>**η**</sub> the latent network, **Δ**<sub>**η**</sub> the scaling of the latent network, and **Σ**<sub>**ε**</sub> the residual variance--covariance structure. These are similarly named in the `lnm()` function.

To showcase the latent network model, I'll take the typical example from the *lavaan* package, which inspired much of the *psychonetrics* package:

``` r
library("lavaan")

## The famous Holzinger and Swineford (1939) example
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)

fit
```

    ## lavaan 0.6-4.1351 ended normally after 35 iterations
    ## 
    ##   Optimization method                           NLMINB
    ##   Number of free parameters                         21
    ## 
    ##   Number of observations                           301
    ## 
    ##   Estimator                                         ML
    ##   Model Fit Test Statistic                      85.306
    ##   Degrees of freedom                                24
    ##   P-value (Chi-square)                           0.000

As a saturated latent network is equivalent to a fully populated variance--covariance structure, we can replicate lavaan exactly:

``` r
Lambda <- matrix(0,9,3)
Lambda[1:3,1] <- Lambda[4:6,2] <- Lambda[7:9,3] <- 1
lnmMod <- lnm(
  HolzingerSwineford1939, 
  vars = paste0("x",1:9), 
  lambda = Lambda, 
  identification = "loadings",
  latents = c("visual", "textual", "speed")
  )
lnmMod <- lnmMod %>% runmodel
```

    ## Estimating baseline model...

    ## Estimating model...

    ## Computing modification indices...

``` r
lnmMod
```

    ##                         _                      _        _          
    ##                        | |                    | |      (_)         
    ##    _ __  ___ _   _  ___| |__   ___  _ __   ___| |_ _ __ _  ___ ___ 
    ##   |  _ \/ __| | | |/ __|  _ \ / _ \|  _ \ / _ \ __|  __| |/ __/ __|
    ##   | |_) \__ \ |_| | (__| | | | (_) | | | |  __/ |_| |  | | (__\__ \
    ##   | .__/|___/\__, |\___|_| |_|\___/|_| |_|\___|\__|_|  |_|\___|___/
    ##   | |         __/ |                                                
    ##   |_|        |___/                                                 
    ##  
    ## 
    ## General: 
    ##  - psychonetrics version: 0.1.2 
    ##  - Model last edited at: 2019-03-05 22:08:18
    ## 
    ## Sample: 
    ##  - Number of cases: 301 
    ##  - Number of groups: 1 
    ##  - Number of observed summary statistics: 54
    ## 
    ## Model: 
    ##  - model used: Latent Network Model (LNM) 
    ##  - Number of parameters: 30
    ## 
    ## Estimation: 
    ##  - Optimizer used: ucminf 
    ##  - Estimator used: Maximum likelihood estimation (ML) 
    ##  - Message: Stopped by small gradient (grtol).
    ## 
    ## Fit: 
    ##  - Model Fit Test Statistic: 85.31 
    ##  - Degrees of freedom: 24 
    ##  - p-value (Chi-square): < 0.0001
    ## 
    ## Tips: 
    ##  - Use 'psychonetrics::compare' to compare psychonetrics models 
    ##  - Use 'psychonetrics::fit' to inspect model fit 
    ##  - Use 'psychonetrics::parameters' to inspect model parameters 
    ##  - Use 'psychonetrics::MIs' to inspect modification indices

We can remove an edge from the latent network, which is different then removing a correlation:

``` r
lnmMod2 <- lnmMod %>% fixpar("omega_eta","speed","textual") %>% runmodel
```

    ## Fixed 1 parameters!

    ## Estimating model...

    ## Computing modification indices...

``` r
compare(lnmMod2, lnmMod)
```

    ##    model DF     AIC     BIC Chisq Chisq_diff DF_diff p_value
    ##  Model 2 24 7535.49 7646.70 85.31                           
    ##  Model 1 25 7534.49 7642.00 86.31       1.00       1    0.32
    ## 
    ## Note: Chi-square difference test assumes models are nested.

Which improves the fit.

Residual network modeling (SEM)
===============================

The residual network model takes the following form:

-   **μ** = **τ**<sub>**y**</sub> + **Λ**(**I** − **B**)<sup>−1</sup>**τ**<sub>**η**</sub>
-   **Σ** = **Λ**(**I** − **B**)<sup>−1</sup>**Σ**<sub>**ζ**</sub>(**I** − **B**)<sup>−1⊤</sup>**Λ**<sup>⊤</sup> + **Σ**<sub>**ε**</sub>
-   **Σ**<sub>**ε**</sub> = **Δ**<sub>**ε**</sub>(**I** − **Ω**<sub>**ε**</sub>)<sup>−1</sup>**Δ**<sub>**ε**</sub>

In which **B** is a matrix of structural effects, **Σ**<sub>**ζ**</sub> a matrix of (residual) latent variance--covariances, **Ω**<sub>**ε**</sub>)<sup>−1</sup> the residual netwokr, and **Δ**<sub>**ε**</sub> the residual scaling. If there is no residual network, the model is equivalent to a structural equation modeling (SEM) model with no residual covariances. Let's again look at lavaan:

``` r
## The industrialization and Political Democracy Example 
## Bollen (1989), page 332
## Adapted to not include residual covariances
model <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
'

fit <- sem(model, data=PoliticalDemocracy)

fit
```

    ## lavaan 0.6-4.1351 ended normally after 40 iterations
    ## 
    ##   Optimization method                           NLMINB
    ##   Number of free parameters                         25
    ##   Number of equality constraints                     3
    ## 
    ##   Number of observations                            75
    ## 
    ##   Estimator                                         ML
    ##   Model Fit Test Statistic                      74.618
    ##   Degrees of freedom                                44
    ##   P-value (Chi-square)                           0.003

We can replicate this analysis in lavaan as follows, making use of integers larger than 1 indicating equality constrains:

``` r
# Lambda with equality constrains:
Lambda <- matrix(0, 11, 3)
Lambda[1:3,1] <- 1
Lambda[4:7,2] <- 2:5
Lambda[8:11,3] <- 2:5

# Beta matrix
beta <- matrix(0,3,3)
beta[2,1] <- beta[3,1] <- beta[3,2] <- 1

vars <- c(paste0("x",1:3),paste0("y",1:8))
latents <- c("ind60","dem60","dem65")

# form RNM model:
rnmMod <- rnm(PoliticalDemocracy, vars = vars, latents = latents, lambda = Lambda, beta = beta)
rnmMod <- rnmMod %>% setoptimizer("ucminf")  %>% runmodel
```

    ## Estimating baseline model...

    ## Estimating model...

    ## Computing modification indices...

``` r
rnmMod
```

    ##                         _                      _        _          
    ##                        | |                    | |      (_)         
    ##    _ __  ___ _   _  ___| |__   ___  _ __   ___| |_ _ __ _  ___ ___ 
    ##   |  _ \/ __| | | |/ __|  _ \ / _ \|  _ \ / _ \ __|  __| |/ __/ __|
    ##   | |_) \__ \ |_| | (__| | | | (_) | | | |  __/ |_| |  | | (__\__ \
    ##   | .__/|___/\__, |\___|_| |_|\___/|_| |_|\___|\__|_|  |_|\___|___/
    ##   | |         __/ |                                                
    ##   |_|        |___/                                                 
    ##  
    ## 
    ## General: 
    ##  - psychonetrics version: 0.1.2 
    ##  - Model last edited at: 2019-03-05 23:09:32
    ## 
    ## Sample: 
    ##  - Number of cases: 75 
    ##  - Number of groups: 1 
    ##  - Number of observed summary statistics: 77
    ## 
    ## Model: 
    ##  - model used: Residual network model (RNM) 
    ##  - Number of parameters: 33
    ## 
    ## Estimation: 
    ##  - Optimizer used: ucminf 
    ##  - Estimator used: Maximum likelihood estimation (ML) 
    ##  - Message: Stopped by small gradient (grtol).
    ## 
    ## Fit: 
    ##  - Model Fit Test Statistic: 74.62 
    ##  - Degrees of freedom: 44 
    ##  - p-value (Chi-square): 0.0027
    ## 
    ## Tips: 
    ##  - Use 'psychonetrics::compare' to compare psychonetrics models 
    ##  - Use 'psychonetrics::fit' to inspect model fit 
    ##  - Use 'psychonetrics::parameters' to inspect model parameters 
    ##  - Use 'psychonetrics::MIs' to inspect modification indices

I made use of the `ucminf` optimizer as it returns the exact same output here as lavaan does (I am not entirely sure why `nlminb` returns a slightly different although comparable *χ*<sup>2</sup> value).

Next, we can find a residual network:

``` r
rnmMod_resid <- rnmMod %>% stepup
```

    ## Estimating model...

    ## Computing modification indices...

    ## Estimating model...

    ## Computing modification indices...

    ## Estimating model...

    ## Computing modification indices...

``` r
compare(rnmMod, rnmMod_resid)
```

    ##    model DF     AIC     BIC Chisq Chisq_diff DF_diff  p_value
    ##  Model 2 41 3173.78 3257.21 44.32                            
    ##  Model 1 44 3198.07 3274.55 74.62      30.29       3 < 0.0001
    ## 
    ## Note: Chi-square difference test assumes models are nested.

The new model fits much better. The residual network can be obtained as follows:

``` r
rnmMod_resid@modelmatrices$`1`$omega_epsilon
```

    ## 11 x 11 sparse Matrix of class "dsCMatrix"
    ##                                                            
    ##  [1,] . . . . .         . .         . .         . .        
    ##  [2,] . . . . .         . .         . .         . .        
    ##  [3,] . . . . .         . .         . .         . .        
    ##  [4,] . . . . .         . .         . .         . .        
    ##  [5,] . . . . .         . 0.3639626 . 0.3874870 . .        
    ##  [6,] . . . . .         . .         . .         . .        
    ##  [7,] . . . . 0.3639626 . .         . .         . .        
    ##  [8,] . . . . .         . .         . .         . .        
    ##  [9,] . . . . 0.3874870 . .         . .         . 0.4095963
    ## [10,] . . . . .         . .         . .         . .        
    ## [11,] . . . . .         . .         . 0.4095963 . .

which has three unique elements. Note that the implied *marginal* covariances are more:

``` r
delta <- rnmMod_resid@modelmatrices$`1`$delta_epsilon
omega <- rnmMod_resid@modelmatrices$`1`$omega_epsilon
delta %*% solve(diag(11) - omega) %*% delta
```

    ## 11 x 11 sparse Matrix of class "dgCMatrix"
    ##                                                                          
    ##  [1,] 0.08049843 .         .         .        .        .        .        
    ##  [2,] .          0.1239368 .         .        .        .        .        
    ##  [3,] .          .         0.4663332 .        .        .        .        
    ##  [4,] .          .         .         1.629821 .        .        .        
    ##  [5,] .          .         .         .        8.479963 .        2.2556198
    ##  [6,] .          .         .         .        .        4.679763 .        
    ##  [7,] .          .         .         .        2.255620 .        3.7120926
    ##  [8,] .          .         .         .        .        .        .        
    ##  [9,] .          .         .         .        3.172186 .        0.8437827
    ## [10,] .          .         .         .        .        .        .        
    ## [11,] .          .         .         .        1.160374 .        0.3086526
    ##                                            
    ##  [1,] .        .         .        .        
    ##  [2,] .        .         .        .        
    ##  [3,] .        .         .        .        
    ##  [4,] .        .         .        .        
    ##  [5,] .        3.1721863 .        1.1603739
    ##  [6,] .        .         .        .        
    ##  [7,] .        0.8437827 .        0.3086526
    ##  [8,] 2.032235 .         .        .        
    ##  [9,] .        5.7060835 .        2.0872640
    ## [10,] .        .         3.614546 .        
    ## [11,] .        2.0872640 .        3.7633262

thus, with only three parameters, *six* residual covariances are added to the model.

References
==========

Chen, Yunxiao, Xiaoou Li, Jingchen Liu, and Zhiliang Ying. 2018. “Robust Measurement via a Fused Latent and Graphical Item Response Theory Model.” *Psychometrika* 83 (3). Springer: 538–62.

Epskamp, Sacha, M.T. Rhemtulla, and Denny Borsboom. 2017. “Generalized Network Psychometrics: Combining Network and Latent Variable Models.” *Psychometrika* 82 (4): 904–27. doi:[10.1007/s11336-017-9557-x](https://doi.org/10.1007/s11336-017-9557-x).

Kan, Kees-Jan, Han LJ van der Maas, and Stephen Z Levine. 2019. “Extending Psychometric Network Analysis: Empirical Evidence Against G in Favor of Mutualism?” *Intelligence* 73. Elsevier: 52–62.
