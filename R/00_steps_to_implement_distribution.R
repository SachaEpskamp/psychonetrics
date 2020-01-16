# These are the steps to implement a new distribution:

# 1. Make a script with loglikelihood functions (04_generalfit_loglikelihood_Ising.R)
# 2. Edit 04_generalFit_logLikelihood.R to link to the script made in 1
# 3. Make a script with the fit function (05_MLestimator_fit_Ising.R)
# 4. Edit 05_MLestimator_fitfunction.R to link to the script made in 3
# 5. Make a script for the gradient (05_MLestimator_gradient_Ising.R)
# 6. Edit 04_generalFit_gradient.R to link to the function made in 5
# 7. Make a script for the information function (05_MLestimator_expected_hessian_Ising.R)
# 8. Edit 04_generalfit_FisherInformation.R to link to function made in 7
# 9. If needed, edit 02_algebrahelpers_expectedmodel.R