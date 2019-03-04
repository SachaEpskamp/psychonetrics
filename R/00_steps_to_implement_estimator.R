# These are the steps to implement a new model:

# 1. If needed, adapt 01_classes.R and o3_modelformation_samplestats.R to include different data/summary stats
# 2. Edit fitfunction (07_FIMLestimator_fitfunction.R) and further scripts if needed (07_FIMLestimator_fitfunction_gauss.R)
# 3. Edit 04_generalFit_fitgunction.R to link to new fit function
# 4. Create a new gradient function (07_FIMLestimator_gradient_Gauss.R)
# 5. Edit 04_generalGit_gradient.R to link to this derivative
# 6. Edit an expected hessian if needed
# 7. Edit 04_generalfit_Fisdherinformation to link to the expeced hessian
# 8. Edit b_modelExpansions_addFit.R if needed