# These are the steps to implement a new model:

# 1. Make a model function, such as a_models_rnm.R (copy and adapt from previous model)
# 2. Make a preperation function (12_rnm_prepare.R), this will require at some point 3
# 3. Make a function that computes the implied model structures (12_rnm_implied.R)
# 3b. Edit 04_generalfit_impliedModel.R to link to the new implied function
# 4. Edit 04_generalit_prepareModel.R to link to the new model preperation function
# 5. If identification is needed, create an identification function (e.g., 12_rnm_identify) and edit b_modelexpansions_identify to link to it
# 5b. If identification is not needed, edit b_modelexpansions_identify to not result in an error
# 6. Create a script with the derivatives funtions (12_rnm_derivatives.R)
# 7. Edit 04_generalFit_gradient.R to link to new derivative function
# 8. Also edit 04_generalfit_FisherInformation to link the derivative function
# 9. Edit a nicer output in f_convenience_printMethod.R
# 10. Edit d_stepup.R to add a default
# 11. Edit e_modelmodifications_prune to add a default
# 12. Edit h_modelsearch.R to add a default