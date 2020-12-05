# Adaptive Rejection Sampler
#
#
# Code is developed based on Adaptive Rejection Sampling for Gibbs Sampling
#
#

# Input:

# Output:

#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

ars <- function(n, dfn, 
                a = -Inf, 
                b = Inf, 
                step, 
                est_mod) {
  
  #### check the inputs
  
  if(is.numeric(n) == FALSE) {
    stop("Number of samples to generate needs to be a numeric value.", call. = FALSE)}
  if(is.function(dfn) == FALSE) {
    stop("Input funciton should be a function.", call. = FALSE)
  }
  if(a == b) {
    stop("Please provide valid lower and upper bounds.", call. = FALSE)
  }
  
}
