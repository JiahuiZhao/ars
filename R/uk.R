

###########################################################
##  Functions related to Uk
###########################################################


## Define the function to calculate uk for a specific x value
## Input: new sample x, corresponding abscissae xj, log concave function h
## Output: a function of uk

ukj = function(xj, x, h) { h(xj) + (x - xj) * fderiv(h, xj) }


## Define uk functions for all abscissae x 
## Input: log concave function h, x abscissae Tk
## Output: a list of ukj functions

calc_uk = function(h, Tk, exp = T) {
  
  uk = function(xj, exp, x) {
    
    # calculate uk 
    if (!exp) {
      fun = function(x) {
        h(xj) + (x - xj) * fderiv(h, xj)
      }
    } else {
      
      # calculate exp(uk)
      fun = function(x) {
        exp(h(xj)) * exp((x - xj) * fderiv(h, xj))
      }
    }
    
    # return either a list of uk or a list of exp(uk) function
    return(fun)
  }
  
  uks = lapply(Tk, function(xj) uk(xj, exp))
  return(uks)
}

## Calculate the value of uk for a specific x value
## Input: new sample x, tangent intersection z vector, list of ukj functions
## Output: a numeric value uk(x)

get_uk_x = function(x, z, uks) {
  bool = sapply(1:(length(z) - 1), function(i) x <= z[i + 1] && x >= z[i])
  ind = which(bool == T)
  return(uks[[ind]](x))
}
