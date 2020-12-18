

###########################################################
##  Functions related to Lk
###########################################################

## Define lk functions for all abscissae x 
## Input: log concave function h, x abscissae Tk
## Output: a list of lk functions

calc_lkj = function(xj, xj1, h) {
  lk = function(x) {
    ((xj1 - x) * h(xj) + (x - xj) * h(xj1))/(xj1 - xj)
  }
  return(lk)
}

calc_lk = function(h, Tk) {
  lks = lapply(1:(length(Tk) - 1), function(j) calc_lkj(Tk[j], Tk[j + 1], h))
  return(lks)
}


## Calculate the value of lk for a specific x value
## Input: new sample x, tangent intersection z vector, list of lk functions
## Output: a numeric value lk(x)

get_lk_x = function(x, Tk, lks) {
  
  bool = sapply(1:(length(Tk) - 1), function(i) x <= Tk[i + 1] && x >= Tk[i])
  if (any(bool) == T) {
    ind = which(bool == T)
    
    # calculate lk(x)
    return(lks[[ind]](x))
  } else {
    
    # set lk equal to 0 when x < x1 and x > xk
    return(-Inf)
  }
}