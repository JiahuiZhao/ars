

###########################################################
##  Functions related to Zk
###########################################################


## Define the function to calculate intersection z given two x abscissae
## Input: corresponding abscissae xj xj1, log concave function h
## Output: a function of zk

calc_zj = function(xj, xj1, h) {
  zj = (h(xj1) - h(xj) - xj1*fderiv(h, xj1) + xj*fderiv(h, xj))/(fderiv(h, xj) - fderiv(h, xj1))
  return(zj)
}


## Function to calculate all tangent intersections
## Input: boundaries, abscissae Tk, logconcave function h
## Output: a vector z, intersections bounded by two boundaries

calc_z = function(bounds, Tk, h) {
  z_Tk = sapply(1:(length(Tk) - 1), function(j) calc_zj(Tk[j], Tk[j + 1], h))
  z = c(bounds[1], z_Tk, bounds[2])
  return(z)
}