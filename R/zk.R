
calc_zj = function(xj, xj1, h) {
  zj = (h(xj1) - h(xj) - xj1*fderiv(h, xj1) + xj*fderiv(h, xj))/(fderiv(h, xj) - fderiv(h, xj1))
  return(zj)
}

calc_z = function(bounds, Tk, h) {
  z_Tk = sapply(1:(length(Tk) - 1), function(j) calc_zj(Tk[j], Tk[j + 1], h)) # can just put in whole vector
  z = c(bounds[1], z_Tk, bounds[2])
  return(z)
}