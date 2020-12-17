
get_uk_x = function(x, z, uks) {
  bool = sapply(1:(length(z) - 1), function(i) x <= z[i + 1] && x >= z[i])
  ind = which(bool == T)
  return(uks[[ind]](x))
}

calc_uk = function(h, Tk, exp = T) {
  uk = function(xj, exp, x) {
    if (!exp) {
      fun = function(x) {
        h(xj) + (x - xj) * fderiv(h, xj)
      }
    } else {
      fun = function(x) {
        exp(h(xj)) * exp((x - xj) * fderiv(h, xj))
      }
    }
    
    return(fun)
  }
  
  uks = lapply(Tk, function(xj) uk(xj, exp))
  return(uks)
}