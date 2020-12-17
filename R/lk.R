
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

get_lk_x = function(x, Tk, lks) {
  bool = sapply(1:(length(Tk) - 1), function(i) x <= Tk[i + 1] && x >= Tk[i])
  if (any(bool) == T) {
    ind = which(bool == T)
    return(lks[[ind]](x))
  } else {
    return(-Inf)
  }
}