

integrate_exp = function(a, b, xj, h, dh) {
  (exp(h(xj) - xj*dh)/dh) * (exp(b*dh) - exp(a*dh))
}

sample_sk = function(Tk, z, h, exp_uks) {
  dH = sapply(Tk, function(xj) fderiv(h, xj))
  # Sample which interval the value should belong to, using the densities of Uk
  integrals = sapply(1:length(Tk), function(j) integrate_exp(z[j], z[j + 1], Tk[j], h, dH[j]))
  norm = sum(integrals)
  pdfs = integrals/norm
  cdfs = cumsum(pdfs)
  
  U1 = runif(1)
  ind = which(cdfs >= U1)[1]
  interval = c(z[ind], c(z[ind + 1]))
  
  dh = dH[ind]
  xj = Tk[ind]
  z_min = interval[1]
  
  U2 = runif(1)
  
  x_star = log(dh * integrals[ind] * U2/exp(h(xj) - xj*dh) + exp(z_min * dh))/dh
  
  return(x_star)
}