

###########################################################
##  Functions related to Sk
###########################################################


## Function to take Integral over a expoential function
## Input: lower bound a, upper bound b, log concave function h, derivative at xj: dh
## Output: the integral value

integrate_exp = function(a, b, xj, h, dh) {
  (exp(h(xj) - xj*dh)/dh) * (exp(b*dh) - exp(a*dh))
  }


## Function to generate sample from sk
## Input: Tk vector, intersection z, function h and list of exponential of uk functions
## Output: a new sample sk

sample_sk = function(Tk, z, h, exp_uks) {

  if(length(Tk) != length(z)-1) {
    stop("Number of Tk needs to one less than the number of Z.", call. = FALSE)}

  x_star = NA

  dH = sapply(Tk, function(xj) pracma::fderiv(h, xj))

  nn = length(Tk)
  integrals = integrate_exp(z[1:nn], z[-1], Tk, h, dH)
  norm = sum(integrals)
  pdfs = integrals/norm
  cdfs = cumsum(pdfs)

  U1 = runif(1)
  ind = which(cdfs >= U1)[1]
  interval = c(z[ind], c(z[ind + 1]))

  dh = dH[ind]
  xj = Tk[ind]
  z_min = interval[1]

  while (is.na(x_star) || x_star <= z[1] || x_star >= z[length(z)]) {

    # Inverse CDF Method to get value of x_star
    U2 = runif(1)
    x_star = log(dh * integrals[ind] * U2/exp(h(xj) - xj*dh) + exp(z_min * dh))/dh

  }

  return(x_star)
}
