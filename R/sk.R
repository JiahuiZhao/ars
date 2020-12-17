

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
  
  dH = sapply(Tk, function(xj) fderiv(h, xj))
  
  # Decide which interval the value should belong to, using the densities of Uk
  integrals = sapply(1:length(Tk), 
                     function(j) integrate_exp(z[j], z[j + 1], Tk[j], h, dH[j]))
  norm = sum(integrals)
  pdfs = integrals/norm
  cdfs = cumsum(pdfs)
  
  U1 = runif(1)
  ind = which(cdfs >= U1)[1]
  interval = c(z[ind], c(z[ind + 1]))
  
  dh = dH[ind]
  xj = Tk[ind]
  z_min = interval[1]
  
  # Inverse CDF Method to get value of x_star
  U2 = runif(1)
  
  x_star = log(dh * integrals[ind] * U2/exp(h(xj) - xj*dh) + exp(z_min * dh))/dh
  
  return(x_star)
}