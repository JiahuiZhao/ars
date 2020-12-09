# Adaptive Rejection Sampler
#
#
# Code is developed based on Adaptive Rejection Sampling for Gibbs Sampling
#
#
# Input:
# @param: n: the number of samples need to generate
# @param: f: the function where samples are generated from
# @param: k: initial numbers of abscissaes
# @param: lower: lower bound of D
# @param: upper: upper bound of D
# @param: step: delta x used in finding valid lower and upper bound

# Output: a vector of samples generated, length n

#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

ars <- function(n,
                f, 
                k = 2,
                lower = -Inf, 
                upper = Inf, 
                step = 0.5) {
  
  ##****************************** Check Inputs **********************************##
  
  if(is.numeric(n) == FALSE) {
    stop("Number of samples to generate needs to be a numeric value.", call. = FALSE)}
  if(is.function(dfn) == FALSE) {
    stop("Input funciton should be a function.", call. = FALSE)
  }
  if(lower == upper) {
    stop("Please provide valid lower and upper bounds.", call. = FALSE)
  }
  
  
  # ADD MORE
  
  
  ##*********************** Some Important functions ************************##
  
  ## Check log concave
  # This function checks if a function is log-concave for a given domain. 
  check_concave <- function(f, lower, upper){
    if (lower == -Inf) {lower = .Machine$double.xmin} 
    if (upper == Inf) {upper = .Machine$double.xmax}
    theta <- seq(0,1, 0.01) 
    test <- f(theta*lower+(1-theta)*upper) >= f(lower)^theta*f(upper)^(1-theta)
    return(all(test))
  }
  
  ## Check function well-defined on interval
  
  ## function that draws sample from sk (inverse CDF?)
  
  ## take care of special distributions, e.g. uniform distribution 
  

  
  ##**************************** Initialization *****************************##
  
  # taking log
  h <- function(x, fun = f) {
    return(log(fun(x, ...)))
  }
  
  # taking derivative
  dh <- function(x, h, dx = 1e-8) {
    return((h(x+dx)-h(x))/dx)
  }
  
  zk <- function(x1, x2, h, dh) {
    return((h(x2)-h(x1)-x2*dh(x2)+x1*h(x1))/(dh(x1)-dh(x2)))
  }
  
  uk <- function(x, h, dh) {
    return(list(slope = dh(x), intercept = h(x) - x*dh(x)))
  }
  
  lk <- function(x1, x2, h) {
    return(list(slope = (h(x2)-h(x1))/(x2-x1),
                intercept = (x2*h(x1)-x1*h(x2))/(x2-x1)))
  }
  
  ### set starting values of abscissaes ###
  
  if (lower == -Inf && upper != Inf) {
    
    start = upper
    while(dh(start) <= 0) { start = start - step }
    x0 = c(start, upper)
    
  } else if (l != -Inf && u == Inf) {
    
    start = lower
    while(dh(start) >= 0) { start = start + step }
    x0 = c(lower, start)
    
  } else if (l == -Inf && u == Inf) {
    
    startl = 0 - step
    startu = 0 + step
    while(dh(start1) <= 0) { startl = startl - step }
    while(dh(startu) >= 0) { startu = startu + step }
    x0 = c(startl, startu)
    
  } else {
    x0 = c(lower, upper)
  }
  
  x_abscissae = seq(x0[1], x0[2], length.out = k) # initial abscissae, with k = 2 by default
  
  
  ##*************************** Sampling Loop, Main Function ******************************##
  
  
  samples <- c()
  ite <- 0
  reinitialize <- FALSE
  
  while(length(samples) < n) {
    
    ###  Updating step, when x* in added Tk+1
    if (reinitialize == TRUE) {
      
      x_vec <- sort(c(x, newx), decreasing = FALSE, index.return=TRUE)
      
    }
    
    ### Sampling step
    
    
  }
  
  return(samples)
}
