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
  
  

  cal_zk <- function(x, h, dh) {
    tmp = c()
    for (j in 1:(length(x)-1)) {
      zj = (h(x[j+1])-h(x[j])-x[j+1]*dh(x[j+1])+x[j]*h(x[j]))/(dh(x[j])-dh(x[j+1]))
      tmp = append(tmp, zj)
    }
    return(tmp)
  }
  
  
  
  cal_uk <- function(x, h, dh) {
    m_tmp = c()
    int_tmp = c()
    for (j in 1:length(x)) {
      m_tmp = append(m_tmp, dh(x[j]))
      int_tmp = append(int_tmp, h(x[j] - x[j]*dh(x[j])))
    }
    return(list(slope = m_tmp, intercept = int_tmp))
  }
  
  
  
  cal_lk <- function(x, h) {
    m_tmp = c()
    int_tmp = c()
    for (j in 1:(length(x)-1)) {
      m_tmp = append(m_tmp, (h(x[j+1])-h(x[j]))/(x[j+1]-x[j]))
      int_tmp = append(int_tmp, (x[j+1]*h(x[j])-x[j]*h([j+1]))/(x[j+1]-x[j]))
    }
    return(list(slope = m_tmp,
                intercept = int_tmp))
  }
  
  
  
  ## ************************* Check Boundary & Set Initial abscissae **************** ###
  
  if (lower == -Inf && upper != Inf) { # valid upper, invalid lower
    
    start = upper
    while(dh(start) <= 0) { start = start - step }
    x0 = c(start, upper)
    
  } else if (l != -Inf && u == Inf) { # valid lower, invalid upper
    
    start = lower
    while(dh(start) >= 0) { start = start + step }
    x0 = c(lower, start)
    
  } else if (l == -Inf && u == Inf) { # both lower, upper invalid
    
    startl = 0 - step
    startu = 0 + step
    while(dh(start1) <= 0) { startl = startl - step } 
    while(dh(startu) >= 0) { startu = startu + step }
    x0 = c(startl, startu)
    
  } else {                            # both lower, upper valid
    
    x0 = c(lower, upper)
  }
  
  x_abscissae = seq(x0[1], x0[2], length.out = k+2)[2:(k+1)] # avoid two boundary points
                                                             # use interior points as abscissae
                                                  
  
  ##******************************** Main Function **********************************##
  
  ## initialize variables that will be updated during iteration
  
  samples <- c()
  zk <- cal_zk(x_abscissae, h, dh)
  uk_slope <- cal_uk(x_abscissae, h, dh)[[1]]
  uk_intercept <- cal_uk(x_abscissae, h, dh)[[2]]
  lk_slope <- cal_lk(x_abscissae, h)[[1]]
  lk_intercept <- cal_lk(x_abscissae, h)[[2]]
  update <- FALSE
  
  ### Iteration ###
  
  while(length(samples) < n) {
    
    ### sampling step
    x_star = get_sample()  ## draw sample from sk
    w = runif(1, 0, 1)
    
    
    if (w <= exp()) {                   # match by first condition, l - u
        samples = append(samples, x_star)
      
      } else if (w <= exp(h(x_star) - )) {       # match by second condition, h - u
          samples = append(samples, x_star)
          update = TRUE
      
      }
    
    ###  Updating step, when x* is added (Tk --> Tk+1)
    if (update == TRUE) {
      
      # update x values
      x_abscissae <- sort(c(x_abscissae, x_star))
      
      # update all relevent variables
      zk <- cal_zk(x_abscissae, h, dh)
      uk_slope <- cal_uk(x_abscissae, h, dh)[[1]]
      uk_intercept <- cal_uk(x_abscissae, h, dh)[[2]]
      lk_slope <- cal_lk(x_abscissae, h)[[1]]
      lk_intercept <- cal_lk(x_abscissae, h)[[2]]
      update <- FALSE
      }
    }
  
  return(samples)
}
