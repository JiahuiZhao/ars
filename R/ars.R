

library(pracma)   

## ********************************************  Initialization  **************************************** ##

calc_init_vals <- function(g=dnorm, bounds = c(-Inf, Inf)) {
  
  h <- function(x) {log(g(x))}
  lower <- bounds[1]
  upper <- bounds[2]
  
  # Case 1
  if (lower == -Inf & upper == Inf) {
    x1 = 0
    x2 = 1
    
    if (fderiv(h, x1)>0) {x1=0} else {
      while (fderiv(h,x1) <= 0) {
        x1 <- x1 - 1
      }
    }
    
    if (fderiv(h, x2)<0) {x2=1} else {
      while (fderiv(h,x2) >= 0) {
        x2 <- x2 + 1
      }
    }
    return(c(x1, x2))
  }
  
  # Case 2
  if (lower == -Inf & upper != Inf) {
    x1 = upper - 0.1
    x2 = upper - 0.01
    
    if (fderiv(h, x1) > 0) {x1 = upper - 0.1} else {
      while (fderiv(h, x1)<=0) {
        x1 <- x1 - 1
      }
    }
    return(c(x1, x2))
  }
  
  # Case 3
  if (lower != -Inf & upper == Inf) {
    x1 = lower + 0.01
    x2 = lower + 0.1
    
    if (fderiv(h, x2) < 0) {x2 = lower + 0.1} else {
      while (fderiv(h, x2) >=0) {
        x2 <- x2 + 1
      }
    }
    return(c(x1, x2))
  } 
  
  # Case 4
  if (lower != Inf & upper != Inf){
    mid_point <- mean(c(lower, upper))
    x1 <- mean(c(mid_point, lower))
    x2 <- mean(c(mid_point, upper))
    return(c(x1,x2))
  }
}


## ********************************************  Main Function **************************************** ##

ars = function(g = dnorm, bounds = c(-Inf, Inf), n = 1000, initial = NULL) {
  if (is.null(initial)) {
    initial = calc_init_vals(g, bounds)
  }
  
  new_sample = rep(NA, n)
  Tk = initial
  while(anyNA(new_sample)) {
    # Functions
    h = function(x) log(g(x))
    z = calc_z(bounds, Tk, h)
    uks = calc_uk(h, Tk, F)
    exp_uks = calc_uk(h, Tk, T)
    lks = calc_lk(h, Tk)
    
    # Sample x_star
    x_star = sample_sk(Tk, z, h, exp_uks)
    
    # Check x_star
    w = runif(1)
    uk_x = get_uk_x(x_star, z, uks)
    lk_x = get_lk_x(x_star, Tk, lks)
    
    if (length(which(!is.na(new_sample))) == 0) {
      ind = 1
    } else {
      ind = max(which(!is.na(new_sample))) + 1
    }
    
    update = F
    
    if (w <= exp(lk_x - uk_x)) {
      # accept x_star
      new_sample[ind] = x_star
    } else {
      update = T
      if (w <= exp(h(x_star) - uk_x)) {
        new_sample[ind] = x_star
      } 
    } 
    
    if (update) {
      Tk = sort(c(Tk, x_star))
    }
    
  }
  
  return(new_sample)
  
}
