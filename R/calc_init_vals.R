# This function generates 2 initial abcisscae values required for ASR. 
# This function separates the bounds into 4 different cases.

library(pracma)

calc_init_vals <- function(g=dnorm, bounds = c(-Inf, Inf)) {
  h <- function(x) {log(g(x))}
  lower <- bounds[1]
  upper <- bounds[2]
  
  # Case 1
  if (lower == -Inf & upper == Inf) {
    x1 = 0
    x2 = 1
    
    if (fderiv(g, x1)>0) {x1=0} else {
      while (fderiv(g,x1) <= 0) {
        x1 <- x1 - 1
      }
    }
    
    if (fderiv(g, x2)<0) {x2=1} else {
      while (fderiv(g,x2) >= 0) {
        x2 <- x2 + 1
      }
    }
    return(c(x1, x2))
  }
  
  # Case 2
  if (lower == -Inf & upper != Inf) {
    x1 = upper - 0.1
    x2 = upper - 0.01
    
    if (fderiv(g, x1) > 0) {x1 = upper - 0.1} else {
      while (fderiv(g, x1)<=0) {
        x1 <- x1 - 1
      }
    }
    return(c(x1, x2))
  }
  
  # Case 3
  if (lower != -Inf & upper == Inf) {
    x1 = lower + 0.01
    x2 = lower + 0.1
    
    if (fderiv(g, x2) < 0) {x2 = lower + 0.1} else {
      while (fderiv(g, x2) >=0) {
        x2 <- x2 + 1
      }
    }
    return(c(x1, x2))
  } 

  # Case 4
  if (lower != Inf & upper != Inf){
    mid_point <- 0.5*(upper - lower)
    x1 <- 0.5*(mid_point - lower)
    x2 <- 0.5*(upper - mid_point)
    return(c(x1,x2))
  }
}
