library(pracma)    # to calculate first derivative using central difference method 


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



get_uk_x = function(x, z, uks) {
  bool = sapply(1:(length(z) - 1), function(i) x <= z[i + 1] && x >= z[i])
  ind = which(bool == T)
  return(uks[[ind]](x))
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


calc_z = function(bounds, Tk, h) {
  z_Tk = sapply(1:(length(Tk) - 1), function(j) calc_zj(Tk[j], Tk[j + 1], h)) # can just put in whole vector
  z = c(bounds[1], z_Tk, bounds[2])
  return(z)
}

calc_zj = function(xj, xj1, h) {
  deriv_xj = fderiv(h, xj)
  deriv_xj1 = fderiv(h, xj1)
  zj = (h(xj1) - h(xj) - xj1*deriv_xj1 + xj*deriv_xj)/(deriv_xj - deriv_xj1)
  return(zj)
}

ukj = function(xj, x, h) {
  h(xj) + (x - xj) * fderiv(h, xj)
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

integrate_exp = function(a, b, xj, h, dh) {
  (exp(h(xj) - xj*dh)/dh) * (exp(b*dh) - exp(a*dh))
}


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
