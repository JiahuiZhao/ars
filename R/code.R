library(pracma)    # to calculate first derivative using central difference method 


get_sample = function(g = dnorm, bounds = c(-5, 5), n = 100, initial = NULL) {
  if (bounds[1] == -Inf) {
    # choose x1 st h'(x_1) = 0
  } 
  
  if (bounds[2] == Inf) {
    # choose h'(x_k) = 0
  }
  
  if (is.null(initial)) {
    initial = calc_init_vals(g, bounds)
  }
  
  Tk = initial
  for (k in length(initial):(n - 1)) {
    Tk = sample_k(g, bounds, Tk, initial)
    #rint(Tk)
  }
  
}



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

sample_k = function(g, bounds, Tk, initial) {
  h = function(x) log(g(x))
  z = calc_z(bounds, Tk, h)
  uks = calc_uk(h, Tk, F)
  exp_uks = calc_uk(h, Tk, T)
  lks = calc_lk(h, Tk)
  
  # Sample x_star from sk and check x_star
 point = NULL
  while (is.null(point)) {
    x_star = sample_sk(Tk, z, h, exp_uks)
    check_x_star = check_x(x_star, z, Tk, uks, lks, h)
    #print(check_x_star)
    if (!is.null(check_x_star)) {
      point = x_star
    }
  }
 
 return(sort(c(Tk, point)))
}

check_x = function(x, z, Tk, uks, lks, h) {
  # Sample from sk(x) to get x_star
  w = runif(1)
  
  uk_x = get_uk_x(x, z, uks)
  lk_x = get_lk_x(x, Tk, lks)
  
  if (w <= exp(lk_x - uk_x)) {
    # accept x_star
    return(x)
  } else if (w <= exp(h(x) - uk_x)) {
    return(x)
  } else {
    return(NULL)
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
  # Sample from uniform to determine where point lies on domain.
  #z_star = runif(1, min = z[1], max = z[length(z)])
  lens = sapply(1:(length(z) - 1), function(i) z[i + 1] - z[i])
  # How to handle infinite case? 
  # infs = which(is.infinite(lens)) 
  # if (length(infs) > 0) {
  #   probs = lens[-infs]/sum(lens[-infs])
  #   if (length(probs) == 0) {
  #     probs = rep(1/length(infs), length(infs))
  #   } else {
  #     probs = 
  #   }
  # }
  probs = lens/sum(lens)
  z_star = rmultinom(1, 1, probs)
  #bool = sapply(1:length(z), function(j) z_star >= z[j - 1] && z_star <=  z[j]) 
  #ind = which(bool == T)
  ind = which(z_star == 1)
  interval = c(z[ind], c(z[ind + 1]))
  
  dH = sapply(Tk, function(xj) fderiv(h, xj))
  C_3 = sum(sapply(1:length(Tk), function(j) integrate_exp(z[j], z[j + 1], Tk[j], h, dH[j])))
  
  dh = dH[ind]
  xj = Tk[ind]
  z_min = interval[1]
  val = -1
  
  while (val < 0) {
    U = runif(1)
    val = dh * C_3 * U/exp(h(xj) - xj*dh) + exp(z_min * dh)
  }
  
  x_star = log(dh * C_3 * U/exp(h(xj) - xj*dh) + exp(z_min * dh))/dh
  return(x_star)
}

integrate_exp = function(a, b, xj, h, dh) {
  (exp(h(xj) - xj*dh)/dh) * (exp(b*dh) - exp(a*dh))
}

sk_x = function(x) {
  exp_uks[[2]](x)/sum(sapply(1:length(exp_uks), function(i) integrate(exp_uks[[i]], lower = z[i], upper = z[i + 1])$value))
}

sk = function(x, z, exp_uks, ind) {
  exp_uks[[ind]](x)/sum(sapply(1:length(exp_uks), function(i) integrate(exp_uks[[i]], lower = z[i], upper = z[i + 1])$value))
}


calc_sk = function(x, z, exp_uks) {
  bool = sapply(1:(length(z) - 1), function(i) x <= z[i + 1] && x >= z[i])
  ind = which(bool == T)
  return(sk(x, z, exp_uks, ind))
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