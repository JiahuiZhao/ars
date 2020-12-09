library(pracma)    # to calculate first derivative using central difference method 


get_sample = function(g = dnorm, bounds = c(-5, 5), n = 100, initial = c(-1, 1)) {
  if (bounds[1] == -Inf) {
    # choose x1 st h'(x_1) = 0
  } 
  
  if (bounds[2] == Inf) {
    # choose h'(x_k) = 0
  }
  
  Tk = initial
  for (k in length(initial):(n - 1)) {
    Tk = sample_k(g, bounds, Tk, initial)
    print(Tk)
  }
  
  
}

sample_k = function(g, bounds, Tk, initial) {
  h = function(x) log(g(x))
  z = calc_z(bounds, Tk, h)
  uks = calc_uk(h, Tk, F)
  exp_uks = calc_uk(h, Tk, T)
  #sk_x = calc_sk(x, z, exp_uks)
  lks = calc_lk(h, Tk)
  
  # Sample x_star from sk 
  # Check x_star
 point = NULL
  while (is.null(point)) {
    x_star = sample(c(runif(5), rnorm(5)), size = 1) # for now 
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