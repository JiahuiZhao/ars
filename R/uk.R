## Define uk functions for all abscissae x
## Input: log concave function h, x abscissae Tk
## Output: a list of ukj functions

ukj = function(xj, x, h) {
  h(xj) + (x - xj) * fderiv(h, xj)
}

calc_uk = function(h, Tk, exp = T) {

  uk = function(xj, exp, x) {

    # calculate uk
    if (!exp) {
      fun = function(x) {
        h(xj) + (x - xj) * pracma::fderiv(h, xj)
      }
    } else {

      # calculate exp(uk)
      fun = function(x) {
        exp(h(xj)) * exp((x - xj) * pracma::fderiv(h, xj))
      }
    }

    # return either a list of uk or a list of exp(uk) function
    return(fun)
  }

  uks = lapply(Tk, function(xj) uk(xj, exp))
  return(uks)
}

## Calculate the value of uk for a specific x value
## Input: new sample x, tangent intersection z vector, list of ukj functions
## Output: a numeric value uk(x)

get_uk_x = function(x, z, uks) {
  ind = findInterval(x, z)
  return(uks[[ind]](x))
}
