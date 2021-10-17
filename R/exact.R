#' The exact solution to the telegraph equation with Dirac IC
#'
#' @param t time input
#' @param x space input
#' @param x0 initial position of particle
#' @param f the initial condition of the PDE, here the Dirac IC
#' @param sp the speed of light
#'
#' @description {See fundamental Green functions for telegrapher PDE.}
#' @return list of time input, space input, and solution matrix
#' @export dtelegraph
dtelegraph <- function(t, x, x0, f, sp = 1)
{
  n <- length(t)
  m <- length(x)
  v <- matrix(0, nrow = n, ncol = m)
  for(i in 1:n)
  {
    for(j in 1:m)
    {

      singularPart <- f(x[j]-(x0-sp*t[i]))+f(x[j]-(x0+sp*t[i]))
      heaveside <- pdes::indicator(sp*t[i]-abs(x[j]-x0) > 0)
      if(heaveside == 1)
      {
        u <- sqrt((sp*t[i])^2-(x[j]-x0)^2)
        i0 <- besselI(sp*u, 0)
        i1 <- besselI(sp*u, 1)
        acPart <- heaveside*(sp*i0+(sp^2)*t[i]*i1/u)
      } else if(heaveside == 0)
      {
        acPart <- 0
      }
      v[i, j] <- 0.5*exp(-t[i]*sp^2)*(singularPart+acPart)
    }
  }
  fk <- list(t = t, x = x, u = v)
  return(fk)
}
