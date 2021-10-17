#' Solve the telegrapher PDE via Kac's Monte-Carlo method
#'
#' @param x grid of space coordinates to compute the solution over
#' @param f the initial condition to the wave-equation, as a function of \eqn{(x,...)}
#' @param s a collection of sample-paths of the Kac random time
#' @param sp speed of light
#' @param ... other arguments to pass to the initial condition
#'
#' @description {Monte-Carlo solver of the telegrapher PDE via averaging
#' Goldstein-Kac processes.}
#'
#' @return list of the time grid 't', space grid 'x' and solution matrix 'u'
#' @export monteCarlo
monteCarlo <- function(x, f, s, sp = 1, ...)
{
  n <- nrow(s)-1
  m <- length(x)-1
  x0 <- s[1, 2]
  tau <- (s[, -1]-x0)/sp # Rescale to the randomized time that moves at +/-1 slope
  u <- function(t, x, ...) (f(x+sp*t, ...)+f(x-sp*t, ...))/2
  p <- matrix(0, nrow = n+1, ncol = m+1)
  for(j in 1:(m+1))
  {
    for(i in 1:(n+1))
    {
      p[i, j] <- mean(u(tau[i, ], x[j]))
    }
  }
  fk <- list(t = s[, 1], x = x, u = p)
  return(fk)
}
