#' Numerical solver for a special case of the Telegrapher PDE
#'
#' @param x0 initial point in space
#' @param tt the time horizon
#' @param N the time resolution (number of time sub-intervals)
#' @param M the space resolution (number of space sub-intervals)
#' @param sp the speed of light
#' @param xScale scale up or down the space interval \eqn{(x_0-ct, x0+ct)}
#' @param f the initial condition as a function of the space coordinate
#'
#' @description {A basic implicit scheme implementation for the formal
#' analog of the Fokker-Planck equation in Minkowski space with one space
#' dimension and sign convention \eqn{-, +}. The initial condition is
#' the dirac delta function and the initial velocity is zero.}
#'
#' @return list of the time-grid, space-grid, and solution matrix.
#' @export pde
pde <- function(x0, tt, N, M, sp = 1, xScale = 1, f = NULL)
{
  b <- sp*tt*xScale # Space-boundary: twice the expected time-drift
  # And then we can use years.
  h <- (2*b)/M
  # Initial Condition: Dirac-Delta function
  if(is.null(f))
  {
    f <- function(y) pdes::indicator(abs(x-x0) <= 10^-8)/h
  }
  g <- function(y) rep(0, length(y)) # Zero initial velocity
  k <- tt/N
  x <- seq(x0-b, x0+b, length.out = M+1)
  tg <- seq(0, tt, length.out = N+1)
  # Coefficients for tridiagonal, under the (-,+) sign convention
  A <- (1+1/(2*k*sp^2)+k/h^2)
  B <- -k/(2*h^2)
  K1 <- (1+1/(k*sp^2))
  K2 <- -1/(2*k*sp^2)
  p <- matrix(0, nrow = N+1, ncol = M+1)
  # Initial condition: position and velocity
  p[1, ] <- f(x)
  p[2, ] <- g(x)*k+p[1, ]
  # BC, already taken care of
  p[, 1] <- rep(f(x[1]), N+1)
  p[, M+1] <- rep(f(x[M+1]), N+1)
  numericInfo <- list(stepsizes = c(k = k, h = h),
                      tridiag = c(B, A, B),
                      coef = c(K1, K2)
  )
  # BC are zero, default
  for(i in 3:(N+1))
  {
    d <- K1*p[i-1, 2:M]+K2*p[i-2, 2:M]
    u <- pdes::const_tridiag_solve(B, A, B, d)
    p[i, 2:M] <- u
  }
  # List for animator
  fk <- list(t = tg, x = x, u = p)
  return(fk)
}
