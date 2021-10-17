#' Generate a sample path of the Telegraph photon process
#'
#' @param tt length of time to run a sample-path for
#' @param x0 the initial position of the particle
#' @param sp the speed of light/velocity of the particle
#' @param n the number of time sub-intervals to use
#' @param startSign the starting direction to use +/- 1, can be NULL for random
#'
#' @description {A basic Euler-scheme to generate a sample path
#' of a Telegraph process, or also known as the Goldstein-Kac process.
#' It travels at a constant absolute velocity changing directions randomly,
#' according to a Poisson process.}
#' @return data.frame containing the time step, the Poisson process,
#' the random sign changing process, Kac's "random time" and finally the telegrapher process proper.
#' @export samplePath
samplePath <- function(tt, x0, sp, n, startSign = NULL)
{
  h <- tt/n
  if((sp^2)*h > 1)
  {
    stop("Need smaller step size")
  }
  N <- matrix(0, nrow = n + 1)
  V <- matrix(0, nrow = n + 1)
  tau <- matrix(0, nrow = n + 1)
  if(is.null(startSign))
  {
    B <- stats::rbinom(1, 1, 0.5)
    startSign <- 2*B-1
  }
  N[1] <- 0
  V[1] <- startSign
  tau[1] <- 0
  for(i in 2:(n+1))
  {
    u <- stats::runif(1, 0, 1)
    arrival <- 0
    if(u <= (sp^2)*h)
    {
      arrival <- 1
    }
    N[i] <- N[i-1] + arrival
    V[i] <- V[i-1] - 2*V[i-1]*(N[i]-N[i-1])
    tau[i] <- tau[i-1] + V[i]*h
  }
  return(data.frame(t = (0:n)*h, N = N, V = V, tau = tau, telegraph = x0+sp*tau))
}



#' Plot wrapper for the Telegraph-process
#'
#' @param s data.frame returned from \code{samplePath}
#' @param sp the speed of light
#'
#' @description {This is to plot all of the processes involved
#' in a single sample-path simulation: the Poisson process,
#' the random sign change, and the photon itself.}
#' @return null
#' @export plotTelegraphSample
plotTelegraphSample <- function(s, sp = 1)
{
  n <- nrow(s)-1
  graphics::par(mfrow = c(2, 2))
  plot(s$t, s$N, xlab = "time", ylab = "# of jumps", main = "N_t", type = "s")
  plot(s$t, s$telegraph, xlab = "time", ylab = "x-position", main = "Photon", type = "l")
  graphics::abline(h = 0, lty = "dashed")
  plot(s$t, s$V, xlab = "time", ylab = "random sign", main = "V_t", type = "s")
  plot(stats::density(diff(s$telegraph)), main = "Increment density", type = "l")
  graphics::abline(v = sp*(s$t[n+1]/n), lty = "dashed")
  graphics::abline(v = -sp*(s$t[n+1]/n), lty = "dashed")
}

#' Simulate a sample path of the Goldstein-Kac process
#'
#' @param tt length of time to simulate under
#' @param x0 the initial position of the particle
#' @param sp the speed of light
#' @param npaths number of paths to generate
#' @param ntime number of time sub-intervals
#' @param startSign the starting direction, use NULL for random, otherwise
#' use \eqn{+/- 1}
#'
#' @description {Simulate a sample path of the Goldstein-Kac
#' process.}
#' @return matrix where the first column is the time, the rest
#' are paths
#' @export sampleEnsemble
sampleEnsemble <- function(tt, x0 = 0, sp = 1, npaths = 100, ntime = 1000, startSign = NULL)
{
  h <- tt/ntime
  pathMatrix <- matrix(0, nrow = ntime + 1, ncol = npaths+1)
  pathMatrix[, 1] <- (0:ntime)*h
  for(j in 2:(npaths+1))
  {
    pathMatrix[, j] <- samplePath(tt, x0, sp, ntime, startSign)$telegraph
  }
  return(pathMatrix)
}

#' Wrapper for plotting all sample paths
#'
#' @param ensemble matrix of sample paths returned from \code{samplePath}
#'
#' @description {Plot multiple sample paths at once.}
#' @return NULL
#' @export plotEnsemble
plotEnsemble <- function(ensemble)
{
  npaths <- ncol(ensemble)-1
  plot(ensemble[, 1], ensemble[, 2],
       type = "l",
       ylim = c(min(ensemble[, -1]), max(ensemble[, -1])),
       xlab = "time", ylab = "position", main = "Sample ensemble")
  for(i in 3:(npaths+1))
  {
    graphics::lines(ensemble[, 1], ensemble[, i], lty = "dashed")
  }
}


