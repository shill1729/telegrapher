#' Compute the PDF of the telegraph process under three methods
#'
#' @param x0 initial point
#' @param tt time horizon
#' @param sp speed of light
#' @param ntime number of time sub-intervals in sample-paths
#' @param N number of time sub-intervals in PDE grid
#' @param M number of space sub-intervals in PDE grip
#' @param npaths number of paths in ensemble for Monte-Carlo scheme
#' @param xscale the scale of the region for the PDE.
#'
#' @description {Compute and plot the PDFs and CDFs of the telegraph process
#' at a given time and initial point, using the exact semi-analytic formula, a Monte-Carlo
#' scheme and a numerical PDE solver.}
#' @return list of the three pdfs plus a list of the input.
#' @export telegraph_model
telegraph_model <- function(x0, tt, sp = 1, ntime = 250, N = 400, M = 400, npaths = 50, xscale = 1.5)
{
  callArgs <- data.frame(x0 = x0, tt = tt, sp  = sp, ntime = ntime, N = N, M = M, npaths = npaths, xscale = xscale)
  # TODO: the MC and PDE functions assume x0 is built in to the
  # IC but the exact solution ought to require as input
  # And I want to test different definitions of the Dirac delta
  # so passing one with x0 built in, will have redundancies
  # causing the exact solution to be incorrectly computed.
  # Therefore two identical functions are declared here, ic and ic2
  ic <- function(x) stats::dnorm(x-x0, 0, 0.01) # for PDE and MC
  ic2 <- function(x) stats::dnorm(x, 0, 0.01) # For the exact
  x <- seq(x0-sp*tt*xscale, x0+sp*tt*xscale, length.out = M+1)
  s <- sampleEnsemble(tt, x0, sp, npaths = npaths, ntime = ntime)
  graphics::par(mfrow = c(2, 2))
  graphics::plot(x, ic(x), type = "l", main = "Dirac-Delta")
  plotEnsemble(s)
  # Computations of transition density
  print("Computing solution via Monte-Carlo")
  mc <- monteCarlo(x, ic, s, sp)
  print("Computing solution via PDE scheme")
  pd <- pde(x0, tt, N, M, sp, xscale, ic)
  print("Computing exact solution")
  ex <- dtelegraph(tt, x, x0, ic2, sp)
  # Extract solutions from returned lists
  p_mc <- mc$u[ntime+1, ]
  p_pd <- pd$u[N+1, ]
  p_ex <- ex$u
  # Get space-step size.
  dx <- diff(x)[1]
  # Plotting density and cdf
  lb <- min(p_mc, p_pd, p_ex)
  ub <- max(p_mc, p_pd, p_ex)
  ylim <- c(lb, ub)
  # PDF
  graphics::plot(x, p_mc, type = "l", ylim = ylim, ylab = "p(t, x)", main = "PDF")
  graphics::lines(x, p_pd, col = "red")
  graphics::lines(x, p_ex, col = "blue")
  graphics::abline(v = x0+sp*tt, lty = "dashed")
  graphics::abline(v = x0-sp*tt, lty = "dashed")
  graphics::legend(x = "topleft",
         legend = c("MC", "PDE", "Exact"),
         col = c("black", "red", "blue", NA, NA),
         lty = 1, cex = 0.5)
  # CDF
  graphics::plot(x, cumsum(p_pd)*dx, type = "l", col = "red", ylab = "F(t, x)", main = "CDF")
  graphics::lines(x, cumsum(p_mc)*dx, col = "black")
  graphics::lines(x, cumsum(p_ex)*dx, col = "blue")
  graphics::legend(x = "topleft",
         legend = c("MC", "PDE", "Exact"),
         col = c("black", "red", "blue", NA, NA),
         lty = 1, cex = 0.5)
  # Check normalization holds
  normCheck <- c(pde = sum(p_pd)*dx, mc = sum(p_mc)*dx, ex = sum(p_ex)*dx)
  print(normCheck)
  return(list(x = x, exact = p_ex, monte_carlo = p_mc, pde = p_pd, callArgs = callArgs))
}
