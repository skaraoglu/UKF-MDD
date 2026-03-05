# =============================================================================
# models.R
# Coupled oscillator ODE models used as the dynamical prior for UKF parameter
# estimation on fMRI resting-state time series.
#
# All models follow the interface:
#   ode_model(t, x, p)  ->  matrix of derivatives (same dims as x)
# where:
#   t  : scalar dummy time (models have no explicit time dependence)
#   x  : (N_y x N_sigma) matrix of state variables (or (N_y x 1) vector)
#   p  : (N_p x N_sigma) matrix of parameters     (or (N_p x 1) vector)
# =============================================================================

# -----------------------------------------------------------------------------
#' coupled_osc_model_abk
#'
#' Two coupled nonlinear oscillators with separate intrinsic frequencies a, b
#' and a symmetric coupling constant k.  This is the primary model used for
#' fitting pairs of fMRI ROI signals.
#'
#' Equations of motion:
#'   theta1_ddot = -a * sin(theta1) - k * (sin(theta1) - sin(theta2))
#'   theta2_ddot = -b * sin(theta2) + k * (sin(theta1) - sin(theta2))
#'
#' Parameters (p rows):
#'   p[1,] = a  : intrinsic oscillation rate of region 1
#'   p[2,] = b  : intrinsic oscillation rate of region 2
#'   p[3,] = k  : coupling strength between the two regions
#'
#' @param t  Scalar dummy time variable (unused in equations).
#' @param x  State matrix (2 x N_sigma): rows are theta1, theta2.
#' @param p  Parameter matrix (3 x N_sigma): rows are a, b, k.
#' @return   Matrix (2 x N_sigma) of second derivatives.
#' @export
coupled_osc_model_abk <- function(t, x, p) {
  if (is.null(dim(x))) x <- matrix(x, ncol = 1)
  if (is.null(dim(p))) p <- matrix(p, ncol = 1)

  a      <- p[1, ]
  b      <- p[2, ]
  k      <- p[3, ]
  theta1 <- x[1, ]
  theta2 <- x[2, ]

  theta1_ddot <- -a * sin(theta1) - k * (sin(theta1) - sin(theta2))
  theta2_ddot <- -b * sin(theta2) + k * (sin(theta1) - sin(theta2))

  rbind(theta1_ddot, theta2_ddot)
}


# -----------------------------------------------------------------------------
#' coupled_osc_model_gLk
#'
#' Physical pendulum form: gravity g, pendulum length L, and coupling k.
#' Kept for reference / comparison against the simpler abk model.
#'
#' Equations of motion:
#'   theta1_ddot = (-g * sin(theta1) - k*L * (sin(theta1) - sin(theta2))) / L
#'   theta2_ddot = (-g * sin(theta2) + k*L * (sin(theta1) - sin(theta2))) / L
#'
#' @param t  Scalar dummy time variable (unused).
#' @param x  State matrix (2 x N_sigma).
#' @param p  Parameter matrix (3 x N_sigma): rows are g, L, k.
#' @return   Matrix (2 x N_sigma) of second derivatives.
#' @export
coupled_osc_model_gLk <- function(t, x, p) {
  if (is.null(dim(x))) x <- matrix(x, ncol = 1)
  if (is.null(dim(p))) p <- matrix(p, ncol = 1)

  g      <- p[1, ]
  L      <- p[2, ]
  k      <- p[3, ]
  theta1 <- x[1, ]
  theta2 <- x[2, ]

  theta1_ddot <- (-g * sin(theta1) - k * L * (sin(theta1) - sin(theta2))) / L
  theta2_ddot <- (-g * sin(theta2) + k * L * (sin(theta1) - sin(theta2))) / L

  rbind(theta1_ddot, theta2_ddot)
}
