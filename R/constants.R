# =============================================================================
# constants.R
# Central definition of all tuning constants used across the UKF pipeline.
# Edit values here; they propagate everywhere automatically.
# =============================================================================

UKF_CONSTANTS <- list(

  # --- Numerical stability (Cholesky / covariance conditioning) -------------
  JITTER_INIT  = 1e-8,   # First jitter added to Pxx/Pyy before Cholesky
  JITTER_MAX   = 1e-2,   # Maximum jitter before falling back to nearPD
  COND_NUM_MAX = 1e12,   # Maximum acceptable condition number for Pyy
  EIGVAL_MIN   = 1e-12,  # Minimum eigenvalue for positive-definiteness check

  # --- Parameter constraints ------------------------------------------------
  PARAM_MIN    = 1e-8,   # Smallest allowed value when forcePositive = TRUE

  # --- Iterative optimisation -----------------------------------------------
  PARAM_TOL_DEFAULT  = 1e-3,   # Default L2 convergence tolerance
  MAXSTEPS_DEFAULT   = 1000,   # Default maximum iterations
  CHISQ_PLATEAU_TOL  = 1e-8,   # Chi-square change below which we declare plateau

  # --- RK4 integration ------------------------------------------------------
  DT_FRACTION        = 0.1,    # dt = DT_FRACTION * dT  (sub-step size)
  STIFFNESS_WARN     = 100     # Warn if any |parameter| exceeds this value
)
