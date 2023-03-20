#' @export
w4 <-
  function(x0,
           func,
           Fa,
           jacfunc,
           dt = 0.5,
           maxiter = 1000L,
           errmax = 1e-4,
           decomposition = c("lu", "lh", "sv")
  ) {

  decomposition <- match.arg(decomposition)

  # Initialization
  dim <- length(x0)
  p0 <- rep(0L, dim)
  p <- p0
  x <- x0

  # Main iteration loop
  for (i in seq_len(maxiter)) {

    jac_x <- jacfunc(x)
    f_x <- func(x)
    Fa0 <- Fa(x)

    # Criterion to stop the iteration
    # Evaluation of source func
    err <- 0
    for (k in seq_len(dim))
      err <- max(err, abs(f_x[k] / Fa0[k]))

    if (err < errmax) break

    with(
      xy(jac = jac_x, decomposition = decomposition),
      {
        srcx <- X %*% p
        srcp <- -2 * p - Y %*% f_x

        # Update x and p.
        x <<- x + srcx * dt
        p <<- p + srcp * dt
      }
    )
  }

  return(list(x = x, i = i, err = err, ok = TRUE))
}
