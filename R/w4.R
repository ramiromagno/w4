#' @importFrom Matrix t
#' @export
w4 <-
  function(x0,
           func,
           Fa,
           jacfunc,
           dt,
           maxiter,
           errmax,
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

    J0 <- jacfunc(x)
    F0 <- func(x)
    Fa0 <- Fa(x)

    # Criterion to stop the iteration
    # Evaluation of source func
    err <- 0
    for (k in seq_len(dim))
      err <- max(err, abs(F0[k] / Fa0[k]))

    if (err < errmax) break

    with(
      xy(jac = J0, decomposition = decomposition),
      {
        srcx <- X %*% p
        srcp <- -2 * p - Y %*% F0

        # Update x and p.
        x <<- x + srcx * dt
        p <<- p + srcp * dt
      }
    )
  }

  return(list(x = x, i = i, err = err, ok = TRUE))
}
