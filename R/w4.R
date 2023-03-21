#' @export
w4 <- function(x0,
               func,
               Fa,
               jacfunc = function(x) numDeriv::jacobian(func, x),
               dtau = 0.5,
               maxiter = 1000L,
               errmax = 1e-4,
               decomposition = c("sv", "lu", "lh"),
               trace = FALSE) {

  # If `x0` is a vector
  if (!is.matrix(x0)) {
    w4_(
      x0 = x0,
      func = func,
      Fa = Fa,
      jacfunc = jacfunc,
      dtau = dtau,
      maxiter = maxiter,
      errmax = errmax,
      decomposition = decomposition,
      trace = trace
    )
  } else {
    lst <- purrr::map(
      purrr::array_branch(x0, 2),
      w4_,
      func = func,
      Fa = Fa,
      jacfunc = jacfunc,
      dtau = dtau,
      maxiter = maxiter,
      errmax = errmax,
      decomposition = decomposition,
      trace = trace
    )

    purrr::list_rbind(lst, names_to = "trace_id")

  }
}

w4_ <-
  function(x0,
           func,
           Fa,
           jacfunc = function(x) numDeriv::jacobian(func, x),
           dtau = 0.5,
           maxiter = 1000L,
           errmax = 1e-4,
           decomposition = c("sv", "lu", "lh"),
           trace = FALSE) {

  decomposition <- match.arg(decomposition)

  n <- length(x0)
  x <- x0
  p <- double(length = n)

  if (trace) {
    m <- matrix(data = NA_real_, nrow = maxiter, ncol = n + 2L)
  }

  for (i in seq_len(maxiter)) {

    jac_x <- jacfunc(x)
    f_x <- func(x)
    Fa0 <- Fa(x)
    error <- max_norm(f_x / Fa0)

    # Preconditioners X and Y
    xy <- xy(jac = jac_x, decomposition = decomposition)
    X <- xy$X
    Y <- xy$Y

    # Eq. 30a
    x <- x + dtau * X %*% p
    # Eq. 30b
    p <- (1 - 2 * dtau) * p - dtau * Y %*% f_x

    # Update trace matrix
    if (trace)
      m[i, ] <- c(i, error, as.vector(x))

    if (error < errmax) break
  }

  xf_names <- c("i", "error", paste0("x", seq_len(n)))

  if (trace) {
    res <- m[seq_len(i),]
    colnames(res) <- xf_names
  } else {
    res <- as.list(c(i, error, as.vector(x)))
    names(res) <- xf_names
  }

  return(tibble::as_tibble(res))
}
