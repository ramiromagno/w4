xy <- function(jac, decomposition = c("lu", "lh", "sv")) {

  decomposition <- match.arg(decomposition)

  switch (decomposition,
    "lu" = xy_lu(jac),
    "lh" = xy_lh(jac),
    "sv" = xy_sv(jac)
  )
}

#' @importFrom Matrix t
xy_lu <- function(jac) {

  # Inverse matrix of the Jacobian matrix at point x
  jac_inv <- try(solve(jac), silent = TRUE)

  # Check if inversion of Jacobian matrix was possible.
  # If not, return X and Y as NULL.
  if (!is.matrix(jac_inv))
    return(list(X = NULL, Y = NULL))

  with(
    # Perform PLU decomposition
    Matrix::expand(Matrix::lu(t(jac_inv))),
    # Return list with the two preconditioning matrices X and Y.
    list(
      X = t(U),
      Y = t(L) %*% P
    )
  )
}

#' @importFrom Matrix t
xy_lh <- function(jac) {

  QR <- qr(t(jac))
  X <- qr.Q(QR)
  Y <- try(solve(t(qr.R(QR))), silent = TRUE)

  if (!is.matrix(Y))
    return(list(X = NULL, Y = NULL))

  list(
    X = X,
    Y = Y
  )
}

#' @importFrom Matrix t
xy_sv <- function(jac, tol = 1e-7) {

  svd <- svd(jac)
  U <- svd$u

  # Singular values
  s <- svd$d
  s_ <- rep(1 , length(s))
  s_[s > tol] <- 1 / s[s > tol]
  S <- diag(s_)

  # Vt is V transposed.
  Vt <- svd$v

  list(
    X = Vt,
    Y = S %*% t(U)
  )

}
