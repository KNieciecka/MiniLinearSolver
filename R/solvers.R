#' @title LU decomposition of a square matrix
#' 
#' @description The function performs an LU decomposition of a matrix with partial pivoting 
#' using C code with OpenMP.
#' 
#' @param A A square numeric matrix.
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{L} - lower triangular matrix,
#'   \item \code{U} - upper triangular matrix,
#'   \item \code{P} - permutation matrix,
#'   \item \code{swaps} - number of row swaps.
#' }
#' 
#' @examples
#' A <- matrix(c(2, 1, 1, 4, -6, 0, -2, 7, 2), 3, 3, byrow=TRUE)
#' lu_decompose(A)
#' 
#' @useDynLib MiniLinearSolver, .registration = TRUE
#' @export
lu_decompose <- function(A) {
  .Call("lu_decompose", as.matrix(A))
}

#' @title Matrix determinant calculation using LU decomposition
#' 
#' @description The function calculates the determinant of a square matrix based on its LU decomposition.
#' 
#' @param A A square numeric matrix.
#' 
#' @return The numeric value of the determinant. If the matrix is singular, it returns 0.
#' 
#' @examples
#' A <- matrix(c(2, 1, 1, 4, -6, 0, -2, 7, 2), 3, 3, byrow=TRUE)
#' lu_det(A)
#' 
#' @export
lu_det <- function(A) {
  .Call("c_determinant", as.matrix(A))
}

#' @title Solving systems of linear equations using LU decomposition
#' 
#' @description The function solves a system of linear equations $Ax = b$ using 
#' LU decomposition.
#' 
#' @param A A square numeric matrix of the system's coefficients.
#' @param b A numeric vector of the right-hand side constants.
#' 
#' @return A numeric vector representing the solution to the system of equations.
#' 
#' @examples
#' A <- matrix(c(2, 1, 1, 4, -6, 0, -2, 7, 2), 3, 3, byrow=TRUE)
#' b <- c(4, -2, 7)
#' lu_solve(A, b)
#' 
#' @export
lu_solve <- function(A, b) {
  .Call("solve_lu", as.matrix(A), as.numeric(b))
}

#' @title Solving systems of linear equations using SOR method
#' 
#' @description The iterative Successive Over-Relaxation (SOR) method for solving
#' systems of linear equations.
#' 
#' @param A A square numeric matrix of the system's coefficients.
#' @param b A numeric vector of the right-hand side constants.
#' @param omega The relaxation parameter. For \code{omega = 1}, the algorithm reduces 
#' to the Gauss-Seidel method. Default is 1.0.
#' @param tol Error tolerance defining the stopping criterion. Default is 1e-7.
#' @param max_iter Maximum number of iterations. Default is 1000.
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{x} - a vector representing the approximate solution of the system,
#'   \item \code{iterations} - the number of iterations performed.
#' }
#' 
#' @examples
#' A <- matrix(c(10, -1, 2, 0, -1, 11, -1, 3, 2, -1, 10, -1, 0, 3, -1, 8), 4, 4, byrow=TRUE)
#' b <- c(6, 25, -11, 15)
#' sor_solve(A, b, omega = 1.1)
#' 
#' @export
sor_solve <- function(A, b, omega = 1.0, tol = 1e-7, max_iter = 1000L) {
  .Call("solve_SOR", 
        as.matrix(A), 
        as.numeric(b), 
        as.numeric(omega), 
        as.numeric(tol), 
        as.integer(max_iter))
}