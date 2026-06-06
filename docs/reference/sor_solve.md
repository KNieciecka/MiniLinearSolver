# Solving systems of linear equations using SOR method

The iterative Successive Over-Relaxation (SOR) method for solving
systems of linear equations.

## Usage

``` r
sor_solve(A, b, omega = 1, tol = 1e-07, max_iter = 1000L)
```

## Arguments

- A:

  A square numeric matrix of the system's coefficients.

- b:

  A numeric vector of the right-hand side constants.

- omega:

  The relaxation parameter. For `omega = 1`, the algorithm reduces to
  the Gauss-Seidel method. Default is 1.0.

- tol:

  Error tolerance defining the stopping criterion. Default is 1e-7.

- max_iter:

  Maximum number of iterations. Default is 1000.

## Value

A list containing:

- `x` - a vector representing the approximate solution of the system,

- `iterations` - the number of iterations performed.

## Examples

``` r
A <- matrix(c(10, -1, 2, 0, -1, 11, -1, 3, 2, -1, 10, -1, 0, 3, -1, 8), 4, 4)
b <- c(6, 25, -11, 15)
sor_solve(A, b, omega = 1.1)
#> $x
#> [1]  1  2 -1  1
#> 
#> $iterations
#> [1] 10
#> 
```
