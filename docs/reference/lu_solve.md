# Solving systems of linear equations using LU decomposition

The function solves a system of linear equations `Ax = b` using LU
decomposition.

## Usage

``` r
lu_solve(A, b)
```

## Arguments

- A:

  A square numeric matrix of the system's coefficients.

- b:

  A numeric vector of the right-hand side constants.

## Value

A numeric vector representing the solution to the system of equations.

## Examples

``` r
A <- matrix(c(2, 4, -2, 1, -6, 7, 1, 0, 2), nrow = 3)
b <- c(4, -2, 7)
lu_solve(A, b)
#> [1] 1 1 1
```
