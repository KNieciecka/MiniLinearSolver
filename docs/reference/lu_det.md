# Matrix determinant calculation using LU decomposition

The function calculates the determinant of a square matrix based on its
LU decomposition.

## Usage

``` r
lu_det(A)
```

## Arguments

- A:

  A square numeric matrix.

## Value

The numeric value of the determinant. If the matrix is singular, it
returns 0.

## Examples

``` r
A <- matrix(c(2, 4, -2, 1, -6, 7, 1, 0, 2), nrow = 3)
lu_det(A)
#> [1] -16
```
