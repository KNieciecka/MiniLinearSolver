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
A <- matrix(c(2, 1, 1, 4, -6, 0, -2, 7, 2), 3, 3, byrow=TRUE)
lu_det(A)
#> [1] -16
```
