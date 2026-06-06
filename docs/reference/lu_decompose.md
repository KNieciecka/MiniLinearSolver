# LU decomposition of a square matrix

The function performs an LU decomposition of a matrix with partial
pivoting using C code with OpenMP.

## Usage

``` r
lu_decompose(A)
```

## Arguments

- A:

  A square numeric matrix.

## Value

A list containing:

- `L` - lower triangular matrix,

- `U` - upper triangular matrix,

- `P` - permutation matrix,

- `swaps` - number of row swaps.

## Examples

``` r
A <- matrix(c(2, 1, 1, 4, -6, 0, -2, 7, 2), 3, 3, byrow=TRUE)
lu_decompose(A)
#> $L
#>      [,1] [,2] [,3]
#> [1,]  1.0    0    0
#> [2,]  0.5    1    0
#> [3,] -0.5    1    1
#> 
#> $U
#>      [,1] [,2] [,3]
#> [1,]    4   -6    0
#> [2,]    0    4    1
#> [3,]    0    0    1
#> 
#> $P
#>      [,1] [,2] [,3]
#> [1,]    0    1    0
#> [2,]    1    0    0
#> [3,]    0    0    1
#> 
#> $swaps
#> [1] 1
#> 
```
