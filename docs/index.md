# MiniLinearSolver

[![R-CMD-check](https://github.com/KNieciecka/Pakiet-R-C/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/KNieciecka/Pakiet-R-C/actions/workflows/R-CMD-check.yaml)

**MiniLinearSolver** is a R package made for solving systems of linear
equations. It uses optimized C code with OpenMP API to provide fast and
reliable numerical solutions.

## Features

- **LU Decomposition**: Efficient matrix factorization with partial
  pivoting.
- **Determinant Calculation**: Fast computation based on LU
  decomposition.
- **LU System Solver**: Solve $Ax = b$ systems using LU decomosition.
- **SOR Method**: Successive Over-Relaxation iterative solver for
  diagonally dominant matrices.
- **Multi-threading support**: Multi-threaded execution via OpenMP for
  large-scale matrices.
