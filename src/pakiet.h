#ifndef PAKIET_PAKIET
#define PAKIET_PAKIET

#include <R.h>
#include <Rinternals.h>

SEXP lu_decompose(SEXP rA);
SEXP c_determinant(SEXP rA);  
SEXP solve_lu(SEXP rA, SEXP rb);
SEXP solve_SOR(SEXP rA, SEXP rb, SEXP rOmega, SEXP rTol, SEXP rMaxIter);

#endif