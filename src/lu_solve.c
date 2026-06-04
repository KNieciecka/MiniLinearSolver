#include "pakiet.h"

SEXP solve_lu(SEXP rA, SEXP rb){

    //zabezpieczenia
    if(!Rf_isReal(rA) || !Rf_isMatrix(rA)) Rf_error("Macierz wejściowa musi być numeryczna.");
    if(!Rf_isReal(rb)) Rf_error("Wektor wejściowy musi być numeryczny.");

    int n = Rf_nrows(rA);
    if(n != Rf_ncols(rA)) Rf_error("Macierz wejściowa musi być kwadratowa.");
    if(n != Rf_length(rb)) Rf_error("Długość wektora wejściowego musi odpowiadać wymiarom macierzy wejściowej.");
    if(n ==0) return Rf_allocVector(REALSXP, 0);

    // rozklad LU
    SEXP lu = PROTECT(lu_decompose(rA));
    SEXP rL = VECTOR_ELT(lu, 0);
    SEXP rU = VECTOR_ELT(lu, 1);
    SEXP rP = VECTOR_ELT(lu,2);
    SEXP rSwaps = VECTOR_ELT(lu, 3);
    int swaps = INTEGER(rSwaps)[0];

    if(swaps == -1){
        UNPROTECT(1);
        Rf_error("Macierz wejściowa jest osobliwa, układ nie ma jednoznacznego rozwiązania.");
    }

    double* L = REAL(rL);
    double* U = REAL(rU);
    double* P = REAL(rP);
    double* b = REAL(rb);

    // alokacja
    SEXP rx = PROTECT(Rf_allocVector(REALSXP, n));
    double* x = REAL(rx);
    double* y = (double*)R_alloc(n, sizeof(double));
    double* Pb = (double*)R_alloc(n, sizeof(double));

    // P*b
    for(int i=0; i<n; i++){
        Pb[i]=0.0;
        for(int j=0; j<n; j++){
            Pb[i] += P[j*n + i] * b[j];
        }
    }

    // L*y = Pb
    for(int i=0; i<n; i++){
        double sum = 0.0;
        for(int j=0; j<i; j++){
            sum += L[j*n +i] * y[j];
        }
        y[i] = Pb[i] - sum;
    }

    // U*x = y
    for(int i=n-1; i>=0; i--) {
        double sum = 0.0;
        for(int j=i+1; j<n; j++) {
            sum += U[j*n + i] * x[j];
        }
        x[i] = (y[i] - sum) / U[i*n + i];
    }
    
    UNPROTECT(2); 
    
    return rx;
}