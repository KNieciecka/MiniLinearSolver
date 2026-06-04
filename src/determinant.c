#include "pakiet.h"

SEXP c_determinant(SEXP rA){

    SEXP lu = PROTECT(lu_decompose(rA));

    if(Rf_length(lu) == 0){
        UNPROTECT(1);
        return Rf_allocVector(REALSXP, 0);
    }

    SEXP rU =VECTOR_ELT(lu, 1);
    SEXP rSwaps = VECTOR_ELT(lu, 3);

    double* U = REAL(rU);
    int swaps = INTEGER(rSwaps)[0];

    int n = Rf_nrows(rU);
    SEXP rDet = PROTECT(Rf_allocVector(REALSXP, 1));

    if (swaps == -1){
        REAL(rDet)[0] = 0.0;
    } else {
        double det = 1.0;
        for(int i=0; i<n; i++){
            det *= U[i*n+i];
        }
        if(swaps % 2 !=0){
            det = -det;
        }
        REAL(rDet)[0] = det;
    }

    UNPROTECT(2);

    return rDet;
}