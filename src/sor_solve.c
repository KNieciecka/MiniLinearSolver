#include "pakiet.h"
#include <math.h>

SEXP solve_SOR(SEXP rA, SEXP rb, SEXP rOmega, SEXP rTol, SEXP rMaxIter){

    // zabezbieczenia 
    if(!Rf_isReal(rA) || !Rf_isMatrix(rA)) Rf_error("Macierz wejściowa musi być numeryczna.");
    if(!Rf_isReal(rb)) Rf_error("Wektor wejściowy musi być numeryczny.");
    if(!Rf_isReal(rOmega) || !Rf_isReal(rTol)) Rf_error("Omega i tolerancja muszą być typu numeric.");
    
    int n = Rf_nrows(rA);
    if(n != Rf_ncols(rA)) Rf_error("Macierz wejściowa musi być kwadratowa.");
    if(n != Rf_length(rb)) Rf_error("Długość wektora wejściowego musi odpowiadać wymiarom macierzy wejściowej.");
    if(n == 0) return Rf_allocVector(REALSXP, 0);

    double* A = REAL(rA);
    double* b = REAL(rb);
    double  omega = REAL(rOmega)[0];
    double tol = REAL(rTol)[0];
    int max_iter = Rf_asInteger(rMaxIter);

    for(int i=0; i<n; i++) {
        if(fabs(A[i*n + i]) < 1e-12) {
            Rf_error("Metoda SOR wymaga niezerowych elementów na głównej przekątnej macierzy wejściowej.");
        }
    }
    
    // alokacja
    SEXP rx = PROTECT(Rf_allocVector(REALSXP, n));
    double* x = REAL(rx);

    // przyb pocz
    for( int i=0; i<n; i++){
        x[i] = 0.0;
    }

    // SOR
    int iter;
    double max_diff;

    for(iter = 0; iter<max_iter; iter++){
        max_diff = 0.0;

        for( int i=0; i<n; i++){
            double sum = 0.0;

            for(int j=0; j<n; j++){
                if(i != j){
                    sum += A[j*n + i] * x[j];
                }
            }
            double xi = x[i];
            
            x[i] = (1.0 - omega) * xi + (omega * (b[i] - sum) / A[i*n +i]);
            
            // roznica wynikow
            double diff = fabs(x[i] - xi);
            if(diff > max_diff){
                max_diff = diff;
            }
        }
        if(max_diff < tol){
            iter++;
            break;
        }
    }

    if(iter == max_iter){
        Rf_warning("Osiągnięto maksymalną liczbę iteracji. Metoda może nie być zbierzna.");
    }

    // wyniki
    SEXP rIter = PROTECT(Rf_allocVector(INTSXP, 1));
    INTEGER(rIter)[0] = iter;
    SEXP list, names;
    PROTECT(list = Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, rx);
    SET_VECTOR_ELT(list, 1, rIter);

    PROTECT(names = Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, Rf_mkChar("x"));
    SET_STRING_ELT(names, 1, Rf_mkChar("iteracje"));
    Rf_setAttrib(list, R_NamesSymbol, names);

    UNPROTECT(4);
    
    return list;
}