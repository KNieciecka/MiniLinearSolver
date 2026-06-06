#include "pakiet.h"
#include <stdio.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

SEXP lu_decompose(SEXP rA){

    // zabezpieczenia
    if(!Rf_isReal(rA)) Rf_error("The input matrix must be numeric.");
    if(!Rf_isMatrix(rA)) Rf_error("The input argument must be a matrix.");
    
    int n = Rf_nrows(rA);
    if(n != Rf_ncols(rA)) Rf_error("The input matrix must be square.");
    
    if(n == 0) return Rf_allocVector(VECSXP, 0);

    // alokacja pamieci
    SEXP rL = PROTECT(Rf_allocMatrix(REALSXP, n, n));
    SEXP rU = PROTECT(Rf_allocMatrix(REALSXP, n, n));
    SEXP rP = PROTECT(Rf_allocMatrix(REALSXP, n, n));
    SEXP rSwaps = PROTECT(Rf_allocVector(INTSXP, 1));

    // wejsciowa macierz
    double* A = REAL(rA);
    // dolno trojkatna
    double* L = REAL(rL);
    // gorno trojkatna
    double* U = REAL(rU);
    // permutacje
    double* P = REAL(rP);
    // liczba zamian wierszy
    int* swaps = INTEGER(rSwaps);
    swaps[0] = 0;
    // tymczasowy wektor permutacji
    int* P_tmp = (int*)R_alloc(n, sizeof(int));

    // zapelnienie macierzy i wektora permutacji
    for(int j = 0; j < n; j++) {
        P_tmp[j] = j;
        for(int i = 0; i < n; i++) {
            int idx = j * n + i;
            U[idx] = A[idx];
            L[idx] = 0.0;      
            P[idx] = 0.0; 
        }
    }

    // rozklad LU
    for(int k = 0; k < n; k++) {
        // szukanie elementu glownego - pivot
        double max_val = 0.0;
        int pivot_row = k;
        for(int i = k; i < n; i++) {
            if(fabs(U[k * n + i]) > max_val) {
                max_val = fabs(U[k * n + i]);
                pivot_row = i;
            }
        }

        // gdy macierz wejsciowa jest osobliwa
        if(max_val < 1e-12) {
            swaps[0] = -1; // błąd
            Rf_warning("The input matrix is singular.");
            break; 
        }

        // zamiana wierszy
        if(pivot_row != k) {

            int temp_p = P_tmp[k];
            P_tmp[k] = P_tmp[pivot_row];
            P_tmp[pivot_row] = temp_p;

            for(int j = k; j < n; j++) {
                double temp = U[j * n + k];
                U[j * n + k] = U[j * n + pivot_row];
                U[j * n + pivot_row] = temp;
            }
            
            for(int j = 0; j < k; j++) {
                double temp = L[j * n + k];
                L[j * n + k] = L[j * n + pivot_row];
                L[j * n + pivot_row] = temp;
            }
            swaps[0]++;
        }

        // gauss
        L[k * n + k] = 1.0; // 1 na przekatnej L
        #pragma omp parallel for
        for(int i = k + 1; i < n; i++) {
            double wsp = U[k * n + i] / U[k * n + k];
            L[k * n + i] = wsp;
            
            for(int j = k; j < n; j++) {
                U[j * n + i] -= wsp * U[j * n + k];
            }
        }
    }

    // permutacje
    for(int i = 0; i < n; i++) {
        int org_row = P_tmp[i];
        P[org_row * n + i] = 1.0; 
    }

    // czyszczenie U
    for(int j = 0; j < n; j++) {
        for(int i = j + 1; i < n; i++) {
            U[j * n + i] = 0.0;
        }
    }

    // wyniki
    SEXP list, names;
    PROTECT(list = Rf_allocVector(VECSXP, 4));
    
    SET_VECTOR_ELT(list, 0, rL);
    SET_VECTOR_ELT(list, 1, rU);
    SET_VECTOR_ELT(list, 2, rP);
    SET_VECTOR_ELT(list, 3, rSwaps);

    PROTECT(names = Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, Rf_mkChar("L"));
    SET_STRING_ELT(names, 1, Rf_mkChar("U"));
    SET_STRING_ELT(names, 2, Rf_mkChar("P"));
    SET_STRING_ELT(names, 3, Rf_mkChar("swaps"));
    Rf_setAttrib(list, R_NamesSymbol, names);

    UNPROTECT(6);

    return list;
}

