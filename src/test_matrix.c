#include <R.h>
#include <Rinternals.h>

#include "matrix.h"

SEXP rcl_test_matrix_product(SEXP A, SEXP B);
SEXP rcl_test_matrix_product(SEXP A, SEXP B)
{
    int nrx = INTEGER(getAttrib(A, R_DimSymbol))[0];
    int ncx = INTEGER(getAttrib(A, R_DimSymbol))[1];
    int nry = INTEGER(getAttrib(B, R_DimSymbol))[0];
    int ncy = INTEGER(getAttrib(B, R_DimSymbol))[1];
    SEXP C = PROTECT(allocVector(REALSXP, nrx * ncy));
    SEXP dims = PROTECT(allocVector(INTSXP, 2));
    INTEGER(dims)[0] = nrx;
    INTEGER(dims)[1] = ncy;
    setAttrib(C, R_DimSymbol, dims);
    matrix_product(REAL(A), nrx, ncx, REAL(B), nry, ncy, REAL(C));
    UNPROTECT(2);
    return C;
}


SEXP rcl_test_matrix_vector_product(SEXP A, SEXP x);
SEXP rcl_test_matrix_vector_product(SEXP A, SEXP x)
{
    int nra = INTEGER(getAttrib(A, R_DimSymbol))[0];
    int nca = INTEGER(getAttrib(A, R_DimSymbol))[1];
    if (LENGTH(x) != nca)
        error("incompatible dimensions");
    SEXP y = PROTECT(allocVector(REALSXP, nra));
    matrix_vector_product(REAL(A), nra, nca, REAL(x), REAL(y));
    UNPROTECT(1);
    return y;
}


// repeated matrix product for square matrices
SEXP rcl_test_matrix_cum(SEXP A, SEXP B);
SEXP rcl_test_matrix_cum(SEXP A, SEXP B)
{
    int nrx = INTEGER(getAttrib(A, R_DimSymbol))[0];
    int ncx = INTEGER(getAttrib(A, R_DimSymbol))[1];
    int nry = INTEGER(getAttrib(B, R_DimSymbol))[0];
    int ncy = INTEGER(getAttrib(B, R_DimSymbol))[1];
    if (nrx != nry)
        error("Non-square matrices");
    if (ncx != ncy)
        error("Non-square matrices");
    if (ncx != nrx)
        error("Non-square matrices");
    if (ncy != nry)
        error("Non-square matrices");
    double *z = Calloc(4 * nrx * nrx, double);
    int zstep = nrx * nrx;
    matrix_product(REAL(A), nrx, nrx, REAL(B), nrx, nrx, z);
    for (int i = 1; i < 4; ++i)
        matrix_product(z + (i-1)*zstep, nrx, nrx, REAL(B), nrx, nrx, z + i*zstep);
    SEXP C = PROTECT(allocVector(REALSXP, nrx * nrx));
    SEXP dims = PROTECT(allocVector(INTSXP, 2));
    INTEGER(dims)[0] = nrx;
    INTEGER(dims)[1] = nrx;
    setAttrib(C, R_DimSymbol, dims);
    memcpy(REAL(C), z+3*zstep, nrx * nrx * sizeof(double));
    Free(z);
    UNPROTECT(2);
    return C;
}


SEXP rcl_test_matrix_exponential(SEXP R, SEXP t);
SEXP rcl_test_matrix_exponential(SEXP R, SEXP t)
{
    int nrow = INTEGER(getAttrib(R, R_DimSymbol))[0];
    SEXP expR = PROTECT(allocVector(REALSXP, nrow * nrow));
    matrix_exponential(REAL(R), REAL(t)[0], nrow, REAL(expR));
    SEXP dim = PROTECT(allocVector(INTSXP, 2));
    INTEGER(dim)[0] = nrow;
    INTEGER(dim)[1] = nrow;
    setAttrib(expR, R_DimSymbol, dim);
    UNPROTECT(2);
    return expR;
}
