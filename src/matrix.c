#include "matrix.h"

double *matrix_alloc(int nrow, int ncol)
{
    int *mat = (int *)malloc((nrow * ncol) * sizeof(double) + 2 * sizeof(int));
    if (!mat)
        error("Could not allocate memory");
    mat[0] = nrow;
    mat[1] = ncol;
    return ((double *)(mat + 2));
}


void matrix_free(double *mat)
{
    free((int *)mat - 2);
}


int *imatrix_alloc(int nrow, int ncol)
{
    int *mat = (int *)malloc((nrow * ncol) * sizeof(int) + 2 * sizeof(int));
    if (!mat)
        error("Could not allocate memory");
    mat[0] = nrow;
    mat[1] = ncol;
    return ((int *)(mat + 2));
}


void imatrix_free(int *mat)
{
    free((int *)mat - 2);
}


/*
** Compute the simple matrix product of `x` and `y` storing the result in `z`
** This is from src/main/array.c, line 807, in the R source directory.
** Use of this function requires the PKG_LIBS line in the src/Makevars file
*/
void matrix_product(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z)
{
    char *transN = "N";
    double one = 1.0, zero = 0.0;
    F77_CALL(dgemm)(transN, transN, &nrx, &ncy, &ncx, &one,
            x, &nrx, y, &nry, &zero, z, &nrx);
}


/*
** compute y = Ax, the product of matrix A with vector x.
** This is from src/main/array.c, line 800, in the R source directory.
*/
void matrix_vector_product(double *A, int nra, int nca, double *x, double *y)
{
    int ione = 1;
    double one = 1.0, zero = 0.0;
    char *transN = "N";
    F77_CALL(dgemv)(transN, &nra, &nca, &one, A, &nra, x, &ione, &zero, y, &ione);
}


/*
** Approximate the matrix exponential exp(Rt) of matrix R (t is a scalar)
** using the Padé approximation and scaling-and-squaring method implemented
** in the Expokit fortran library function dgpadm.
*/
void matrix_exponential(double *R, double t, int nrow, double *expR)
{
    int ideg = 6,
        m = nrow,
        lwsp = 4*nrow*nrow + ideg + 1,
        ipiv[nrow],
        iexph, ns, iflag;

    double wsp[lwsp], alpha = t;

    F77_CALL(dgpadm)(&ideg, &m, &alpha, R, &m, wsp, &lwsp, ipiv, &iexph, &ns, &iflag);

    memcpy(expR, wsp+iexph-1, nrow*nrow*sizeof(double));
}


SEXP rcl_matrix_exponential(SEXP mat, SEXP t)
{
    int nrow = INTEGER(getAttrib(mat, R_DimSymbol))[0];
    SEXP expR = PROTECT(allocVector(REALSXP, nrow*nrow));
    matrix_exponential(REAL(mat), REAL(t)[0], nrow, REAL(expR));
    SEXP dim = PROTECT(allocVector(INTSXP, 2));
    INTEGER(dim)[0] = INTEGER(dim)[1] = nrow;
    setAttrib(expR, R_DimSymbol, dim);
    UNPROTECT(2);
    return expR;
}
