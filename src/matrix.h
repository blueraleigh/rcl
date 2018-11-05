#ifndef RCL_MATRIX_H_
#define RCL_MATRIX_H_

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>  // dgemm

/*
** Functions and macros for creating and accessing matrices
** stored as a flat array in column major order (following R convention).
*/

double *matrix_alloc(int, int);
void matrix_free(double*);
int *imatrix_alloc(int, int);
void imatrix_free(int*);
void matrix_product(double*, int, int, double*, int, int, double*);
void matrix_vector_product(double*, int, int, double*, double*);
void matrix_exponential(double*, double, int, double*);

SEXP rcl_matrix_exponential(SEXP, SEXP);

// DGPADM subroutine from expokit for computing matrix exponential
BLAS_extern void
F77_NAME(dgpadm)(int*, int*, double*, double*, int*,
    double*, int*, int*, int*, int*, int*);

/*

#define mat_raw(m) ((int *)(m) - 2)
#define mat_nrow(m) mat_raw(m)[0]
#define mat_ncol(m) mat_raw(m)[1]

// macro to retrieve the i, j-th item of a matrix
// stored as a double* in column-major order
#define ELEM(m, i, j) (m)[(i) + (j)*mat_nrow(m)]

// macro to retrieve the j-th column of a matrix
// stored as a double* in column-major order
#define COLUMN(m, j) (m) + (j)*mat_nrow(m)

*/

#endif
