#include <complex.h>
/****************************************************************
PAL_csr_matvec computes y = A*x, where A is stored in CSR format.
A is real, and x and y are double complex arrays.
****************************************************************/
void PAL_csr_matvec( int n, int *ptrrow, int *indrow, double *A, 
                     double complex *x, double complex *y ) 
{
    int i, j, k;
    double complex zero={0.0+0.0*I};

    for ( i=0; i<n; i++ ) {
        y[i] = zero;
    }
    for ( i=0; i<n; i++ ) { 
        for ( j=ptrrow[i]; j<ptrrow[i+1]; j++ ) {
            k = indrow[j];
            y[k] += A[j]*x[i];
        }
    }

    return;
}
